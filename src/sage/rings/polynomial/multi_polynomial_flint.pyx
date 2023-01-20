"""
Sparse multivariate polynomials over `\ZZ`, implemented using FLINT
"""

import sys
import os
import re

import itertools

import threading

from collections import Counter

from sage.rings.integer_ring import ZZ

from sage.rings.polynomial.polydict cimport ETuple

from sage.rings.fraction_field_element import is_FractionFieldElement

from sage.structure.element import coerce_binop
from sage.structure.element cimport Element, CommutativeRingElement
from sage.structure.richcmp import rich_to_bool
from sage.structure.category_object cimport CategoryObject

from sage.structure.factorization import Factorization

from sage.libs.gmp.mpz cimport mpz_get_ui, mpz_get_si

from sage.libs.flint.mpoly cimport *
from sage.libs.flint.fmpz cimport *
from sage.libs.flint.fmpz_mpoly cimport *
from sage.libs.flint.fmpz_mpoly_factor cimport *

from libc.stdlib cimport malloc, realloc, free
from libc.stdint cimport UINT64_MAX
from libc.signal cimport raise_, SIGSEGV
from libc.errno cimport errno, EAGAIN
from posix.stdio cimport FILE, popen, fileno, pclose
from posix.strings cimport bcopy
from posix.fcntl cimport creat, open, O_RDONLY
from posix.unistd cimport read, write, close, dup
from posix.stat cimport S_IRWXU, S_IRWXG, S_IRWXO, S_IRUSR, S_IWUSR, S_IRGRP, S_IWGRP, S_IROTH, S_IWOTH
from posix.wait cimport WIFEXITED, WEXITSTATUS

from sage.misc.sage_eval import sage_eval
from sage.misc.misc_c import prod

from sage.functions.other import binomial

from sage.cpython.string cimport char_to_str, str_to_bytes

from functools import reduce

# Some system include files that Cython should probably provide, but doesn't

cdef extern from "<semaphore.h>" nogil:
    ctypedef union sem_t:
        pass
    int sem_init(sem_t *sem, int pshared, unsigned int value)
    int sem_wait(sem_t *sem)
    int sem_trywait(sem_t *sem)
    int sem_post(sem_t *sem)
    int sem_destroy(sem_t *sem)
    int sem_getvalue(sem_t *sem, int *sval)

from posix.time cimport timeval, timespec
from posix.types cimport clockid_t

cdef extern from "<pthread.h>" nogil:
    ctypedef union pthread_mutex_t:
        pass
    ctypedef union pthread_mutexattr_t:
        pass
    ctypedef union pthread_cond_t:
        pass
    ctypedef union pthread_condattr_t:
        pass
    ctypedef union pthread_rwlock_t:
        pass
    ctypedef union pthread_rwlockattr_t:
        pass

    ctypedef enum pthread_rwlock_prefer:
        PTHREAD_RWLOCK_PREFER_READER_NP,
        PTHREAD_RWLOCK_PREFER_WRITER_NP,
        PTHREAD_RWLOCK_PREFER_WRITER_NONRECURSIVE_NP,
        PTHREAD_RWLOCK_DEFAULT_NP = PTHREAD_RWLOCK_PREFER_READER_NP

    int pthread_mutex_init (pthread_mutex_t *__mutex, const pthread_mutexattr_t *__mutexattr)
    int pthread_mutex_destroy (pthread_mutex_t *__mutex)
    int pthread_mutex_trylock (pthread_mutex_t *__mutex)
    int pthread_mutex_lock (pthread_mutex_t *__mutex)
    int pthread_mutex_timedlock (pthread_mutex_t * __mutex, const timespec * __abstime)
    int pthread_mutex_clocklock (pthread_mutex_t * __mutex, clockid_t __clockid, const timespec * __abstime)
    int pthread_mutex_unlock (pthread_mutex_t *__mutex)
    int pthread_mutex_getprioceiling (const pthread_mutex_t * __mutex, int * __prioceiling)
    int pthread_mutex_setprioceiling (pthread_mutex_t * __mutex, int __prioceiling, int * __old_ceiling)
    int pthread_mutex_consistent (pthread_mutex_t *__mutex)
    int pthread_mutex_consistent_np (pthread_mutex_t *__mutex)

    int pthread_mutexattr_init (pthread_mutexattr_t *__attr)
    int pthread_mutexattr_destroy (pthread_mutexattr_t *__attr)
    int pthread_mutexattr_getpshared (const pthread_mutexattr_t * __attr, int * __pshared)
    int pthread_mutexattr_setpshared (pthread_mutexattr_t *__attr, int __pshared)
    int pthread_mutexattr_gettype (const pthread_mutexattr_t * __attr, int * __kind)
    int pthread_mutexattr_settype (pthread_mutexattr_t *__attr, int __kind)
    int pthread_mutexattr_getprotocol (const pthread_mutexattr_t * __attr, int * __protocol)
    int pthread_mutexattr_setprotocol (pthread_mutexattr_t *__attr, int __protocol)
    int pthread_mutexattr_getprioceiling (const pthread_mutexattr_t * __attr, int * __prioceiling)
    int pthread_mutexattr_setprioceiling (pthread_mutexattr_t *__attr, int __prioceiling)
    int pthread_mutexattr_getrobust (const pthread_mutexattr_t *__attr, int *__robustness)
    int pthread_mutexattr_getrobust_np (const pthread_mutexattr_t *__attr, int *__robustness)
    int pthread_mutexattr_setrobust (pthread_mutexattr_t *__attr, int __robustness)
    int pthread_mutexattr_setrobust_np (pthread_mutexattr_t *__attr, int __robustness)

    int pthread_cond_init (pthread_cond_t * __cond, const pthread_condattr_t * __cond_attr)
    int pthread_cond_destroy (pthread_cond_t *__cond)
    int pthread_cond_signal (pthread_cond_t *__cond)
    int pthread_cond_broadcast (pthread_cond_t *__cond)
    int pthread_cond_wait (pthread_cond_t * __cond, pthread_mutex_t * __mutex)
    int pthread_cond_timedwait (pthread_cond_t * __cond, pthread_mutex_t * __mutex, const timespec * __abstime)
    int pthread_cond_clockwait (pthread_cond_t * __cond, pthread_mutex_t * __mutex,
                                clockid_t __clock_id, const timespec * __abstime)
    int pthread_condattr_init (pthread_condattr_t *__attr)
    int pthread_condattr_destroy (pthread_condattr_t *__attr)
    int pthread_condattr_getpshared (const pthread_condattr_t * __attr, int * __pshared)
    int pthread_condattr_setpshared (pthread_condattr_t *__attr, int __pshared)
    int pthread_condattr_getclock (const pthread_condattr_t * __attr, clockid_t * __clock_id)
    int pthread_condattr_setclock (pthread_condattr_t *__attr, clockid_t __clock_id)

    int pthread_rwlock_init (pthread_rwlock_t * __rwlock, const pthread_rwlockattr_t * __attr)
    int pthread_rwlock_destroy (pthread_rwlock_t *__rwlock)
    int pthread_rwlock_rdlock (pthread_rwlock_t *__rwlock)
    int pthread_rwlock_tryrdlock (pthread_rwlock_t *__rwlock)
    int pthread_rwlock_timedrdlock (pthread_rwlock_t * __rwlock, const timespec * __abstime)
    int pthread_rwlock_clockrdlock (pthread_rwlock_t * __rwlock, clockid_t __clockid, const timespec * __abstime)
    int pthread_rwlock_wrlock (pthread_rwlock_t *__rwlock)
    int pthread_rwlock_trywrlock (pthread_rwlock_t *__rwlock)
    int pthread_rwlock_timedwrlock (pthread_rwlock_t * __rwlock, const timespec * __abstime)
    int pthread_rwlock_clockwrlock (pthread_rwlock_t * __rwlock, clockid_t __clockid, const timespec * __abstime)
    int pthread_rwlock_unlock (pthread_rwlock_t *__rwlock)

    int pthread_rwlockattr_init (pthread_rwlockattr_t *__attr)
    int pthread_rwlockattr_destroy (pthread_rwlockattr_t *__attr)
    int pthread_rwlockattr_getpshared (const pthread_rwlockattr_t * __attr, int * __pshared)
    int pthread_rwlockattr_setpshared (pthread_rwlockattr_t *__attr, int __pshared)
    int pthread_rwlockattr_getkind_np (const pthread_rwlockattr_t * __attr, int * __pref)
    int pthread_rwlockattr_setkind_np (pthread_rwlockattr_t *__attr, int __pref)


status_string = ""
status_string_encode = status_string.encode()
cdef char * status_string_ptr = status_string_encode

# Raising a Python exception in a Cython callback from FLINT does nothing other than print a message,
# so deal with fatal errors by generating a seg fault, which will be caught by gdb if running with "sage --gdb"

# Functions to compute the coefficents used to encode deglex exponents
#
# We use a C lookup table for speed.  The deglex coefficient (based on binomial coefficients) is computed
# using Sage bigints, to avoid any integer overflow issues, then converted to a 64-bit ulong
# to be stored in the lookup table.  Some care is still required if the end result doesn't
# fit into 64 bits; in this case the largest 64-bit integer is stored; the table entry won't be correct,
# but this is detected and handled.

cdef ulong choose_with_replacement(ulong setsize, ulong num):
    cdef ulong retval
    value = binomial(setsize + num - 1, num)
    try:
        retval = <ulong> value
    except OverflowError:
        retval = UINT64_MAX
    return retval

# deglex(len_exps, 0, -1) = 0 will be called for constant terms; otherwise offset is non-zero
# and has to be >= num-1

cpdef ulong deglex_coeff_slow(ulong len_exps, ulong num, ulong offset):
    cdef ulong i
    cdef ulong choose
    cdef ulong retval = 0
    for i in range(0, num):
        choose = choose_with_replacement(len_exps, offset-i)
        if choose == UINT64_MAX: raise_(SIGSEGV)
        if (UINT64_MAX - choose < retval): raise_(SIGSEGV)
        retval += choose
    return retval

ctypedef struct deglex_table_entry:
    ulong * table
    ulong size

ctypedef struct deglex_table_entry2:
    deglex_table_entry * table
    ulong size

cdef deglex_table_entry2 * deglex_table = NULL
cdef ulong deglex_table_size = 0

# Create a rwlock to protect the deglex coefficient table.
# Set its locking policy to block readers whenever a writer is waiting;
# i.e, we stop everything whenever a thread wants to extend the table

cdef pthread_rwlock_t deglex_table_rwlock
cdef pthread_rwlockattr_t deglex_table_rwlockattr

pthread_rwlockattr_init(& deglex_table_rwlockattr)
# pthread_rwlockattr_setkind_np(& deglex_table_rwlockattr, PTHREAD_RWLOCK_PREFER_WRITER_NONRECURSIVE_NP)
pthread_rwlock_init(& deglex_table_rwlock, & deglex_table_rwlockattr)

cdef int require_deglex_table_prefilled = 1

cpdef void deglex_fill_table(ulong setsize, ulong num, ulong offset) nogil:
    global deglex_table, deglex_table_size, require_deglex_table_prefilled
    cdef ulong size
    if require_deglex_table_prefilled or pthread_rwlock_wrlock(& deglex_table_rwlock) != 0:
        raise_(SIGSEGV)
    with gil:
        if setsize >= deglex_table_size:
            deglex_table = <deglex_table_entry2 *> realloc(deglex_table,
                                                              (setsize+1) * sizeof(deglex_table_entry2))
            while deglex_table_size <= setsize:
                deglex_table[deglex_table_size].table = NULL
                deglex_table[deglex_table_size].size = 0
                deglex_table_size += 1

        if num >= deglex_table[setsize].size:
            deglex_table[setsize].table = <deglex_table_entry *> realloc(deglex_table[setsize].table,
                                                             (num+1) * sizeof(deglex_table_entry))
            while deglex_table[setsize].size <= num:
                size = deglex_table[setsize].size
                deglex_table[setsize].table[size].table = NULL
                deglex_table[setsize].table[size].size = 0
                deglex_table[setsize].size += 1

        if offset - (num-1) >= deglex_table[setsize].table[num].size:
            deglex_table[setsize].table[num].table = <ulong *> realloc(deglex_table[setsize].table[num].table,
                                                             (offset - (num-1) + 1) * sizeof(ulong))
            while deglex_table[setsize].table[num].size <= offset - (num-1):
                size = deglex_table[setsize].table[num].size
                deglex_table[setsize].table[num].table[size] = deglex_coeff_slow(setsize, num, size + (num-1))
                deglex_table[setsize].table[num].size += 1
    if pthread_rwlock_unlock(& deglex_table_rwlock) != 0:
        raise_(SIGSEGV)

def deglex_prefill_table(num_exps, max_degree):
    global require_deglex_table_prefilled
    save_require_deglex_table_prefilled = require_deglex_table_prefilled
    require_deglex_table_prefilled = 0
    for i in range(num_exps+1):
        for j in range(max_degree+1):
            deglex_fill_table(i, j, max_degree)
    require_deglex_table_prefilled = save_require_deglex_table_prefilled

cpdef ulong deglex_coeff(ulong setsize, ulong num, ulong offset) nogil:
    """
        TESTS::
            sage: from sage.rings.polynomial.multi_polynomial_flint import deglex_coeff
            sage: deglex_coeff(60, 0, 3)
            0
    """
    global deglex_table, deglex_table_size, require_deglex_table_prefilled
    cdef ulong retval

    if not require_deglex_table_prefilled:
        if pthread_rwlock_rdlock(& deglex_table_rwlock) != 0:
            raise_(SIGSEGV)

    if setsize >= deglex_table_size or num >= deglex_table[setsize].size \
       or offset - (num-1) >= deglex_table[setsize].table[num].size:
        if require_deglex_table_prefilled:
            raise_(SIGSEGV)
        if pthread_rwlock_unlock(& deglex_table_rwlock) != 0:
            raise_(SIGSEGV)
        deglex_fill_table(setsize, num, offset)
        if pthread_rwlock_rdlock(& deglex_table_rwlock) != 0:
            raise_(SIGSEGV)

    retval = deglex_table[setsize].table[num].table[offset - (num-1)]
    if not require_deglex_table_prefilled:
        if pthread_rwlock_unlock(& deglex_table_rwlock) != 0:
            raise_(SIGSEGV)
    return retval

# Functions to encode and decode deglex exponents

# encoding raises SIGSEGV in an overflow situation

cdef ulong encode_deglex(ulong * exps, ulong len_exps, int rev=0) nogil:
    cdef ulong delta = 0
    cdef ulong i
    cdef ulong exp
    for i in range(len_exps): delta += exps[i]
    cdef ulong retval = deglex_coeff(len_exps, delta, delta-1)
    for i in range(0,len_exps-1):
        exp = exps[i] if not rev else exps[len_exps - i - 1]
        retval += deglex_coeff(len_exps-i-1, exp, delta)
        delta -= exp
    # no need to encode the last exponent, since we encoded the total degree first
    return retval

def encode_deglex_test(exps):
    cdef ulong len_exps = len(exps)
    cdef ulong * cexps = <ulong *> malloc(len_exps)
    for i in range(len_exps): cexps[i] = exps[i]
    cdef retval = encode_deglex(cexps, len_exps)
    free(cexps)
    return retval

# decoding never raises SIGSEGV because a number between 0 and UINT64_MAX
# will always decode to some exponent vector
#
# Maybe these loops could be replaced with binary search for speed

cdef void decode_deglex(ulong ind, ulong * exps, ulong len_exps, int rev=0) nogil:
    cdef ulong total_degree = 0
    cdef ulong ind_saved = ind

    while True:
        if ind < deglex_coeff(len_exps, total_degree+1, total_degree): break
        total_degree += 1
    ind -= deglex_coeff(len_exps, total_degree, total_degree-1)

    cdef ulong d = total_degree
    cdef unsigned char this_exp
    cdef ulong i
    for i in range(0, len_exps-1):
        this_exp = 0
        while True:
            if ind < deglex_coeff(len_exps-i-1, this_exp+1, d): break
            this_exp += 1
        exps[i if not rev else len_exps-i-1] = this_exp
        if this_exp > 0:
            ind -= deglex_coeff(len_exps-i-1, this_exp, d)
            d -= this_exp
    exps[len_exps-1 if not rev else 0] = d
    #if encode_deglex(exps, len_exps) != ind_saved:
    #    raise_(SIGSEGV)

def decode_deglex_test(ind, len_exps):
    """
    TESTS::
        sage: from sage.rings.polynomial.multi_polynomial_flint import decode_deglex_test
        sage: decode_deglex_test(517, 4)
        [0, 2, 3, 4]
        sage: decode_deglex_test(168, 4)
        [1, 2, 3, 0]
        sage: decode_deglex_test(803, 4)
        [1, 2, 3, 4]
    """
    cdef ulong * cexps = <ulong *> malloc(len_exps*sizeof(ulong))
    decode_deglex(ind, cexps, len_exps)
    retval = []
    for i in range(len_exps): retval.append(cexps[i])
    free(cexps)
    return retval

# Encode/decode to/from a pointer to an address in memory

cdef void encode_to_mem(encoding_format * format, ulong * dest, flint_bitcnt_t bits, const ulong * exp, const fmpz_t coeff, const fmpz_mpoly_ctx_t ctx) nogil:
    cdef int nvarsencoded = 0
    cdef int i
    cdef ulong user_exps[256]
    if ctx.minfo.nvars > 256:
        # NotImplementedError, but we can't raise Python exceptions in a callback function
        # should have some way to dynamically size user_exps, but this is easier for now
        raise_(SIGSEGV)
    if fmpz_is_mpz(coeff):
        # Can't currently encode bigints
        raise_(SIGSEGV)
    mpoly_get_monomial_ui_sp(user_exps, exp, bits, ctx.minfo)
    for i in range(format.words):
        if format.variables[i] > 0:
            # ctx.minfo.rev indicates if our ordering is reversed in the mathematical sense (degrevlex)
            # Since we list our variables from MSV to LSV, but our byte ordering is LSB to MSB, FLINT
            # reversed the variables in the array if "not ctx.minfo.rev", so we need to reverse the
            # order than we encode and decode if "not ctx.minfo.rev".
            if ctx.minfo.rev:
                dest[i] = encode_deglex(user_exps + nvarsencoded, format.variables[i])
            else:
                dest[i] = encode_deglex(user_exps + ctx.minfo.nvars - nvarsencoded - format.variables[i], format.variables[i], rev=1)
            nvarsencoded += format.variables[i]
        else:
            dest[i] = (<ulong *>coeff)[0]

cdef void decode_from_mem(encoding_format * format, const ulong * src, flint_bitcnt_t bits, ulong * exp, fmpz_t coeff, const fmpz_mpoly_ctx_t ctx) nogil:
    cdef slong N = mpoly_words_per_exp(bits, ctx.minfo)
    cdef int i
    cdef int nvarsencoded = 0
    cdef ulong user_exps[256]
    if ctx.minfo.nvars > 256:
        # NotImplementedError, but we can't raise Python exceptions in a callback function
        # should have some way to dynamically size user_exps, but this is easier for now
        raise_(SIGSEGV)
    for i in range(format.words):
        if format.variables[i] > 0:
            # see comment in encode_to_mem
            if ctx.minfo.rev:
                decode_deglex(src[i], user_exps + nvarsencoded, format.variables[i])
            else:
                decode_deglex(src[i], user_exps + ctx.minfo.nvars - nvarsencoded - format.variables[i], format.variables[i], rev=1)
            nvarsencoded += format.variables[i]
        else:
            fmpz_set_si(coeff, src[i])
    if nvarsencoded != ctx.minfo.nvars:
        # some kind of internal error caused the number of vars in the format to differ from the number of vars in the ctx
        raise_(SIGSEGV)
    mpoly_set_monomial_ui(exp, user_exps, bits, ctx.minfo)


# Encode/decode to/from a memory buffer

ctypedef struct Buffer_structure:
    encoding_format * format
    ulong * buffer
    ulong buffer_size
    ulong count

cdef class Buffer:
    cdef MPolynomialRing_flint R
    cdef Py_ssize_t shape[2]
    cdef Py_ssize_t strides[2]
    cdef Buffer_structure buffer

    def __init__(self, R):
        self.R = R
        self.buffer.format = & self.R._encoding_format
        self.buffer.buffer = NULL
        self.buffer.buffer_size = 0
        self.buffer.count = 0

    def __dealloc__(self):
        if (self.buffer.buffer != NULL):
            free(self.buffer.buffer)

    def __getbuffer__(self, Py_buffer *buffer, int flags):

        self.shape[0] = self.buffer.count
        self.shape[1] = self.buffer.format.words
        self.strides[1] = sizeof(ulong)
        self.strides[0] = self.buffer.format.words * sizeof(ulong)

        buffer.buf = <char *> self.buffer.buffer
        buffer.format = 'q'                     # native aligned signed 64-bit integer
        buffer.internal = NULL
        buffer.itemsize = sizeof(ulong)
        buffer.len = self.buffer.count * self.buffer.format.words * sizeof(ulong)
        buffer.ndim = 2
        buffer.obj = self
        buffer.readonly = 1
        buffer.shape = self.shape
        buffer.strides = self.strides
        buffer.suboffsets = NULL

    def __releasebuffer__(self, Py_buffer *buffer):
        pass

cdef void encode_to_buffer(void * ptr, slong index, flint_bitcnt_t bits, ulong * exp, fmpz_t coeff, const fmpz_mpoly_ctx_t ctx) nogil:
    cdef Buffer_structure * buffer = <Buffer_structure *> ptr

    if index == 0:
        if buffer.buffer != NULL:
            free(buffer.buffer)
        buffer.buffer_size = 1024
        buffer.buffer = <ulong *>malloc(buffer.format.words * buffer.buffer_size * sizeof(ulong))
        buffer.count = 0
    if index == -1:
        return
    if <ulong>index >= buffer.buffer_size:
        buffer.buffer_size += 1024
        buffer.buffer = <ulong *>realloc(buffer.buffer, buffer.format.words * buffer.buffer_size * sizeof(ulong))

    encode_to_mem(buffer.format, buffer.buffer + buffer.format.words*buffer.count, bits, exp, coeff, ctx)

    buffer.count += 1

cdef void decode_from_buffer(void * ptr, ulong index, flint_bitcnt_t bits, ulong * exp, fmpz_t coeff, const fmpz_mpoly_ctx_t ctx) nogil:
    cdef Buffer_structure * buffer = <Buffer_structure *> ptr

    if index >= buffer.count:
        fmpz_set_ui(coeff, 0)
    else:
        decode_from_mem(buffer.format, buffer.buffer + buffer.format.words*index, bits, exp, coeff, ctx)

def copy_to_buffer(p):
    """
    TESTS::
        sage: from sage.rings.polynomial.multi_polynomial_flint import copy_to_buffer, copy_from_buffer
        sage: R.<x,y,z> = PolynomialRing(ZZ, implementation="FLINT", encoding="deglex64(3),sint64")
        sage: p = x^2 + y^2 + z^2
        sage: copy_to_buffer(p)
        Traceback (most recent call last):
        ...
        NotImplementedError: copy_to_buffer currently only works with 8 bit exponents
    """
    cdef MPolynomial_flint np = p
    cdef MPolynomialRing_flint parent = p.parent()
    cdef void * fptr[1]
    fptr[0] = <void *>np._poly
    buffer = Buffer(parent)
    cdef Buffer cbuffer = buffer
    with nogil:
        fmpz_mpoly_abstract_add(<void *> &cbuffer.buffer, fptr, 1, np._poly.bits, parent._ctx, NULL, encode_to_buffer)
    return buffer

def copy_from_buffer(buffer):
    """
    TESTS::
        sage: from sage.rings.polynomial.multi_polynomial_flint import copy_to_buffer, copy_from_buffer
        sage: R.<a,b,c,d,e,x,y,z> = PolynomialRing(ZZ, implementation="FLINT", order="lex", encoding="deglex64(8),sint64")
        sage: p = x^2 + y^2 + z^2
        sage: buffer = copy_to_buffer(p)
        sage: copy_from_buffer(buffer) == p
        True
        sage: import numpy
        sage: numpy.asarray(buffer)
        array([[14,  1],
               [11,  1],
               [ 9,  1]], dtype=int64)

        sage: R.<a,b,c,d,e,x,y,z> = PolynomialRing(ZZ, implementation="FLINT", order="lex", \
        ....:                                      encoding="deglex64(5),deglex64(3),sint64")
        sage: p = x^2 + y^2 + z^2
        sage: buffer = copy_to_buffer(p)
        sage: copy_from_buffer(buffer) == p
        True
        sage: numpy.asarray(buffer)
        array([[0,  9,  1],
               [0,  6,  1],
               [0,  4,  1]], dtype=int64)
        sage: numpy.asarray(copy_to_buffer(x))
        array([[0, 3, 1]], dtype=int64)

    Check that this is degree lexicographic encoding::

        sage: numpy.asarray(copy_to_buffer(x^2))
        array([[0, 9, 1]], dtype=int64)
        sage: numpy.asarray(copy_to_buffer(x*y))
        array([[0, 8, 1]], dtype=int64)
        sage: numpy.asarray(copy_to_buffer(x*z))
        array([[0, 7, 1]], dtype=int64)
        sage: numpy.asarray(copy_to_buffer(y^2))
        array([[0, 6, 1]], dtype=int64)
        sage: numpy.asarray(copy_to_buffer(y*z))
        array([[0, 5, 1]], dtype=int64)
        sage: numpy.asarray(copy_to_buffer(z^2))
        array([[0, 4, 1]], dtype=int64)
    """
    cdef Buffer cbuffer = buffer
    cdef MPolynomialRing_flint parent = cbuffer.R
    cdef MPolynomial_flint np = MPolynomial_flint.__new__(MPolynomial_flint)
    np._parent = parent
    cdef void * fptr[1]
    fptr[0] = <void *> &cbuffer.buffer
    with nogil:
        fmpz_mpoly_abstract_add(np._poly, fptr, 1, np._poly.bits, parent._ctx, decode_from_buffer, NULL)
    return np

# Encode/decode to/from a file
#
# Use stdio functions because they're faster, and because they don't require us to hold
# the GIL, so they can be used multi-threaded.
#
# Global variables, so we can only write to one file at a time :-(

ctypedef struct encode_to_file_struct:
    encoding_format * format
    FILE * FILE
    int fd
    ulong * buffer
    ulong buffer_size
    ulong count
    ulong total

cdef void encode_to_file(void * ptr, slong index, flint_bitcnt_t bits, ulong * exp, fmpz_t coeff, const fmpz_mpoly_ctx_t ctx) nogil:
    cdef encode_to_file_struct * state = <encode_to_file_struct *> ptr
    if index == -1:
        if state.count != 0:
            write(state.fd, state.buffer, state.format.words * state.count * sizeof(ulong))
        return
    elif index == 0:
        state.total = 0
    elif index % state.buffer_size == 0:
        write(state.fd, state.buffer, state.format.words * state.buffer_size * sizeof(ulong))
        state.count = 0
    encode_to_mem(state.format, state.buffer + state.format.words*state.count, bits, exp, coeff, ctx)
    state.count += 1
    state.total += 1

cdef open_file_for_encoding(encode_to_file_struct *state, filename):
    if filename.endswith('.gz'):
        command = "gzip > {}".format(filename)
        command_type = "w"
        state.FILE = popen(command.encode(), command_type.encode())
        if state.FILE == NULL:
            raise Exception("popen() failed on gzip > " + filename)
        # I'd prefer to dup the file descriptor and close the buffered popen
        # to avoid any kind of conflict between the buffered file I/O and
        # the raw I/O that we use, but pclose'ing (or fclose'ing) a popen'ed FILE
        # will wait for the process to terminate.
        state.fd = fileno(state.FILE)
        if state.fd == -1:
            raise Exception("fileno() failed on " + filename)
    else:
        state.fd = creat(filename.encode(), S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH)
        if state.fd == -1:
            raise Exception("creat() failed")
        state.FILE = NULL
    state.buffer_size = 1024
    state.count = 0
    state.buffer = <ulong *>malloc(state.format.words * state.buffer_size * sizeof(ulong))

cdef close_file_for_encoding(encode_to_file_struct *state):
    free(state.buffer)
    if state.FILE != NULL:
        retval = pclose(state.FILE)
        if retval == -1:
            raise Exception("pclose() failed")
        elif not WIFEXITED(retval):
            raise Exception("pclose() indicated abnormal exit of gzip")
        elif WEXITSTATUS(retval) != 0:
            raise Exception("gzip exit status " + str(WEXITSTATUS(retval)))
    else:
        close(state.fd)

# Decoding is actually a bit different from encoding, because we can decode in parallel.
#
# Encoding in parallel is not as useful, both because encoding is a bit faster (table
# lookup, but no table search), but mostly because encoding comes at the end of our
# processing pipeline and other things (like decoding 40 files being summed
# into only one, or multiplying polynomials) tend to dominate our runtime.

ctypedef struct decode_from_file_struct:
    encoding_format * format
    FILE * FILE
    int fd
    fmpz_mpoly_ctx_struct * ctx

    flint_bitcnt_t bits

    # `buffer` contains the bytes as they are read in from disk.
    #
    # 'count' polynomial terms (of format.words bytes each) have been read in so far total.
    # 'start' is the index of the first term in 'buffer'
    # There are (count-start) mod buffer_size terms in the buffer,
    #     and (count-start) mod buffer_size * format.words + trailing_bytes bytes in the buffer,
    #     which can never be more than buffer_size.
    # `trailing_bytes`, which can never be more than format.words, is the number of bytes
    #     in the partial term at the end of the buffer (if any).

    ulong * buffer
    ulong buffer_size
    ulong count
    ulong start
    ulong trailing_bytes

    # `exps` and `coeffs` are the decoded terms, ready to be fed to FLINT.
    #
    # They are divided into equal sized segments.  We allow for multiple threads to decode multiple
    # segments simultaneously
    #
    # We use two semaphores for each segment: one to indicate that it is free and we can begin writing
    # decoded terms into it, and one to indicate that it is complete and we can begin sending it to FLINT.
    #
    #                                Free  Disk  Ready
    #
    # Empty state (initial state)     +     -     -
    # Reading from disk               -     -     -
    # Ready to decode                 -     +     -
    # Decoding buffer                 -     -     -
    # Ready for FLINT                 -     -     +
    # Reading terms out               -     -     -
    # Finished reading out            +     -     - (back to empty state)
    #
    # We initialize everything as Free +, Disk -, Ready -
    #
    # READ LOOP:
    #
    #     Lock the mutex.
    #     If there is no disk reading thread running, no eof condition, and we're at the start of
    #         a segment, which means that (count-start) mod segment_size == 0, then
    #         sem_trywait that segment's Free and if we get it, set disk_reading_thread_active,
    #         release the mutex and try to read buffer_size/num_segments bytes from disk.
    #     If we read the entire segment, lock the mutex, post the segment's Disk,
    #         signal one thread (if any) waiting on the condition variable,
    #         unlock the mutex and move on to reading the next segment from disk.
    #     If we hit eof, write zeros into the next word, lock the mutex, set eof, post Disk,
    #         signal all threads waiting on the condition variable, and unlock the mutex.
    #     If we didn't read the entire buffer, keep reading until we finish it or get eof.
    #          (reading across buffers is more complicated, simplest is to do short reads)
    #     Goto START of READ LOOP
    #
    # DECODE LOOP:
    #
    #     Otherwise, run through the segments, ideally starting at the one after (count-start) mod buffer_size,
    #     sem_trywait'ing on Disk.  Keep going until we get one.  If we get none and eof is set,
    #     then unlock the mutex return.
    #     Otherwise, (the mutex is already locked) wait on the condition variable.  Afterwards, unlock the mutex.
    #     (or just go to the beginning of the decode loop with the mutex already locked)
    #
    #     Once we obtain a Disk segment and sem_trywait it, unlock the mutex and start decoding it.
    #
    #     Once a segment has been completed decoded, we post Ready.
    #     Goto START of DECODE LOOP
    #
    # The FLINT input_function knows where it is in the buffer because of its index variable being passed in.
    # It sem_wait's on a Ready every time it crosses a segment boundary.
    #
    # Since the FLINT function is blocking on Ready, we need at least one DECODE LOOP thread running above
    # and beyond the FLINT thread.
    #
    # The FLINT input_function could sem_trywait on a segment boundary and run decode loop if it couldn't
    # get Ready, but that would stall the entire FLINT input_function while a buffer was read or decoded.
    #
    # Best to have at least two threads running.
    #
    # If num_segments == 1, we'd be better off with the single-threaded code that decodes `buffer` right
    # into the locations provided to the FLINT input_function, and we simply don't use the next block
    # of variables, except for num_segments being one.

    int eof

    ulong num_segments
    ulong segment_size
    pthread_mutex_t mutex
    pthread_cond_t condvar
    sem_t * segments_free
    sem_t * segments_disk
    sem_t * segments_ready
    ulong * exps
    fmpz * coeffs

    ulong status_array[1024]
    ulong status_index


def decode_from_file_disk_reading_thread(ulong arg):
    cdef decode_from_file_struct * state = <decode_from_file_struct *> arg
    cdef ulong offset_in_segment
    cdef ulong segment_number
    cdef int at_end_of_segment
    cdef int retval

    with nogil:
        while not state.eof:
            offset_in_segment = state.count % state.segment_size
            segment_number = (state.count % state.buffer_size) / state.segment_size
            if offset_in_segment == 0:
                sem_wait(& state.segments_free[segment_number])
            if state.count == 0:
                state.trailing_bytes = 0
            while not state.eof:
                offset_in_segment = state.count % state.segment_size
                start_of_read = (<char *>state.buffer) + state.format.words * sizeof(ulong) * (state.count % state.buffer_size) + state.trailing_bytes
                length_of_read = state.format.words * sizeof(ulong) * (state.segment_size - offset_in_segment) - state.trailing_bytes
                retval = read(state.fd, start_of_read, length_of_read)
                if retval == -1:
                    raise_(SIGSEGV)
                if retval == 0:
                    state.eof = 1
                    if state.trailing_bytes != 0:
                        raise_(SIGSEGV)
                else:
                    retval += state.trailing_bytes
                    state.trailing_bytes = retval % (state.format.words * sizeof(ulong))
                    state.count += retval / (state.format.words * sizeof(ulong))
                    at_end_of_segment = ((state.count % state.segment_size) == 0)
                    if at_end_of_segment:
                        break
            pthread_mutex_lock(& state.mutex)
            # state.status_array[state.status_index] = segment_number + 1000*state.count
            # state.status_index += 1
            # with gil:
            #     print("posting segments_disk", segment_number, "count", state.count, "eof", state.eof)
            sem_post(& state.segments_disk[segment_number])
            if state.eof:
                pthread_cond_broadcast(& state.condvar)
            else:
                # Wake a thread to decode the segment we just read
                pthread_cond_signal(& state.condvar)
            pthread_mutex_unlock(& state.mutex)

cpdef decode_from_file_decoding_thread(ulong arg):
    cdef decode_from_file_struct * state = <decode_from_file_struct *> arg
    cdef ulong offset
    cdef ulong offset_in_segment
    cdef ulong segment
    cdef ulong segment_number
    cdef ulong decode_count
    cdef ulong i
    cdef slong N = mpoly_words_per_exp(state.bits, state.ctx.minfo)

    with nogil:
        pthread_mutex_lock(& state.mutex)

        while True:
            offset_in_segment = state.count % state.segment_size
            segment_number = (state.count % state.buffer_size) / state.segment_size
            for offset in range(state.num_segments):
                # pick the first segment ready to decode after state.count (location of disk read)
                # would be better to pick the first segment ready after last FLINT read
                segment = (segment_number + 1 + offset) % state.num_segments
                if sem_trywait(& state.segments_disk[segment]) == 0:
                    if not state.eof or segment_number != segment:
                        decode_count = state.segment_size
                    else:
                        decode_count = offset_in_segment
                    # state.status_array[state.status_index] = segment + 100
                    # state.status_index += 1
                    pthread_mutex_unlock(& state.mutex)
                    # decode the segment
                    # with gil:
                    #     print("decoding", segment, decode_count)
                    for i in range(decode_count):
                        decode_from_mem(state.format, state.buffer + state.format.words * (i + segment * state.segment_size),
                                        state.bits, state.exps + N * (i + segment * state.segment_size),
                                        state.coeffs + (i + segment * state.segment_size), state.ctx)

                    sem_post(& state.segments_ready[segment])
                    pthread_mutex_lock(& state.mutex)
                    break
                elif errno != EAGAIN:
                    raise_(SIGSEGV)
            else:
                if state.eof:
                    break
                pthread_cond_wait(& state.condvar, & state.mutex)

        pthread_mutex_unlock(& state.mutex)

cdef void decode_from_file(void * ptr, ulong index, flint_bitcnt_t bits, ulong * exp, fmpz_t coeff, const fmpz_mpoly_ctx_t ctx) nogil:
    cdef decode_from_file_struct * state = <decode_from_file_struct *> ptr
    cdef slong N = mpoly_words_per_exp(bits, ctx.minfo)
    cdef unsigned char * exps = <unsigned char *> exp
    cdef int retval

    if state.num_segments > 1:
        if (index % state.segment_size == 0):
            sem_wait(& state.segments_ready[(index % state.buffer_size) / state.segment_size])
        if index == state.count:
            if not state.eof:
                raise_(SIGSEGV)
            fmpz_set_ui(coeff, 0)
        else:
            fmpz_set(coeff, state.coeffs + (index % state.buffer_size))
            mpoly_monomial_set(exp, state.exps + N*(index % state.buffer_size), N)
        if ((index + 1) % state.segment_size == 0):
            # with gil:
            #     print("posting segments_free", (index % state.buffer_size) / state.segment_size)
            pthread_mutex_lock(& state.mutex)
            # state.status_array[state.status_index] = ((index % state.buffer_size) / state.segment_size) + 200
            # state.status_index += 1
            sem_post(& state.segments_free[(index % state.buffer_size) / state.segment_size])
            pthread_mutex_unlock(& state.mutex)
        return

    while index == state.count:
        if index == 0:
            state.trailing_bytes = 0
        if state.trailing_bytes > 0:
            bcopy(state.buffer + state.format.words * (state.count - state.start), state.buffer, state.trailing_bytes)
        state.start = state.count
        retval = read(state.fd, (<char *>state.buffer) + state.trailing_bytes,
                      state.format.words * state.buffer_size * sizeof(ulong) - state.trailing_bytes)
        if retval == -1:
            raise_(SIGSEGV)
        if retval == 0:
            fmpz_set_si(coeff, 0)
            return
        retval += state.trailing_bytes
        state.trailing_bytes = retval % (state.format.words * sizeof(ulong))
        state.count += retval / (state.format.words * sizeof(ulong))

    decode_from_mem(state.format, state.buffer + state.format.words*(index-state.start), bits, exp, coeff, ctx)

cdef open_file_for_decoding(decode_from_file_struct *state, filename, num_decoding_threads=1):
    """
    `state` needs to have its `format`, `bits`, and `ctx` fields set correctly.  Everything else
    gets filled in here, including malloc'ing a buffer.
    """
    cdef ulong i

    if filename.endswith('.gz'):
        command = "zcat {}".format(filename)
        command_type = "r"
        state.FILE = popen(command.encode(), command_type.encode())
        if state.FILE == NULL:
            raise Exception("popen() failed on zcat " + filename)
        # I'd prefer to dup the file descriptor and close the buffered popen
        # to avoid any kind of conflict between the buffered file I/O and
        # the raw I/O that we use, but pclose'ing (or fclose'ing) a popen'ed FILE
        # will wait for the process to terminate.
        state.fd = fileno(state.FILE)
        if state.fd == -1:
            raise Exception("fileno() failed on " + filename)
    else:
        state.fd = open(filename.encode(), O_RDONLY)
        if state.fd == -1:
            raise Exception("open() failed on " + filename)
        state.FILE = NULL
    state.num_segments = num_decoding_threads
    state.segment_size = 1024
    state.buffer_size = state.num_segments * state.segment_size
    state.start = 0
    state.count = 0
    state.buffer = <ulong *>malloc(state.format.words * state.buffer_size * sizeof(ulong))

    cdef slong N = mpoly_words_per_exp(state.bits, state.ctx.minfo)

    state.status_index = 0

    state.eof = 0
    state.segments_free = <sem_t *> malloc(state.num_segments * sizeof(sem_t))
    state.segments_disk = <sem_t *> malloc(state.num_segments * sizeof(sem_t))
    state.segments_ready = <sem_t *> malloc(state.num_segments * sizeof(sem_t))
    state.exps = <ulong *>malloc(state.buffer_size * N * sizeof(ulong))
    state.coeffs = <fmpz *>malloc(state.buffer_size * sizeof(fmpz))
    for i in range(state.buffer_size):
        fmpz_init(state.coeffs + i)
    pthread_mutex_init(& state.mutex, NULL)
    pthread_cond_init(& state.condvar, NULL)
    for i in range(state.num_segments):
        sem_init(& state.segments_free[i], 0, 1)
        sem_init(& state.segments_disk[i], 0, 0)
        sem_init(& state.segments_ready[i], 0, 0)
    if state.num_segments > 1:
        for i in range(state.num_segments):
            # XXX need to set state.bits and state.ctx before these threads start
            th = threading.Thread(target = decode_from_file_decoding_thread, args = (<ulong> state,))
            th.start()
        # XXX Ideally, I'd now like to block until all of the decoding threads are waiting on state.condvar
        th = threading.Thread(target = decode_from_file_disk_reading_thread, args = (<ulong> state,))
        th.start()

cdef close_file_for_decoding(decode_from_file_struct *state):
    """
    `state`'s buffer will also get free'd by this function.
    """
    free(state.buffer)
    state.buffer = NULL
    for i in range(state.num_segments):
        sem_destroy(& state.segments_free[i])
        sem_destroy(& state.segments_disk[i])
        sem_destroy(& state.segments_ready[i])
    pthread_cond_destroy(& state.condvar)
    pthread_mutex_destroy(& state.mutex)
    free(state.segments_free)
    free(state.segments_disk)
    free(state.segments_ready)
    if state.exps != NULL:
        free(state.exps)
    free(state.coeffs)
    state.segments_free = NULL
    state.segments_disk = NULL
    state.segments_ready = NULL
    state.exps = NULL
    state.coeffs = NULL

    if state.FILE != NULL:
        retval = pclose(state.FILE)
        if retval == -1:
            raise Exception("pclose() failed")
        elif not WIFEXITED(retval):
            raise Exception("pclose() indicated abnormal exit of gzip")
        elif WEXITSTATUS(retval) != 0:
            raise Exception("gzip exit status " + str(WEXITSTATUS(retval)))
    else:
        close(state.fd)

def copy_to_file(p, filename="bigflint.out"):
    cdef MPolynomial_flint np = p
    cdef MPolynomialRing_flint parent = p.parent()
    cdef const void * fptr[1]
    cdef encode_to_file_struct state
    fptr[0] = <void *>np._poly
    state.format = & parent._encoding_format
    open_file_for_encoding(& state, filename)
    with nogil:
        fmpz_mpoly_abstract_add(& state, fptr, 1, 8, parent._ctx, NULL, encode_to_file)
    close_file_for_encoding(& state)

def copy_from_file(R, filename="bigflint.out"):
    """
    TESTS::
        sage: from sage.rings.polynomial.multi_polynomial_flint import copy_to_file, copy_from_file
        sage: R.<x,y,z> = PolynomialRing(ZZ, implementation="FLINT")
        sage: copy_from_file(R, filename=''.join(chr(randrange(65,91)) for _ in range(250)))
        Traceback (most recent call last):
        ...
        Exception: open() failed on ...
        sage: copy_from_file(R, filename=''.join(chr(randrange(65,91)) for _ in range(250)) + '.gz')
        Traceback (most recent call last):
        ...
        Exception: gzip exit status 1

        sage: R.<a,b,c,d,e,x,y,z> = PolynomialRing(ZZ, implementation="FLINT", order="lex", encoding="deglex64(8),sint64")
        sage: p = x^2 + y^2 + z^2
        sage: filename='/tmp/' + ''.join(chr(randrange(65,91)) for _ in range(250))
        sage: copy_to_file(p, filename=filename)
        sage: p2 = copy_from_file(R, filename=filename)
        sage: p == p2
        True
        sage: os.unlink(filename)

        sage: filename='/tmp/' + ''.join(chr(randrange(65,91)) for _ in range(250)) + '.gz'
        sage: copy_to_file(p, filename=filename)
        sage: p2 = copy_from_file(R, filename=filename)
        sage: p == p2
        True
        sage: os.unlink(filename)
    """
    cdef MPolynomialRing_flint parent = R
    cdef MPolynomial_flint np = MPolynomial_flint.__new__(MPolynomial_flint)
    np._parent = R
    cdef const void * fptr[1]
    cdef decode_from_file_struct state

    state.format = & parent._encoding_format
    state.bits = np._poly.bits
    state.ctx = & parent._ctx[0]
    open_file_for_decoding(& state, filename)
    fptr[0] = <void *> &state
    with nogil:
        fmpz_mpoly_abstract_add(np._poly, fptr, 1, np._poly.bits, parent._ctx, decode_from_file, NULL)
    close_file_for_decoding(& state)
    return np

# Substitute variables by reading a polynomial in an input file in
# blocks with each block having identical substituted variables.

cdef ulong substitute_last_radii = UINT64_MAX
cdef ulong radii_block_size = 0
cdef ulong radii_count = 0
cdef ulong * radii_exp_block = NULL
cdef ulong * radii_coeff_block = NULL

cdef const char * encode_to_file_returning_status(void * ptr, slong index, flint_bitcnt_t bits, ulong * exp, fmpz_t coeff, const fmpz_mpoly_ctx_t ctx) nogil:
    cdef encode_to_file_struct * state = <encode_to_file_struct *> ptr
    encode_to_file(ptr, index, bits, exp, coeff, ctx)
    return status_string_ptr

# Once these blocks have been passed to this function, they get processed and then freed.

ctypedef struct output_block_data:
    ulong radii
    ulong * exp
    fmpz * coeffs
    ulong count
    flint_bitcnt_t bits
    fmpz_mpoly_ctx_t ctx
    encoding_format * format

radii_info = {}

substitute_threads = []
substitute_thread_limit = 12
substitute_thread_semaphore = threading.Semaphore(substitute_thread_limit)

# We can't pass C pointers to Python functions, but we need a Python
# function to pass to threading.Thread, so we convert the pointer to a
# ulong, pass it to this function, and then cast it back to a pointer.

cpdef substitute_output_block(ulong arg1):
    cdef output_block_data * data = <output_block_data *> arg1;
    global r1poly, r2poly, r12poly
    cdef fmpz_mpoly_struct poly1
    cdef MPolynomial_flint poly2
    cdef void * fptr[2]
    cdef slong iptr[1]
    cdef encode_to_file_struct state
    cdef int r1_power, r2_power, r12_power

    filename = "radii-{}.out.gz".format(data.radii)

    if not os.path.isfile(filename):

        poly1.coeffs = data.coeffs
        poly1.exps = data.exp
        poly1.length = data.count
        poly1.bits = data.bits

        # discard LSBs; they were presevered below when the exponents were stored
        r1_power = ((data.radii >> 16) & 254) / 2
        r2_power = ((data.radii >> 8) & 254) / 2
        r12_power = (data.radii & 254) / 2
        poly2 = r1poly**r1_power * r2poly**r2_power * r12poly**r12_power

        fptr[0] = <void *> & poly1
        fptr[1] = <void *> poly2._poly
        iptr[0] = 2

        state.format = data.format
        open_file_for_encoding(& state, filename)

        radii_info[data.radii] = (data.count, poly2._poly.length)

        with nogil:
            fmpz_mpoly_addmul_multi_threaded_abstract(<void *> &state, <const fmpz_mpoly_struct **> fptr, iptr, 1, data.ctx, encode_to_file_returning_status)

        radii_info[data.radii] += (state.total,)

        close_file_for_encoding(& state)

    else:

        print("Skipping radii", data.radii, "(file exists)", file=sys.stderr)

    free(data.coeffs)
    free(data.exp)
    free(data)
    substitute_thread_semaphore.release()

# Output function for the substitution routine

cdef void substitute_output_function(void * format, slong index, flint_bitcnt_t bits, ulong * exp, fmpz_t coeff, const fmpz_mpoly_ctx_t ctx) nogil:
    global substitute_last_radii
    global substitute_threads, substitute_thread_limit, substitute_thread_semaphore
    global radii_block_size, radii_count, radii_exp_block, radii_coeff_block
    cdef unsigned char * exps
    cdef output_block_data * data

    cdef slong N = mpoly_words_per_exp(bits, ctx.minfo)
    cdef ulong current_radii

    if index == -1:
        current_radii = UINT64_MAX
    else:
        current_radii = (exp[15] >> 56) | (exp[16] << 8)

    if (current_radii != substitute_last_radii):
        if substitute_last_radii != UINT64_MAX:
            # start a new thread to multiply and output polynomial
            with gil:
                data = <output_block_data *> malloc(sizeof(output_block_data))
                data.radii = substitute_last_radii
                data.exp = radii_exp_block
                data.coeffs = <fmpz *> radii_coeff_block
                data.count = radii_count
                data.bits = bits
                data.ctx = ctx
                data.format = <encoding_format *> format
                substitute_thread_semaphore.acquire()
                th = threading.Thread(target = substitute_output_block, args = (<ulong> data,))
                th.start()
                substitute_threads.append(th)
            radii_exp_block = NULL
            radii_coeff_block = NULL
            radii_block_size = 0
        substitute_last_radii = current_radii
        radii_count = 0

    if radii_count >= radii_block_size:
        if radii_block_size == 0:
            radii_block_size = 1024*1024
            radii_exp_block = <ulong *>malloc(radii_block_size * N * sizeof(ulong))
            radii_coeff_block = <ulong *>malloc(radii_block_size * sizeof(ulong))
        else:
            radii_block_size += 1024*1024
            radii_exp_block = <ulong *>realloc(radii_exp_block, radii_block_size * N * sizeof(ulong))
            radii_coeff_block = <ulong *>realloc(radii_coeff_block, radii_block_size * sizeof(ulong))

    if index == -1:
        if radii_exp_block != NULL: free(radii_exp_block)
        if radii_coeff_block != NULL: free(radii_coeff_block)
        radii_exp_block = NULL
        radii_coeff_block = NULL
        radii_block_size = 0
    else:
        mpoly_monomial_set(radii_exp_block + radii_count*N, exp, N)
        radii_coeff_block[radii_count] = (<ulong *>coeff)[0]

        # mask off all but LSBs of r1, r2, and r12
        exps = <unsigned char *> (radii_exp_block + radii_count*N)
        exps[127] &= 1
        exps[128] &= 1
        exps[129] &= 1

        radii_count += 1

def substitute_file(R, filename="bigflint.out"):
    global r1poly, r2poly, r12poly
    global substitute_threads
    x1 = R.gens_dict()['x1']
    y1 = R.gens_dict()['y1']
    z1 = R.gens_dict()['z1']
    x2 = R.gens_dict()['x2']
    y2 = R.gens_dict()['y2']
    z2 = R.gens_dict()['z2']
    r1poly = (x1**2 + y1**2 + z1**2)
    r2poly = (x2**2 + y2**2 + z2**2)
    r12poly = ((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
    cdef MPolynomialRing_flint parent = R

    cdef const void * fptr[1]
    cdef decode_from_file_struct state
    state.format = & parent._encoding_format
    state.bits = 8
    state.ctx = & parent._ctx[0]
    open_file_for_decoding(& state, filename)
    fptr[0] = <void *> &state

    with nogil:
        fmpz_mpoly_abstract_add(<void *> & parent._encoding_format, fptr, 1, 8, parent._ctx, decode_from_file, substitute_output_function)

    close_file_for_decoding(& state)

    for th in substitute_threads: th.join()
    substitute_threads = []

cdef ulong of2_count = 0

cdef ulong last_radii = 0
cdef ulong last_radii_count = 0
cdef ulong max_radii_count = 0
cdef ulong radii_blocks = 0

cdef ulong * last_exp = NULL
cdef ulong * max_exp = NULL

cdef int max_vdeg = 0
cdef int max_cdeg = 0

cdef const char * verify_fcn_returning_status(void * ptr, slong index, flint_bitcnt_t bits, ulong * exp, fmpz_t coeff, const fmpz_mpoly_ctx_t ctx) nogil:
    global status_string, status_string_encode, status_string_ptr
    global last_radii, last_radii_count, radii_blocks, max_radii_count
    global last_exp, max_exp
    global max_cdeg, max_vdeg
    global of2_count

    if index == -1:
        if last_radii_count > max_radii_count:
            max_radii_count = last_radii_count
        with gil:
            status_string = "{}/{} vdeg={} cdeg={}".format(radii_blocks, max_radii_count, max_vdeg, max_cdeg)
            status_string_encode = status_string.encode()
            status_string_ptr = status_string_encode
        return status_string_ptr

    if of2_count != <ulong>index:
        raise_(SIGSEGV)
    of2_count += 1

    cdef slong N = mpoly_words_per_exp(bits, ctx.minfo)

    # Check to see if monomial exponents are ordered correctly
    # XXX doesn't resize properly if verify_file is called twice
    if last_exp == NULL:
        last_exp = <ulong *>malloc(N * sizeof(ulong))
        mpoly_monomial_set(last_exp, exp, N)
    else:
        if mpoly_monomial_lt_nomask(last_exp, exp, N):
            with gil:
                print("mpoly_monomial_lt_nomask failed at index", index)
            raise_(SIGSEGV)
        else:
            mpoly_monomial_set(last_exp, exp, N)

    if max_exp == NULL:
        max_exp = <ulong *>malloc(N * sizeof(ulong))
        mpoly_monomial_set(max_exp, exp, N)
    else:
        mpoly_monomial_max(max_exp, max_exp, exp, bits, N, mpoly_overflow_mask_sp(bits))

    cdef unsigned char * exps = <unsigned char *> exp
    cdef int vdeg = 0
    cdef int cdeg = 0
    cdef int i
    for i in range(11):
        vdeg += exps[129-i]
    for i in range(11,130):
        cdeg += exps[129-i]

    cdef ulong current_radii = (exp[15] >> 56) | (exp[16] << 8)
    if (current_radii != last_radii):
        if (last_radii != 0) and (current_radii > last_radii):
            with gil:
                print("current_radii > last_radii failed at index", index)
            raise_(SIGSEGV)
        last_radii = current_radii
        radii_blocks += 1
        if last_radii_count > max_radii_count:
            max_radii_count = last_radii_count
        last_radii_count = 1
    else:
        last_radii_count += 1

    if (last_radii_count == 1) or (vdeg > max_vdeg) or (cdeg > max_cdeg):
        if vdeg > max_vdeg: max_vdeg = vdeg
        if cdeg > max_cdeg: max_cdeg = cdeg
        exp1 = exps[129]
        exp2 = exps[128]
        exp3 = exps[127]
        with gil:
            exps = <unsigned char *> max_exp
            max_exps = ''.join(map(str, exps[0:119]))
            status_string = "{},{},{} {}/{} vdeg={} cdeg={} {}".format(exp1, exp2, exp3, radii_blocks, max_radii_count, max_vdeg, max_cdeg, max_exps)
            status_string_encode = status_string.encode()
            status_string_ptr = status_string_encode

    return status_string_ptr

cdef void verify_fcn(void * ptr, slong index, flint_bitcnt_t bits, ulong * exp, fmpz_t coeff, const fmpz_mpoly_ctx_t ctx) nogil:
    global of2_count
    verify_fcn_returning_status(NULL, index, bits, exp, coeff, ctx)
    if ((index > 0) and (index % 1000000 == 0)) or (index == -1):
        with gil:
            if index == -1:
                print("Output length", of2_count, status_string)
            else:
                print("Output length", index, status_string, end='\r')

def verify_file(R, filename="bigflint.out"):
    cdef MPolynomialRing_flint parent = R

    cdef void * fptr[1]
    cdef decode_from_file_struct state
    state.format = & parent._encoding_format
    # XXX bits is hardwired here
    state.bits = 8
    state.ctx = & parent._ctx[0]
    open_file_for_decoding(& state, filename)
    fptr[0] = <void *> &state

    with nogil:
        fmpz_mpoly_abstract_add(NULL, fptr, 1, 8, parent._ctx, decode_from_file, verify_fcn)

    close_file_for_decoding(&state)

# The buffer is divided into equal sized segments
#
# We use two semaphores - one to indicate that the load function can advance to the next segment,
# and another to indicate that the load function is done with a segment and it can be freed.
#
# We initialize segments_free = num_segments and segments_ready_to_load = 0
#
# Once a segment has been completed decoded, we post segments_ready_to_load.  The load
# function waits on it and advances into the next segment once it gets posted.
#
# Once a segment has been completed loaded, we post segments_free.  The decode function
# waits on it and starts decoding the next segment once it gets posted.

ctypedef struct load_from_decoded_buffer_struct:
    sem_t segments_ready_to_load
    sem_t segments_free
    ulong * exps
    fmpz * coeffs
    ulong buffer_size
    ulong num_segments

cdef void load_from_decoded_buffer(void * poly, ulong index, flint_bitcnt_t bits, ulong * exp, fmpz_t coeff, const fmpz_mpoly_ctx_t ctx) nogil:
    cdef load_from_decoded_buffer_struct * state = <load_from_decoded_buffer_struct *> poly
    cdef slong N = mpoly_words_per_exp(bits, ctx.minfo)

    if index % (state.buffer_size / state.num_segments) == 0:
        sem_wait(& state.segments_ready_to_load)

    mpoly_monomial_set(exp, state.exps + N*(index % state.buffer_size), N)
    fmpz_set(coeff, state.coeffs + (index % state.buffer_size))

    if (index+1) % (state.buffer_size / state.num_segments) == 0:
        #sem_post(& state.segments_free[(index % state.buffer_size) % state.num_segments])
        sem_post(& state.segments_free)

ctypedef struct read_and_decode_file_data:
    load_from_decoded_buffer_struct * loader_state
    fmpz_mpoly_ctx_t ctx
    encoding_format * format
    int fd
    ulong buffer_size
    flint_bitcnt_t bits
    slong N
    int eof

    ulong * buffer
    ulong count
    ulong start
    ulong index
    ulong trailing_bytes

cdef void read_and_decode_one_buffer(read_and_decode_file_data * state) nogil:
    cdef int retval
    cdef ulong * exp

    while not state.eof and state.index == state.count:
        if state.trailing_bytes > 0:
            bcopy(state.buffer + state.format.words * (state.count - state.start), state.buffer, state.trailing_bytes)
        state.start = state.count
        retval = read(state.fd, (<char *>state.buffer) + state.trailing_bytes,
                      state.format.words * state.buffer_size * sizeof(ulong) - state.trailing_bytes)
        if retval == -1:
            raise_(SIGSEGV)
        if retval == 0:
            state.eof = 1
        else:
            retval += state.trailing_bytes
            state.trailing_bytes = retval % (state.format.words * sizeof(ulong))
            state.count += retval / (state.format.words * sizeof(ulong))

    while state.index < state.count:

        if state.index % (state.loader_state.buffer_size / state.loader_state.num_segments) == 0:
            sem_wait(& state.loader_state.segments_free)

        exp = state.loader_state.exps + state.N*(state.index % state.loader_state.buffer_size)

        decode_from_mem(state.format, state.buffer + state.format.words*(state.index-state.start), state.bits, exp,
                        state.loader_state.coeffs + (state.index % state.loader_state.buffer_size), state.ctx)

        if (state.index+1) % (state.loader_state.buffer_size / state.loader_state.num_segments) == 0:
            sem_post(& state.loader_state.segments_ready_to_load)

        state.index += 1

    if state.eof:

        if state.index % (state.loader_state.buffer_size / state.loader_state.num_segments) == 0:
            sem_wait(& state.loader_state.segments_free)

        fmpz_set_si(state.loader_state.coeffs + (state.index % state.loader_state.buffer_size), 0)

        sem_post(& state.loader_state.segments_ready_to_load)


cdef read_and_decode_file(read_and_decode_file_data * state):
    cdef int retval
    cdef ulong * exp
    cdef unsigned char * exps

    state.eof = 0
    state.trailing_bytes = 0
    state.start = 0
    state.count = 0
    state.index = 0
    state.buffer = <ulong *> malloc(3 * state.buffer_size * sizeof(ulong))

    with nogil:
        while not state.eof:
            read_and_decode_one_buffer(state)
        free(state.buffer)

cpdef read_and_decode_file_py(ulong arg1):
    cdef read_and_decode_file_data * data = <read_and_decode_file_data *> arg1;
    read_and_decode_file(data)

def sum_files(R, filename_list=[], filename=None):
    cdef MPolynomialRing_flint parent = R

    cdef int nfiles = len(filename_list)
    if nfiles == 0: return

    cdef const fmpz_mpoly_struct ** fptr = <const fmpz_mpoly_struct **>malloc(sizeof(fmpz_mpoly_struct *) * nfiles)
    cdef load_from_decoded_buffer_struct * state = <load_from_decoded_buffer_struct *> malloc(sizeof(load_from_decoded_buffer_struct) * nfiles)
    cdef read_and_decode_file_data * data = <read_and_decode_file_data *> malloc(sizeof(read_and_decode_file_data) * nfiles)
    cdef FILE ** popen_FILE = <FILE **> malloc(sizeof(FILE *) * nfiles)

    cdef encode_to_file_struct encoding_state
    cdef FILE * encoding_FILE = NULL
    cdef ulong j

    if filename != None:
        if filename.endswith('.gz'):
            command = "gzip > {}".format(filename)
            command_type = "w"
            encoding_FILE = popen(command.encode(), command_type.encode())
            if encoding_FILE == NULL:
                raise Exception("popen() failed on gzip > " + filename)
            encoding_state.fd = fileno(encoding_FILE)
            if encoding_state.fd == -1:
                raise Exception("fileno() failed on " + filename)
        else:
            encoding_state.fd = creat(filename.encode(), S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH)
            if encoding_state.fd == -1:
                raise Exception("creat() failed on " + filename)
        encoding_state.buffer = <ulong *>malloc(3 * 1024 * sizeof(ulong))
        encoding_state.buffer_size = 1024
        encoding_state.count = 0
    else:
        encoding_state.fd = -1
        encoding_state.buffer = NULL
        encoding_state.buffer_size = 0
        encoding_state.count = 0

    # bits is hardwired
    N = mpoly_words_per_exp(8, parent._ctx.minfo)

    threads = []

    for i, input_filename in enumerate(filename_list):
        if input_filename.endswith('.gz'):
            command = "zcat {}".format(input_filename)
            command_type = "r"
            popen_FILE[i] = popen(command.encode(), command_type.encode())
            if popen_FILE[i] == NULL:
                # XXX doesn't free linguring buffers
                raise Exception("popen() failed on zcat " + input_filename)
            data[i].fd = fileno(popen_FILE[i])
        else:
            data[i].fd = open(input_filename.encode(), O_RDONLY)
            if data[i].fd == -1:
                raise Exception("open() failed on " + input_filename)
            popen_FILE[i] = NULL

        state[i].buffer_size = 4 * 1024
        state[i].num_segments = 4
        state[i].exps = <ulong *>malloc(N * state[i].buffer_size * sizeof(ulong))
        state[i].coeffs = <fmpz *>malloc(state[i].buffer_size * sizeof(fmpz))
        sem_init(& state[i].segments_ready_to_load, 0, 0)
        sem_init(& state[i].segments_free, 0, 4)

        for j in range(state[i].buffer_size):
            fmpz_init(state[i].coeffs + j)

        data[i].loader_state = & state[i]
        data[i].buffer_size = 1024
        # XXX bits is hardwired
        data[i].bits = 8
        data[i].N = N

        fptr[i] = <fmpz_mpoly_struct *> & state[i]

        th = threading.Thread(target = read_and_decode_file_py, args = (<ulong> &(data[i]),))
        th.start()
        threads.append(th)

    with nogil:
        # XXX bits is hardwired
        if encoding_state.buffer == NULL:
            fmpz_mpoly_abstract_add(NULL, <void **> fptr, nfiles, 8, parent._ctx, load_from_decoded_buffer, verify_fcn)
        else:
            fmpz_mpoly_abstract_add(<void *> &encoding_state, <void **> fptr, nfiles, 8, parent._ctx, load_from_decoded_buffer, encode_to_file)

    for th in threads:
        th.join()

    if encoding_state.buffer != NULL:
        free(encoding_state.buffer)
        if encoding_FILE != NULL:
            pclose(encoding_FILE)
        else:
            close(encoding_state.fd)

    for i in range(nfiles):
        if popen_FILE[i] == NULL:
            close(data[i].fd)
        else:
            pclose(popen_FILE[i])
        free(state[i].exps)
        free(state[i].coeffs)

    free(popen_FILE)
    free(state)
    free(data)
    free(fptr)

cdef class MPolynomialRing_flint(MPolynomialRing_base):

    def __init__(self, base_ring, n, names, order='degrevlex', encoding=None):
        """
        Construct a multivariate polynomial ring subject to the
        following conditions:

        INPUT:

        - ``base_ring`` - base ring (must be ZZ)

        - ``n`` - number of variables (must be at least 1)

        - ``names`` - names of ring variables, may be string of list/tuple

        - ``order`` - term order (must be ``degrevlex``; FLINT also supports ``lex`` and ``deglex``)

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: P
            Multivariate Polynomial Ring in x, y, z over Rational Field

            sage: f = 27/113 * x^2 + y*z + 1/2; f
            27/113*x^2 + y*z + 1/2

            sage: P.term_order()
            Degree reverse lexicographic term order

            sage: P = PolynomialRing(GF(127),3,names='abc', order='lex')
            sage: P
            Multivariate Polynomial Ring in a, b, c over Finite Field of size 127

            sage: a,b,c = P.gens()
            sage: f = 57 * a^2*b + 43 * c + 1; f
            57*a^2*b + 43*c + 1

            sage: P.term_order()
            Lexicographic term order

            sage: z = QQ['z'].0
            sage: K.<s> = NumberField(z^2 - 2)
            sage: P.<x,y> = PolynomialRing(K, 2)
            sage: 1/2*s*x^2 + 3/4*s
            (1/2*s)*x^2 + (3/4*s)

            sage: P.<x,y,z> = ZZ[]; P
            Multivariate Polynomial Ring in x, y, z over Integer Ring

            sage: P.<x,y,z> = Zmod(2^10)[]; P
            Multivariate Polynomial Ring in x, y, z over Ring of integers modulo 1024

            sage: P.<x,y,z> = Zmod(3^10)[]; P
            Multivariate Polynomial Ring in x, y, z over Ring of integers modulo 59049

            sage: P.<x,y,z> = Zmod(2^100)[]; P
            Multivariate Polynomial Ring in x, y, z over Ring of integers modulo 1267650600228229401496703205376

            sage: P.<x,y,z> = Zmod(2521352)[]; P
            Multivariate Polynomial Ring in x, y, z over Ring of integers modulo 2521352
            sage: type(P)
            <type 'sage.rings.polynomial.multi_polynomial_libsingular.MPolynomialRing_libsingular'>

            sage: P.<x,y,z> = Zmod(25213521351515232)[]; P
            Multivariate Polynomial Ring in x, y, z over Ring of integers modulo 25213521351515232
            sage: type(P)
            <class 'sage.rings.polynomial.multi_polynomial_ring.MPolynomialRing_polydict_with_category'>

            sage: P.<x,y,z> = PolynomialRing(Integers(2^32),order='lex')
            sage: P(2^32-1)
            4294967295

        TESTS:

        Make sure that a faster coercion map from the base ring is used;
        see :trac:`9944`::

            sage: R.<x,y> = PolynomialRing(ZZ)
            sage: R.coerce_map_from(R.base_ring())
            Polynomial base injection morphism:
              From: Integer Ring
              To:   Multivariate Polynomial Ring in x, y over Integer Ring

        Check some invalid input::

            sage: from sage.rings.polynomial.multi_polynomial_libsingular import MPolynomialRing_libsingular
            sage: MPolynomialRing_libsingular(Zmod(1), 1, ["x"], "lex")
            Traceback (most recent call last):
            ...
            NotImplementedError: polynomials over Ring of integers modulo 1 are not supported in Singular
            sage: MPolynomialRing_libsingular(SR, 1, ["x"], "lex")
            Traceback (most recent call last):
            ...
            NotImplementedError: polynomials over Symbolic Ring are not supported in Singular
            sage: MPolynomialRing_libsingular(QQ, 0, [], "lex")
            Traceback (most recent call last):
            ...
            NotImplementedError: polynomials in 0 variables are not supported in Singular
            sage: MPolynomialRing_libsingular(QQ, -1, [], "lex")
            Traceback (most recent call last):
            ...
            NotImplementedError: polynomials in -1 variables are not supported in Singular
        """

        if base_ring != ZZ:
            raise NotImplementedError("FLINT base_ring must be ZZ")

        if encoding != None:
            encoded_vars = []
            encoded_coeff = False
            for encoding_item in encoding.split(','):
                match = re.fullmatch('deglex64\\(([0-9]+)\\)', encoding_item)
                if match:
                    if encoded_coeff:
                        raise NotImplementedError('FLINT coefficient encodings must currently be the last thing encoded')
                    if int(match[1]) == 0:
                        raise ValueError('FLINT variable encodings can not encode zero variables')
                    encoded_vars.append(int(match[1]))
                elif encoding_item == 'sint64':
                    if encoded_coeff:
                        raise ValueError('FLINT encodings must include exactly one coefficient encoding')
                    encoded_coeff = True
                else:
                    raise ValueError('FLINT encodings must be deglex64 or sint64')
            if not encoded_coeff:
                raise ValueError('FLINT encodings must include exactly one coefficient encoding')
            if sum(encoded_vars) != n:
                raise ValueError('FLINT encodings must encode exactly the number of variables in the ring')
            self._encoding_format.variables = <int *> malloc(sizeof(int) * (len(encoded_vars) + 1))
            for i in range(len(encoded_vars)):
                self._encoding_format.variables[i] = encoded_vars[i]
            self._encoding_format.variables[len(encoded_vars)] = 0
            self._encoding_format.words = len(encoded_vars) + 1
        else:
            self._encoding_format.words = 0
            self._encoding_format.variables = NULL

        if order == 'degrevlex':
            fmpz_mpoly_ctx_init(self._ctx, n, ORD_DEGREVLEX)
        elif order == 'deglex':
            fmpz_mpoly_ctx_init(self._ctx, n, ORD_DEGLEX)
        elif order == 'lex':
            fmpz_mpoly_ctx_init(self._ctx, n, ORD_LEX)
        else:
            raise NotImplementedError("FLINT orderings must be degrevlex, deglex or lex")

        self._cnames = <const char **>malloc(sizeof(char *) * n)
        self._bnames = []
        for i in range(n):
            self._bnames.append(str_to_bytes(names[i]))
            self._cnames[i] = self._bnames[i]

        MPolynomialRing_base.__init__(self, base_ring, n, names, order)

    def __dealloc__(self):
        free(<void *>self._cnames)

    cpdef _coerce_map_from_(self, other):
        """
        Return True if and only if there exists a coercion map from
        ``other`` to ``self``.

        TESTS::

            sage: R.<x,y> = QQ[]
            sage: type(R)
            <type 'sage.rings.polynomial.multi_polynomial_libsingular.MPolynomialRing_libsingular'>
            sage: R.has_coerce_map_from(ZZ['t'])
            False
            sage: R.coerce_map_from(ZZ['x'])
            Coercion map:
              From: Univariate Polynomial Ring in x over Integer Ring
              To:   Multivariate Polynomial Ring in x, y over Rational Field

        """
        if other == ZZ:
            return True
        else:
            return False
#        base_ring = self.base_ring()
#        if other is base_ring:
#            # Because this parent class is a Cython class, the method
#            # UnitalAlgebras.ParentMethods.__init_extra__(), which normally
#            # registers the coercion map from the base ring, is called only
#            # when inheriting from this class in Python (cf. Trac #26958).
#            return self._coerce_map_from_base_ring()
#        f = self._coerce_map_via([base_ring], other)
#        if f is not None:
#            return f
#
#        if isinstance(other, MPolynomialRing_libsingular):
#            if self is other:
#                return True
#            n = other.ngens()
#            if(other.base_ring is base_ring and self.ngens() >= n and
#               self.variable_names()[:n] == other.variable_names()):
#                return True
#            elif base_ring.has_coerce_map_from(other._mpoly_base_ring(self.variable_names())):
#                return True
#        elif isinstance(other, MPolynomialRing_polydict):
#            if self == other:
#                return True
#            elif other.ngens() == 0:
#                return True
#            elif base_ring.has_coerce_map_from(other._mpoly_base_ring(self.variable_names())):
#                return True
#        elif is_PolynomialRing(other):
#            if base_ring.has_coerce_map_from(other._mpoly_base_ring(self.variable_names())):
#                return True

    Element = MPolynomial_flint

    def _element_constructor_(self, element, check=True):
        """
        Construct a new element in this polynomial ring by converting
        ``element`` into ``self`` if possible.

        INPUT:

        - ``element`` -- several types are supported, see below

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]

        We can coerce elements of self to self::

            sage: P._coerce_(x*y + 1/2)
            x*y + 1/2

        We can coerce elements for a ring with the same algebraic properties::

            sage: from sage.rings.polynomial.multi_polynomial_libsingular import MPolynomialRing_libsingular
            sage: R.<x,y,z> = MPolynomialRing_libsingular(QQ,3)
            sage: P == R
            True

            sage: P is R
            False

            sage: P._coerce_(x*y + 1)
            x*y + 1

        We can coerce base ring elements::

            sage: P._coerce_(3/2)
            3/2

        and all kinds of integers::

            sage: P._coerce_(ZZ(1))
            1

            sage: P._coerce_(int(1))
            1

            sage: k.<a> = GF(2^8)
            sage: P.<x,y> = PolynomialRing(k,2)
            sage: P._coerce_(a)
            (a)

            sage: z = QQ['z'].0
            sage: K.<s> = NumberField(z^2 - 2)
            sage: P.<x,y> = PolynomialRing(K, 2)
            sage: P._coerce_(1/2*s)
            (1/2*s)

        TESTS::

            sage: P.<x,y> = PolynomialRing(GF(127))
            sage: P("111111111111111111111111111111111111111111111111111111111")
            21
            sage: P.<x,y> = PolynomialRing(QQ)
            sage: P("111111111111111111111111111111111111111111111111111111111")
            111111111111111111111111111111111111111111111111111111111
            sage: P("31367566080")
            31367566080

        Check if :trac:`7582` is fixed::

            sage: R.<x,y,z> = PolynomialRing(CyclotomicField(2),3)
            sage: R.coerce(1)
            1

        Check if :trac:`6160` is fixed::

            sage: x=var('x')
            sage: K.<j> = NumberField(x-1728)
            sage: R.<b,c> = K[]
            sage: R.coerce(1)
            1

        Check if coercion from zero variable polynomial rings work
        (:trac:`7951`)::

            sage: P = PolynomialRing(QQ,0,'')
            sage: R.<x,y> = QQ[]
            sage: P(5)*x
            5*x
            sage: P = PolynomialRing(ZZ,0,'')
            sage: R.<x,y> = GF(127)[]
            sage: R.coerce(P(5))
            5

        Conversion from strings::

            sage: P.<x,y,z> = QQ[]
            sage: P('x+y + 1/4')
            x + y + 1/4

        Coercion from SINGULAR elements::

            sage: P._singular_()
            polynomial ring, over a field, global ordering
            //   coefficients: QQ
            //   number of vars : 3
            //        block   1 : ordering dp
            //                  : names    x y z
            //        block   2 : ordering C

            sage: P._singular_().set_ring()
            sage: P(singular('x + 3/4'))
            x + 3/4

        Coercion from symbolic variables::

            sage: R = QQ['x,y,z']
            sage: var('x')
            x
            sage: R(x)
            x

        Coercion from 'similar' rings, which maps by index::

            sage: P.<x,y,z> = QQ[]
            sage: R.<a,b,c> = ZZ[]
            sage: P(a)
            x

        ::

            sage: P.<x,y> = QQ[]
            sage: R.<a,b,c> = QQ[]
            sage: R(x)
            a

        Coercion from PARI objects::

            sage: P.<x,y,z> = QQ[]
            sage: P(pari('x^2 + y'))
            x^2 + y
            sage: P(pari('x*y'))
            x*y

        Coercion from boolean polynomials, also by index::

            sage: B.<x,y,z> = BooleanPolynomialRing(3)
            sage: P.<x,y,z> = QQ[]
            sage: P(B.gen(0))
            x

        If everything else fails, we try to coerce to the base ring::

            sage: R.<x,y,z> = GF(3)[]
            sage: R(1/2)
            -1

        Finally, conversions from other polynomial rings which are not
        coercions are provided. Variables are mapped as follows. Say,
        we are mapping an element from `P` to `Q` (this ring). If the
        variables of `P` are a subset of `Q`, we perform a name
        preserving conversion::

            sage: P.<y_2, y_1, z_3, z_2, z_1> = GF(3)[]
            sage: Q = GF(3)['y_4', 'y_3', 'y_2', 'y_1', 'z_5', 'z_4', 'z_3', 'z_2', 'z_1']
            sage: Q(y_1*z_2^2*z_1)
            y_1*z_2^2*z_1

        Otherwise, if `P` has less than or equal the number of
        variables as `Q`, we perform a conversion by index::

            sage: P.<a,b,c> = GF(2)[]
            sage: Q = GF(2)['c','b','d','e']
            sage: f = Q.convert_map_from(P)
            sage: f(a), f(b), f(c)
            (c, b, d)

        ::

            sage: P.<a,b,c> = GF(2)[]
            sage: Q = GF(2)['c','b','d']
            sage: f = Q.convert_map_from(P)
            sage: f(a),f(b),f(c)
            (c, b, d)

        In all other cases, we fail::

            sage: P.<a,b,c,f> = GF(2)[]
            sage: Q = GF(2)['c','d','e']
            sage: f = Q.convert_map_from(P)
            sage: f(a)
            Traceback (most recent call last):
            ...
            TypeError: Could not find a mapping of the passed element to this ring.

        Coerce in a polydict where a coefficient reduces to 0 but isn't 0. ::

            sage: R.<x,y> = QQ[]; S.<xx,yy> = GF(5)[]; S( (5*x*y + x + 17*y)._mpoly_dict_recursive() )
            xx + 2*yy

        Coerce in a polynomial one of whose coefficients reduces to 0. ::

            sage: R.<x,y> = QQ[]; S.<xx,yy> = GF(5)[]; S(5*x*y + x + 17*y)
            xx + 2*yy

        Some other examples that illustrate the same coercion idea::

            sage: R.<x,y> = ZZ[]
            sage: S.<xx,yy> = GF(25,'a')[]
            sage: S(5*x*y + x + 17*y)
            xx + 2*yy

            sage: S.<xx,yy> = Integers(5)[]
            sage: S(5*x*y + x + 17*y)
            xx + 2*yy

        See :trac:`5292`::

            sage: R.<x> = QQ[]; S.<q,t> = QQ[]; F = FractionField(S)
            sage: x in S
            False
            sage: x in F
            False

        Check if :trac:`8228` is fixed::

            sage: P.<x,y> = Zmod(10)[]; P(0)
            0
            sage: P.<x,y> = Zmod(2^10)[]; P(0)
            0

        And :trac:`7597` is fixed if this does not segfault::

            sage: F2 = GF(2)
            sage: F.<x> = GF(2^8)
            sage: R4.<a,b> = PolynomialRing(F)
            sage: R.<u,v> = PolynomialRing(F2)
            sage: P = a
            sage: (P(0,0).polynomial()[0])*u
            0
            sage: P(a,b)
            a

        Check that :trac:`15746` is fixed::

            sage: R.<x,y> = GF(7)[]
            sage: R(2^31)
            2

        Check that :trac:`17964` is fixed::

            sage: K.<a> = QuadraticField(17)
            sage: Q.<x,y> = K[]
            sage: f = (-3*a)*y + (5*a)
            sage: p = K.primes_above(5)[0]
            sage: R = K.residue_field(p)
            sage: S = R['x','y']
            sage: S(f)
            (2*abar)*y

        """

        cdef MPolynomial_flint p

        base_ring = self.base_ring()

        if isinstance(element, MPolynomial_flint):
            n = (<MPolynomial_flint>element)._parent.ngens()
            if element.parent() is self:
                return element
#            elif(base_ring is element.base_ring() and
#                 self.ngens() >= n and
#                 self.variable_names()[:n] == (<MPolynomial_libsingular>element)._parent.variable_names()):
#                if self.term_order() == (<MPolynomial_libsingular>element)._parent.term_order():
#                    _p = prCopyR_NoSort((<MPolynomial_libsingular>element)._poly,
#                                        (<MPolynomial_libsingular>element)._parent_ring,
#                                        _ring)
#                else:
#                    _p = prCopyR((<MPolynomial_libsingular>element)._poly,
#                                 (<MPolynomial_libsingular>element)._parent_ring, _ring)
#                return new_MP(self, _p)
#            elif base_ring.has_coerce_map_from(element.parent()._mpoly_base_ring(self.variable_names())):
#                return self(element._mpoly_dict_recursive(self.variable_names(), base_ring))
#
#        elif isinstance(element, MPolynomial_polydict):
#            if element.parent() == self:
#                bucket = sBucketCreate(_ring)
#                try:
#                    for (m,c) in element.element().dict().iteritems():
#                        mon = p_Init(_ring)
#                        p_SetCoeff(mon, sa2si(c, _ring), _ring)
#                        for pos in m.nonzero_positions():
#                            overflow_check(m[pos], _ring)
#                            p_SetExp(mon, pos+1, m[pos], _ring)
#                        p_Setm(mon, _ring)
#                        #we can use "_m" because we're merging a monomial and
#                        #"Merge" because this monomial is different from the rest
#                        sBucket_Merge_m(bucket, mon)
#                    e=0
#                    #we can use "Merge" because the monomials are distinct
#                    sBucketClearMerge(bucket, &_p, &e)
#                    sBucketDestroy(&bucket)
#                except:
#                     sBucketDeleteAndDestroy(&bucket)
#                     raise
#                return new_MP(self, _p)
#            elif element.parent().ngens() == 0:
#                # zero variable polynomials
#                _p = p_NSet(sa2si(base_ring(element[tuple()]), _ring),
#                        _ring)
#                return new_MP(self, _p)
#            elif base_ring.has_coerce_map_from(element.parent()._mpoly_base_ring(self.variable_names())):
#                return self(element._mpoly_dict_recursive(self.variable_names(), base_ring))
#
#        elif isinstance(element, polynomial_element.Polynomial):
#            if base_ring.has_coerce_map_from(element.parent()._mpoly_base_ring(self.variable_names())):
#                return self(element._mpoly_dict_recursive(self.variable_names(), base_ring))
#
#        if isinstance(element, (SingularElement, cypari2.gen.Gen)):
#            element = str(element)
#
#        if isinstance(element, MPolynomial_libsingular) and element.parent() is not self and element.parent() != self:
#            variable_names_s = element.parent().variable_names()
#            variable_names_t = self.variable_names()
#
#            if set(variable_names_s).issubset(variable_names_t):
#                for v in variable_names_s:
#                    ind_map.append(variable_names_t.index(v)+1)
#            else:
#                ind_map = [i+1 for i in range(_ring.N)]
#
#            if element.parent().ngens() <= self.ngens():
#                # Map the variables by indices
#                _p = p_ISet(0, _ring)
#                Element = <MPolynomial_libsingular>element
#                El_poly = Element._poly
#                El_parent = Element._parent
#                El_ring = Element._parent_ring
#                El_base = El_parent._base
#
#                #this loop needs improvement
#                while El_poly:
#                    c = si2sa(p_GetCoeff(El_poly, El_ring), El_ring, El_base)
#                    if check:
#                        try:
#                            c = base_ring(c)
#                        except TypeError:
#                            p_Delete(&_p, _ring)
#                            raise
#                    if c:
#                        mon = p_Init(_ring)
#                        p_SetCoeff(mon, sa2si(c, _ring), _ring)
#                        for j from 1 <= j <= El_ring.N:
#                            e = p_GetExp(El_poly, j, El_ring)
#                            if e:
#                                p_SetExp(mon, ind_map[j-1], e, _ring)
#                        p_Setm(mon, _ring)
#                        _p = p_Add_q(_p, mon, _ring)
#                    El_poly = pNext(El_poly)
#                return new_MP(self, _p)
#
#        if isinstance(element, MPolynomial_polydict):
#            variable_names_s = element.parent().variable_names()
#            variable_names_t = self.variable_names()
#
#            if set(variable_names_s).issubset(variable_names_t):
#                for v in variable_names_s:
#                    ind_map.append(variable_names_t.index(v)+1)
#            else:
#                ind_map = [i+1 for i in range(_ring.N)]
#
#            if element.parent().ngens() <= self.ngens():
#                bucket = sBucketCreate(_ring)
#                try:
#                    for (m,c) in element.element().dict().iteritems():
#                        if check:
#                            c = base_ring(c)
#                        if not c:
#                            continue
#                        mon = p_Init(_ring)
#                        p_SetCoeff(mon, sa2si(c , _ring), _ring)
#                        for pos in m.nonzero_positions():
#                            overflow_check(m[pos], _ring)
#                            p_SetExp(mon, ind_map[pos], m[pos], _ring)
#                        p_Setm(mon, _ring)
#                        sBucket_Merge_m(bucket, mon)
#                    e=0
#                    sBucketClearMerge(bucket, &_p, &e)
#                    sBucketDestroy(&bucket)
#                except TypeError:
#                    sBucketDeleteAndDestroy(&bucket)
#                    raise
#                return new_MP(self, _p)
#
#        from sage.rings.polynomial.pbori import BooleanPolynomial
#        if isinstance(element, BooleanPolynomial):
#            if element.constant():
#                if element:
#                    return self._one_element
#                else:
#                    return self._zero_element
#
#            variable_names_s = set(element.parent().variable_names())
#            variable_names_t = self.variable_names()
#
#            if variable_names_s.issubset(variable_names_t):
#                return eval(str(element),self.gens_dict(copy=False))
#
#            elif element.parent().ngens() <= self.ngens():
#                Q = element.parent()
#                gens_map = dict(zip(Q.variable_names(),self.gens()[:Q.ngens()]))
#                return eval(str(element),gens_map)

        if isinstance(element, str):
            # let python do the parsing
            d = self.gens_dict()
            if self.base_ring().gen() != 1:
                d[str(self.base_ring().gen())]=self.base_ring().gen()
            try:
                if '/' in element:
                    element = sage_eval(element,d)
                else:
                    element = element.replace("^","**")
                    element = eval(element, d, {})
            except (SyntaxError, NameError):
                raise TypeError("Could not find a mapping of the passed element to this ring.")

            # we need to do this, to make sure that we actually get an
            # element in self.
            return self._coerce_c(element)

#        if isinstance(element, dict):
#            if len(element)==0:
#                _p = p_ISet(0, _ring)
#            else:
#                bucket = sBucketCreate(_ring)
#                try:
#                    for (m,c) in element.iteritems():
#                        if check:
#                            c = base_ring(c)
#                        if not c:
#                            continue
#                        mon = p_Init(_ring)
#                        p_SetCoeff(mon, sa2si(c , _ring), _ring)
#                        if len(m) != self.ngens():
#                            raise TypeError("tuple key must have same length as ngens")
#                        for pos from 0 <= pos < len(m):
#                            if m[pos]:
#                                overflow_check(m[pos], _ring)
#                                p_SetExp(mon, pos+1, m[pos], _ring)
#                        p_Setm(mon, _ring)
#                        sBucket_Merge_m(bucket, mon)
#                    e=0
#                    sBucketClearMerge(bucket, &_p, &e)
#                    sBucketDestroy(&bucket)
#                except TypeError:
#                    sBucketDeleteAndDestroy(&bucket)
#                    raise
#            return new_MP(self, _p)

        try: #if hasattr(element,'_polynomial_'):
            # SymbolicVariable
            return element._polynomial_(self)
        except AttributeError:
            pass

#        if is_Macaulay2Element(element):
#            return self(element.external_string())
#        try:
#            return self(str(element))
#        except TypeError:
#            pass

        try:
            # now try calling the base ring's __call__ methods
            element = self.base_ring()(element)
            elementint = Integer(element)

            p = MPolynomial_flint.__new__(MPolynomial_flint)
            fmpz_mpoly_init(p._poly, self._ctx)
            p._parent = self

            fmpz_mpoly_set_si(p._poly, mpz_get_si(elementint.value), self._ctx)

            return p
        except (TypeError, ValueError):
            raise TypeError("Could not find a mapping of the passed element to this ring.")

    def gen(self, int n=0):
        """
        Returns the ``n``-th generator of this multivariate polynomial
        ring.

        INPUT:

        - ``n`` -- an integer ``>= 0``

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: P.gen(),P.gen(1)
            (x, y)

            sage: P = PolynomialRing(GF(127),1000,'x')
            sage: P.gen(500)
            x500

            sage: P.<SAGE,SINGULAR> = QQ[] # weird names
            sage: P.gen(1)
            SINGULAR
        """

        if n < 0 or n >= self.__ngens:
            raise ValueError("Generator not defined.")

        cdef MPolynomial_flint p = MPolynomial_flint.__new__(MPolynomial_flint)
        p._parent = self

        fmpz_mpoly_gen(p._poly, n, self._ctx)

        return p

    def addmul_multi(self, terms, verbose=False):
        """
        Accepts a list of terms, each of which is a list of factors, and constructs the sum
        of the products.  The factors can be either polynomials or rational functions in
        the polynomial ring's fraction field.
        """

        cdef MPolynomial_flint p = MPolynomial_flint.__new__(MPolynomial_flint)
        p._parent = self

        cdef encode_to_file_struct state

        if verbose: print("fmpz_mpoly_addmul_multi entry", len(terms), file=sys.stderr)

        # Construct a list of polynomials that we will refer to by their index number.
        # Factoring each polynomial seems a bit too time consuming, so I just ensure
        # that they're relatively prime.

        # Each list element is either a polynomial or a list of indices indicating that
        # that polynomial has been factored.

        polys = list()

        # Add a polynomial to the polys list (if needed), ensuring that polys is maintained
        # relatively prime and return a list of indices into polys that when multiplied
        # together form the new polynomial.  Returns a list of indices whose polynomials,
        # when multiplied together, form the input polynomial.

        def add_to_polys(newpoly, startpoly=0):
            assert newpoly.parent() is self
            if verbose: print("add_to_polys: adding polynomial of length", len(newpoly))
            result = []
            factorization = newpoly.factor()
            for poly, power in factorization:
                if poly not in polys:
                    polys.append(poly)
                result += [polys.index(poly)] * power
            return (result, Integer(factorization.unit()))

        # Lists of Counters, one pair for each term, mapping polynomial indices
        # to counts of polynomials in each numerator and denominator
        n = list()
        d = list()

        # Lists of integer coefficients, one pair for each term
        nc = list()
        dc = list()

        # build a list of lists of indices into polys for numerators and denominators
        # next_len = 10
        for term in terms:
            n.append(Counter())
            d.append(Counter())
            nc.append(Integer(1))
            dc.append(Integer(1))
            if verbose: print("addmul_multi: processing term", len(n))
            for factor in term:
               if is_FractionFieldElement(factor) or factor.parent() is self.base_ring().fraction_field():
                   assert(factor.numerator().parent() is self or factor.numerator().parent() is self.base_ring())
                   assert(factor.denominator().parent() is self or factor.denominator().parent() is self.base_ring())
                   l,u = add_to_polys(self(factor.numerator()))
                   n[-1].update(l)
                   nc[-1] *= u
                   # if verbose: print("added polynomial to numerator", len(n), n[-1], nc[-1])
                   l,u = add_to_polys(self(factor.denominator()))
                   d[-1].update(l)
                   dc[-1] *= u
                   #if verbose: print("added polynomial to denominator", len(d), d[-1], dc[-1])
               else:
                   if not (factor.parent() is self or factor.parent() is self.base_ring()):
                       print("factor", factor, factor.parent())
                   assert(factor.parent() is self or factor.parent() is self.base_ring())
                   l,u = add_to_polys(self(factor))
                   n[-1].update(l)
                   nc[-1] *= u
                   #if verbose: print("added polynomial to numerator", len(n), n[-1], nc[-1])
               #if verbose and len(polys) >= next_len:
               #    next_len += 10
               #    print("len(polys) =", len(polys))
            common = n[-1] & d[-1]
            if len(common) > 0:
               n[-1].subtract(common)
               d[-1].subtract(common)

        def flatten(t):
            return [item for sublist in t for item in sublist]

        # construct the LCM of all the denominators, in the form of a list of indices

        def common_elements(lists):
            common_set = set(itertools.chain.from_iterable(lists))
            return flatten([[e]*max(map(lambda l: l[e], lists)) for e in common_set])

        #print("fmpz_mpoly_addmul_multi term_lengths =", term_lengths, file=sys.stderr)

        lcm = common_elements(d)
        lcmc = reduce(lambda a,b: a.lcm(b), dc)

        if verbose: print("fmpz_mpoly_addmul_multi lcm =", lcm, "lcmc =", lcmc, file=sys.stderr)

        #print("fmpz_mpoly_addmul_multi n =", n, file=sys.stderr)
        #print("fmpz_mpoly_addmul_multi d =", d, file=sys.stderr)

        if verbose:
            sys.modules['__main__'].addmul_lcm = lcm
            sys.modules['__main__'].polys = polys
            #raise Exception("addmul")

        if verbose: print("fmpz_mpoly_addmul_multi array building", file=sys.stderr)

        cdef slong * iptr = <slong *>malloc(sizeof(slong) * len(terms))

        num_polys = 0
        for t in range(len(terms)):
           term_len = sum(n[t].values())
           lcm2 = Counter(lcm)
           lcm2.subtract(d[t])
           term_len += sum(lcm2.values())
           if nc[t] * lcmc != dc[t]:
               term_len += 1
           iptr[t] = term_len
           num_polys += term_len

        if verbose: print("fmpz_mpoly_addmul_multi num_polys =", num_polys, file=sys.stderr)

        cdef const fmpz_mpoly_struct ** fptr = <const fmpz_mpoly_struct **>malloc(sizeof(fmpz_mpoly_struct *) * num_polys)
        cdef MPolynomial_flint ni
        constants = []
        k = 0
        maxdeg = 0
        maxlen = 0
        vardeg = {var:0 for var in self.gens()}
        try:
            for t in range(len(terms)):
                thisdeg = 0
                thislen = 1
                thisvardeg = {var:0 for var in self.gens()}
                for i in n[t].elements():
                    thisdeg += polys[i].degree()
                    thislen *= len(polys[i])
                    for var in self.gens():
                        thisvardeg[var] += polys[i].degree(var)
                    ni = polys[i]
                    fptr[k] = <const fmpz_mpoly_struct *>ni._poly
                    k += 1
                lcm2 = Counter(lcm)
                lcm2.subtract(d[t])
                for i in lcm2.elements():
                    thisdeg += polys[i].degree()
                    thislen *= len(polys[i])
                    for var in self.gens():
                        thisvardeg[var] += polys[i].degree(var)
                    ni = polys[i]
                    fptr[k] = <const fmpz_mpoly_struct *>ni._poly
                    k += 1
                if nc[t] * lcmc != dc[t]:
                    constant = Integer(nc[t] * lcmc / dc[t])
                    #if verbose: print("constant", constant)
                    if constant not in constants:
                        constants.append(self(constant))
                    ni = constants[constants.index(constant)]
                    #if verbose: print("constant ni", ni)
                    fptr[k] = <const fmpz_mpoly_struct *>ni._poly
                    k += 1
                if thisdeg > maxdeg:
                    maxdeg = thisdeg
                maxlen += thislen
                for var in self.gens():
                    if thisvardeg[var] > vardeg[var]:
                        vardeg[var] = thisvardeg[var]
        except Exception as ex:
            # make sure we see an exception instead of letting it get caught
            print(ex, file=sys.stdout)
            raise

        if verbose: print("fmpz_mpoly_addmul_multi maxlen =", max(len(p) for p in polys), file=sys.stderr)
        if verbose: print("fmpz_mpoly_addmul_multi maxlen result =", maxlen, file=sys.stderr)
        if verbose: print("fmpz_mpoly_addmul_multi maxdeg =", maxdeg, file=sys.stderr)
        if verbose: print("fmpz_mpoly_addmul_multi vardeg =", vardeg, file=sys.stderr)
        if verbose: print("fmpz_mpoly_addmul_multi len(terms) =", len(terms), file=sys.stderr)

        cdef slong len_terms = len(terms)

        # XXX we use 'verbose' to figure out which polynomials we need to dump to disk
        # XXX obviously need a better way
        #if verbose:
        if False:
            filename = "bigflint.out.gz"
            state.format = & self._encoding_format
            open_file_for_encoding(& state, filename)

            with nogil:
                fmpz_mpoly_addmul_multi_threaded_abstract(<void *> &state, fptr, iptr, len_terms, self._ctx, encode_to_file_returning_status)

            close_file_for_encoding(& state)
        else:
            fmpz_mpoly_addmul_multi_threaded(p._poly, fptr, iptr, len(terms), self._ctx)

        #if verbose: raise Exception("fmpz_mpoly_addmul_multi")
        if verbose: print("fmpz_mpoly_addmul_multi done len(result) =", len(p), file=sys.stderr)
        if len(lcm) == 0:
            return p
        else:
            denom = prod(polys[i] for i in lcm)
            return p/denom

    # It is required in cython to redefine __hash__ when __richcmp__ is
    # overloaded. Also just writing
    #         __hash__ = CategoryObject.__hash__
    # doesn't work.
    def __hash__(self):
        """
        Return a hash for this ring, that is, a hash of the string
        representation of this polynomial ring.

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: hash(P)      # somewhat random output
            967902441410893180 # 64-bit
            -1767675994        # 32-bit
        """
        return CategoryObject.__hash__(self)

    def __reduce__(self):
        """
        Serializes self.

        EXAMPLES::

            sage: P.<x,y,z> = PolynomialRing(QQ, order='degrevlex')
            sage: P == loads(dumps(P))
            True

            sage: P.<x,y,z> = PolynomialRing(ZZ, order='degrevlex')
            sage: P == loads(dumps(P))
            True

            sage: P = PolynomialRing(GF(127), names='abc')
            sage: P == loads(dumps(P))
            True

            sage: P = PolynomialRing(GF(2^8,'F'), names='abc')
            sage: P == loads(dumps(P))
            True

            sage: P = PolynomialRing(GF(2^16,'B'), names='abc')
            sage: P == loads(dumps(P))
            True
            sage: z = QQ['z'].0
            sage: P = PolynomialRing(NumberField(z^2 + 3,'B'), names='abc')
            sage: P == loads(dumps(P))
            True
        """
        return unpickle_MPolynomialRing_flint, \
            (self.base_ring(), self.variable_names(), self.term_order())

def unpickle_MPolynomialRing_flint(base_ring, names, term_order):
    """
    inverse function for ``MPolynomialRing_flint.__reduce__``

    EXAMPLES::

        sage: P.<x,y> = PolynomialRing(QQ)
        sage: loads(dumps(P)) is P # indirect doctest
        True
    """
    from sage.rings.polynomial.polynomial_ring_constructor import _multi_variate
    # If libsingular would be replaced by a different implementation in future
    # sage version, the unpickled ring will belong the new implementation.
    return _multi_variate(base_ring, tuple(names), None, term_order, "FLINT")


cdef class MPolynomial_flint(MPolynomial):
    """
    A multivariate polynomial implemented using FLINT.
    """

    def __cinit__(self):
        fmpz_mpoly_init(self._poly, NULL)

    def __init__(self, MPolynomialRing_flint parent):
        """
        Construct a zero element in parent.

        EXAMPLES::

            sage: from sage.rings.polynomial.multi_polynomial_libsingular import MPolynomial_libsingular
            sage: P = PolynomialRing(GF(32003),3,'x')
            sage: MPolynomial_libsingular(P)
            0
        """
        # fmpz_mpoly_init(self._poly, parent._ctx)
        
        self._parent = parent
        self._ctx = & parent._ctx[0]

    def __dealloc__(self):
        # WARNING: the Cython class self._parent is now no longer accessible!
        fmpz_mpoly_clear(self._poly, self._ctx)

    cpdef _richcmp_(left, right, int op):
        """
        Compare left and right.

        EXAMPLES::

            sage: P.<x,y,z> = PolynomialRing(QQ,order='degrevlex')
            sage: x == x
            True

            sage: x > y
            True
            sage: y^2 > x
            True

            sage: (2/3*x^2 + 1/2*y + 3) > (2/3*x^2 + 1/4*y + 10)
            True

        TESTS::

            sage: P.<x,y,z> = PolynomialRing(QQ, order='degrevlex')
            sage: x > P(0)
            True

            sage: P(0) == P(0)
            True

            sage: P(0) < P(1)
            True

            sage: x > P(1)
            True

            sage: 1/2*x < 3/4*x
            True

            sage: (x+1) > x
            True

            sage: f = 3/4*x^2*y + 1/2*x + 2/7
            sage: f > f
            False
            sage: f < f
            False
            sage: f == f
            True

            sage: P.<x,y,z> = PolynomialRing(GF(127), order='degrevlex')
            sage: (66*x^2 + 23) > (66*x^2 + 2)
            True
        """
        assert (<MPolynomial_flint>left)._parent == (<MPolynomial_flint>right)._parent

        return rich_to_bool(op, fmpz_mpoly_cmp((<MPolynomial_flint>left)._poly, (<MPolynomial_flint>right)._poly, (<MPolynomialRing_flint>left._parent)._ctx))

    def number_of_terms(self):
        """
        Return the number of non-zero coefficients of this polynomial.

        This is also called weight, :meth:`hamming_weight` or sparsity.

        EXAMPLES::

            sage: R.<x, y> = CC[]
            sage: f = x^3 - y
            sage: f.number_of_terms()
            2
            sage: R(0).number_of_terms()
            0
            sage: f = (x+y)^100
            sage: f.number_of_terms()
            101

        The method :meth:`hamming_weight` is an alias::

            sage: f.hamming_weight()
            101
        """
        return len(self)

    hamming_weight = number_of_terms

    cpdef _add_(left, right):
        """
        Add left and right.

        EXAMPLES::

            sage: P.<x,y,z>=PolynomialRing(QQ,3)
            sage: 3/2*x + 1/2*y + 1 #indirect doctest
            3/2*x + 1/2*y + 1
        """

        assert (<MPolynomial_flint>left)._parent == (<MPolynomial_flint>right)._parent

        cdef MPolynomial_flint p = MPolynomial_flint.__new__(MPolynomial_flint)
        p._parent = (<MPolynomial_flint>left)._parent

        fmpz_mpoly_add(p._poly, (<MPolynomial_flint>left)._poly, (<MPolynomial_flint>right)._poly, (<MPolynomialRing_flint>left._parent)._ctx)

        return p

    cpdef _sub_(left, right):
        """
        Subtract left and right.

        EXAMPLES::

            sage: P.<x,y,z>=PolynomialRing(QQ,3)
            sage: 3/2*x - 1/2*y - 1 #indirect doctest
            3/2*x - 1/2*y - 1
        """

        assert (<MPolynomial_flint>left)._parent == (<MPolynomial_flint>right)._parent

        cdef MPolynomial_flint p = MPolynomial_flint.__new__(MPolynomial_flint)
        p._parent = (<MPolynomial_flint>left)._parent

        fmpz_mpoly_sub(p._poly, (<MPolynomial_flint>left)._poly, (<MPolynomial_flint>right)._poly, (<MPolynomialRing_flint>left._parent)._ctx)

        return p

    cpdef _lmul_(self, Element left):
        """
        Multiply self with a base ring element.

        EXAMPLES::

            sage: P.<x,y,z>=PolynomialRing(QQ,3)
            sage: 3/2*x # indirect doctest
            3/2*x

        ::

            sage: P.<x,y,z>=PolynomialRing(QQ,3)
            sage: (3/2*x - 1/2*y - 1) * (3/2) # indirect doctest
            9/4*x - 3/4*y - 3/2
        """

        cdef MPolynomial_flint p = MPolynomial_flint.__new__(MPolynomial_flint)
        p._parent = (<MPolynomial_flint>self)._parent

        leftint = Integer(left)

        fmpz_mpoly_scalar_mul_si(p._poly, (<MPolynomial_flint>self)._poly, mpz_get_si(leftint.value), (<MPolynomialRing_flint>self._parent)._ctx)

        return p


    cpdef _mul_(left, right):
        """
        Multiply left and right.

        EXAMPLES::

            sage: P.<x,y,z>=PolynomialRing(QQ,3)
            sage: (3/2*x - 1/2*y - 1) * (3/2*x + 1/2*y + 1) # indirect doctest
            9/4*x^2 - 1/4*y^2 - y - 1

            sage: P.<x,y> = PolynomialRing(QQ,order='lex')
            sage: (x^2^15) * x^2^15
            Traceback (most recent call last):
            ...
            OverflowError: exponent overflow (...)
        """
        # all currently implemented rings are commutative
        assert (<MPolynomial_flint>left)._parent == (<MPolynomial_flint>right)._parent

        cdef MPolynomial_flint p = MPolynomial_flint.__new__(MPolynomial_flint)
        p._parent = (<MPolynomial_flint>left)._parent

        fmpz_mpoly_mul(p._poly, (<MPolynomial_flint>left)._poly, (<MPolynomial_flint>right)._poly, (<MPolynomialRing_flint>left._parent)._ctx)

        return p

    cpdef _div_(left, right_ringelement):
        """
        Divide left by right

        EXAMPLES::

            sage: R.<x,y>=PolynomialRing(QQ,2)
            sage: f = (x + y)/3 # indirect doctest
            sage: f.parent()
            Multivariate Polynomial Ring in x, y over Rational Field

        Note that / is still a constructor for elements of the
        fraction field in all cases as long as both arguments have the
        same parent and right is not constant. ::

            sage: R.<x,y>=PolynomialRing(QQ,2)
            sage: f = x^3 + y
            sage: g = x
            sage: h = f/g; h
            (x^3 + y)/x
            sage: h.parent()
            Fraction Field of Multivariate Polynomial Ring in x, y over Rational Field

        If we divide over `\ZZ` the result is the same as multiplying
        by 1/3 (i.e. base extension). ::

            sage: R.<x,y> = ZZ[]
            sage: f = (x + y)/3
            sage: f.parent()
            Multivariate Polynomial Ring in x, y over Rational Field
            sage: f = (x + y) * 1/3
            sage: f.parent()
            Multivariate Polynomial Ring in x, y over Rational Field

        But we get a true fraction field if the denominator is not in
        the fraction field of the base ring.""

            sage: f = x/y
            sage: f.parent()
            Fraction Field of Multivariate Polynomial Ring in x, y over Integer Ring

        Division will fail for non-integral domains::

            sage: P.<x,y> = Zmod(1024)[]
            sage: x/P(3)
            Traceback (most recent call last):
            ...
            TypeError: self must be an integral domain.

            sage: x/3
            Traceback (most recent call last):
            ...
            TypeError: self must be an integral domain.

        TESTS::

            sage: R.<x,y>=PolynomialRing(QQ,2)
            sage: x/0
            Traceback (most recent call last):
            ...
            ZeroDivisionError: rational division by zero
        """
        assert (<MPolynomial_flint>left)._parent == (<MPolynomial_flint>right_ringelement)._parent

        return (left._parent).fraction_field()(left,right_ringelement)

    cpdef quo_rem(left, right):
        """
        Returns quotient and remainder of self and right.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: f = y*x^2 + x + 1
            sage: f.quo_rem(x)
            (x*y + 1, 1)
            sage: f.quo_rem(y)
            (x^2, x + 1)

            sage: R.<x,y> = ZZ[]
            sage: f = 2*y*x^2 + x + 1
            sage: f.quo_rem(x)
            (2*x*y + 1, 1)
            sage: f.quo_rem(y)
            (2*x^2, x + 1)
            sage: f.quo_rem(3*x)
            (0, 2*x^2*y + x + 1)

        TESTS::

            sage: R.<x,y> = QQ[]
            sage: R(0).quo_rem(R(1))
            (0, 0)
            sage: R(1).quo_rem(R(0))
            Traceback (most recent call last):
            ...
            ZeroDivisionError

        """
        if right.is_zero():
            raise ZeroDivisionError

        assert (<MPolynomial_flint>left)._parent == (<MPolynomial_flint>right)._parent

        cdef MPolynomial_flint q = MPolynomial_flint.__new__(MPolynomial_flint)
        cdef MPolynomial_flint r = MPolynomial_flint.__new__(MPolynomial_flint)
        q._parent = (<MPolynomial_flint>left)._parent
        r._parent = (<MPolynomial_flint>left)._parent

        fmpz_mpoly_divrem(q._poly, r._poly, (<MPolynomial_flint>left)._poly, (<MPolynomial_flint>right)._poly, (<MPolynomialRing_flint>left._parent)._ctx)

        return q, r

    def _repr_(self):
        """
        EXAMPLES::

            sage: R.<x,y>=PolynomialRing(QQ,2)
            sage: f = x^3 + y
            sage: f # indirect doctest
            x^3 + y
        """
        cdef const char * cstr
        cstr = fmpz_mpoly_get_str_pretty(self._poly, (<MPolynomialRing_flint>self._parent)._cnames, (<MPolynomialRing_flint>self._parent)._ctx)
        pstr = char_to_str(cstr)
        free(<void *>cstr)
        return pstr

    def __len__(self):
        cdef fmpz_mpoly_struct A = self._poly[0]
        return Integer(A.length)

    def degree(self, x=None):
        """
        Return the degree of this polynomial.

        INPUT:

        - ``x`` -- (default: ``None``) a generator of the parent ring

        OUTPUT:

        If ``x`` is not given, return the maximum degree of the monomials of
        the polynomial. Note that the degree of a monomial is affected by the
        gradings given to the generators of the parent ring. If ``x`` is given,
        it is (or coercible to) a generator of the parent ring and the output
        is the maximum degree in ``x``. This is not affected by the gradings of
        the generators.

        EXAMPLES::

            sage: R.<x, y> = QQ[]
            sage: f = y^2 - x^9 - x
            sage: f.degree(x)
            9
            sage: f.degree(y)
            2
            sage: (y^10*x - 7*x^2*y^5 + 5*x^3).degree(x)
            3
            sage: (y^10*x - 7*x^2*y^5 + 5*x^3).degree(y)
            10

        The term ordering of the parent ring determines the grading of the
        generators. ::

            sage: T = TermOrder('wdegrevlex', (1,2,3,4))
            sage: R = PolynomialRing(QQ, 'x', 12, order=T+T+T)
            sage: [x.degree() for x in R.gens()]
            [1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4]

        A matrix term ordering determines the grading of the generators by the
        first row of the matrix. ::

            sage: m = matrix(3, [3,2,1,1,1,0,1,0,0])
            sage: m
            [3 2 1]
            [1 1 0]
            [1 0 0]
            sage: R.<x,y,z> = PolynomialRing(QQ, order=TermOrder(m))
            sage: x.degree(), y.degree(), z.degree()
            (3, 2, 1)
            sage: f = x^3*y + x*z^4
            sage: f.degree()
            11

        If the first row contains zero, the grading becomes the standard one. ::

            sage: m = matrix(3, [3,0,1,1,1,0,1,0,0])
            sage: m
            [3 0 1]
            [1 1 0]
            [1 0 0]
            sage: R.<x,y,z> = PolynomialRing(QQ, order=TermOrder(m))
            sage: x.degree(), y.degree(), z.degree()
            (1, 1, 1)
            sage: f = x^3*y + x*z^4
            sage: f.degree()
            5

        To get the degree with the standard grading regardless of the term
        ordering of the parent ring, use ``std_grading=True``. ::

            sage: f.degree(std_grading=True)
            5

        TESTS::

            sage: P.<x, y> = QQ[]
            sage: P(0).degree(x)
            -1
            sage: P(1).degree(x)
            0

        The following example is inspired by :trac:`11652`::

            sage: R.<p,q,t> = ZZ[]
            sage: poly = p + q^2 + t^3
            sage: poly = poly.polynomial(t)[0]
            sage: poly
            q^2 + p

        There is no canonical coercion from ``R`` to the parent of ``poly``, so
        this doesn't work::

            sage: poly.degree(q)
            Traceback (most recent call last):
            ...
            TypeError: argument is not coercible to the parent

        Using a non-canonical coercion does work, but we require this
        to be done explicitly, since it can lead to confusing results
        if done automatically::

            sage: poly.degree(poly.parent()(q))
            2
            sage: poly.degree(poly.parent()(p))
            1
            sage: T.<x,y> = ZZ[]
            sage: poly.degree(poly.parent()(x))  # noncanonical coercions can be confusing
            1

        The argument to degree has to be a generator::

            sage: pp = poly.parent().gen(0)
            sage: poly.degree(pp)
            1
            sage: poly.degree(pp+1)
            Traceback (most recent call last):
            ...
            TypeError: argument is not a generator

        Canonical coercions are used::

            sage: S = ZZ['p,q']
            sage: poly.degree(S.0)
            1
            sage: poly.degree(S.1)
            2
        """
        cdef slong deg
        if x == None:
            deg = fmpz_mpoly_total_degree_si(self._poly, (<MPolynomialRing_flint>self._parent)._ctx)
            return Integer(deg)
        else:
            for i,gen in enumerate(self.parent().gens()):
                if gen == x:
                    deg = fmpz_mpoly_degree_si(self._poly, i, (<MPolynomialRing_flint>self._parent)._ctx)
                    return Integer(deg)
            raise TypeError('argument is not a generator')

    def dict(self):
        """
        Return a dictionary representing self. This dictionary is in
        the same format as the generic MPolynomial: The dictionary
        consists of ``ETuple:coefficient`` pairs.

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: f=2*x*y^3*z^2 + 1/7*x^2 + 2/3
            sage: f.dict()
            {(0, 0, 0): 2/3, (1, 3, 2): 2, (2, 0, 0): 1/7}
        """
        cdef fmpz_mpoly_struct A = self._poly[0]

        n = (<MPolynomialRing_flint>self._parent).ngens()
        cdef ulong *exp = <ulong *>malloc(sizeof(ulong) * n)

        result = {}

        for i in range(A.length):
            fmpz_mpoly_get_term_exp_ui(exp, self._poly, i, (<MPolynomialRing_flint>self._parent)._ctx)
            explist = []
            for j in range(n):
                explist.append(Integer(exp[j]))
            result[ETuple(explist)] = Integer(fmpz_mpoly_get_coeff_si_ui(self._poly, exp, (<MPolynomialRing_flint>self._parent)._ctx))

        free(exp)

        return result

    def __iter__(self):
        """
        Iterate over ``self`` respecting the term order.
        """

        cdef fmpz_mpoly_struct A = self._poly[0]

        n = (<MPolynomialRing_flint>self._parent).ngens()
        cdef ulong *exp = <ulong *>malloc(sizeof(ulong) * n)

        result = {}

        for i in range(A.length):
            fmpz_mpoly_get_term_exp_ui(exp, self._poly, i, (<MPolynomialRing_flint>self._parent)._ctx)
            explist = []
            for j in range(n):
                explist.append(Integer(exp[j]))
            coeff = Integer(fmpz_mpoly_get_coeff_si_ui(self._poly, exp, (<MPolynomialRing_flint>self._parent)._ctx))
            # FIXME - should return a monomial, not an ETuple
            yield (coeff, ETuple(explist))

        free(exp)

    def iterator_exp_coeff(self):
        """
        Iterate over "self" as pairs of (ETuple, coefficient)
        """

        cdef fmpz_mpoly_struct A = self._poly[0]

        n = (<MPolynomialRing_flint>self._parent).ngens()
        cdef ulong *exp = <ulong *>malloc(sizeof(ulong) * n)

        result = {}

        for i in range(A.length):
            fmpz_mpoly_get_term_exp_ui(exp, self._poly, i, (<MPolynomialRing_flint>self._parent)._ctx)
            explist = []
            for j in range(n):
                explist.append(Integer(exp[j]))
            coeff = Integer(fmpz_mpoly_get_coeff_si_ui(self._poly, exp, (<MPolynomialRing_flint>self._parent)._ctx))
            yield (ETuple(explist), coeff)

        free(exp)

    def exponents(self, as_ETuples=True):
        """
        Return the exponents of the monomials appearing in this
        polynomial.

        INPUT:

        - ``as_ETuples`` -- (default: ``True``) if ``True`` returns the
          result as an list of ETuples, otherwise returns a list of tuples

        EXAMPLES::

            sage: R.<a,b,c> = QQ[]
            sage: f = a^3 + b + 2*b^2
            sage: f.exponents()
            [(3, 0, 0), (0, 2, 0), (0, 1, 0)]
            sage: f.exponents(as_ETuples=False)
            [(3, 0, 0), (0, 2, 0), (0, 1, 0)]
        """
        cdef fmpz_mpoly_struct A = self._poly[0]

        n = (<MPolynomialRing_flint>self._parent).ngens()
        cdef ulong *exp = <ulong *>malloc(sizeof(ulong) * n)

        result = []

        for i in range(A.length):
            fmpz_mpoly_get_term_exp_ui(exp, self._poly, i, (<MPolynomialRing_flint>self._parent)._ctx)
            explist = []
            for j in range(n):
                explist.append(Integer(exp[j]))
            if as_ETuples:
                result.append(ETuple(explist))
            else:
                result.append(tuple(explist))

        free(exp)

        return result

    def inverse_of_unit(self):
        """
        Return the inverse of this polynomial if it is a unit.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: x.inverse_of_unit()
            Traceback (most recent call last):
            ...
            ArithmeticError: Element is not a unit.

            sage: R(1/2).inverse_of_unit()
            2
        """


        if fmpz_mpoly_equal_si((<MPolynomial_flint>self)._poly, 1, (<MPolynomialRing_flint>self._parent)._ctx) or \
           fmpz_mpoly_equal_si((<MPolynomial_flint>self)._poly, -1, (<MPolynomialRing_flint>self._parent)._ctx):
            return self
        else:
            raise ArithmeticError("Element is not a unit.")

    def constant_coefficient(self):
        """
        Return the constant coefficient of this multivariate
        polynomial.

        EXAMPLES::

            sage: P.<x, y> = QQ[]
            sage: f = 3*x^2 - 2*y + 7*x^2*y^2 + 5
            sage: f.constant_coefficient()
            5
            sage: f = 3*x^2
            sage: f.constant_coefficient()
            0
        """
        n = (<MPolynomialRing_flint>self._parent).ngens()
        cdef ulong *exp = <ulong *>malloc(sizeof(ulong) * n)

        for i in range(n):
            exp[i] = 0

        cdef slong coeff = fmpz_mpoly_get_coeff_si_ui(self._poly, exp, (<MPolynomialRing_flint>self._parent)._ctx)

        free(exp)

        return Integer(coeff)

    cpdef _floordiv_(self, right):
        """
        Perform division with remainder and return the quotient.

        INPUT:

        - ``right`` - something coercible to an MPolynomial_flint
          in ``self.parent()``

        EXAMPLES::

            sage: R.<x,y,z> = GF(32003)[]
            sage: f = y*x^2 + x + 1
            sage: f//x
            x*y + 1
            sage: f//y
            x^2

            sage: P.<x,y> = ZZ[]
            sage: x//y
            0
            sage: (x+y)//y
            1

            sage: P.<x,y> = QQ[]
            sage: (x+y)//y
            1
            sage: (x)//y
            0

            sage: P.<x,y> = Zmod(1024)[]
            sage: (x+y)//x
            1
            sage: (x+y)//(2*x)
            Traceback (most recent call last):
            ...
            NotImplementedError: Division of multivariate polynomials over non fields by non-monomials not implemented.

        TESTS::

            sage: P.<x,y> = ZZ[]
            sage: p = 3*(-x^8*y^2 - x*y^9 + 6*x^8*y + 17*x^2*y^6 - x^3*y^2)
            sage: q = 7*(x^2 + x*y + y^2 + 1)
            sage: p*q//q == p
            True
            sage: p*q//p == q
            True
        """

        assert (<MPolynomial_flint>self)._parent == (<MPolynomial_flint>right)._parent

        cdef MPolynomial_flint p = MPolynomial_flint.__new__(MPolynomial_flint)
        p._parent = (<MPolynomial_flint>self)._parent

        fmpz_mpoly_div(p._poly, (<MPolynomial_flint>self)._poly, (<MPolynomial_flint>right)._poly, (<MPolynomialRing_flint>self._parent)._ctx)

        return p

#        if _right._poly == NULL:
#            raise ZeroDivisionError
#        elif p_IsOne(_right._poly, r):
#            return self

#        if r.cf.type != n_unknown:
#            if (singular_polynomial_length_bounded(_right._poly, 2) == 1
#                    and r.cf.cfIsOne(p_GetCoeff(_right._poly, r), r.cf)):
#                p = self._poly
#                quo = p_ISet(0,r)
#                while p:
#                    if p_DivisibleBy(_right._poly, p, r):
#                        temp = p_MDivide(p, _right._poly, r)
#                        p_SetCoeff0(temp, n_Copy(p_GetCoeff(p, r), r), r)
#                        quo = p_Add_q(quo, temp, r)
#                    p = pNext(p)
#                return new_MP(parent, quo)
#            if r.cf.type == n_Znm or r.cf.type == n_Zn or r.cf.type == n_Z2m :
#                raise NotImplementedError("Division of multivariate polynomials over non fields by non-monomials not implemented.")

#        count = singular_polynomial_length_bounded(self._poly, 15)

#        # fast in the most common case where the division is exact; returns zero otherwise
#        if count >= 15:  # note that _right._poly must be of shorter length than self._poly for us to care about this call
#            sig_on()
#        quo = p_Divide(p_Copy(self._poly, r), p_Copy(_right._poly, r), r)
#        if count >= 15:
#            sig_off()

#        if quo == NULL:
#            if r.cf.type == n_Z:
#                P = parent.change_ring(QQ)
#                f = (<MPolynomial_libsingular>P(self))._floordiv_(P(right))
#                return parent(sum([c.floor() * m for c, m in f]))
#            else:
#                sig_on()
#                quo = singclap_pdivide(self._poly, _right._poly, r)
#                sig_off()

#        return new_MP(parent, quo)


    @coerce_binop
    def gcd(left, right, algorithm=None, **kwds):
        """
        Return the greatest common divisor of left and right.

        INPUT:

        - ``right`` - polynomial
        - ``algorithm``
          - ``ezgcd`` - EZGCD algorithm
          - ``modular`` - multi-modular algorithm (default)
        - ``**kwds`` - ignored

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: f = (x*y*z)^6 - 1
            sage: g = (x*y*z)^4 - 1
            sage: f.gcd(g)
            x^2*y^2*z^2 - 1
            sage: GCD([x^3 - 3*x + 2, x^4 - 1, x^6 -1])
            x - 1

            sage: R.<x,y> = QQ[]
            sage: f = (x^3 + 2*y^2*x)^2
            sage: g = x^2*y^2
            sage: f.gcd(g)
            x^2

        We compute a gcd over a finite field::

            sage: F.<u> = GF(31^2)
            sage: R.<x,y,z> = F[]
            sage: p = x^3 + (1+u)*y^3 + z^3
            sage: q = p^3 * (x - y + z*u)
            sage: gcd(p,q)
            x^3 + (u + 1)*y^3 + z^3
            sage: gcd(p,q)  # yes, twice -- tests that singular ring is properly set.
            x^3 + (u + 1)*y^3 + z^3

        We compute a gcd over a number field::

            sage: x = polygen(QQ)
            sage: F.<u> = NumberField(x^3 - 2)
            sage: R.<x,y,z> = F[]
            sage: p = x^3 + (1+u)*y^3 + z^3
            sage: q = p^3 * (x - y + z*u)
            sage: gcd(p,q)
            x^3 + (u + 1)*y^3 + z^3

        TESTS::

            sage: Q.<x,y,z> = QQ[]
            sage: P.<x,y,z> = QQ[]
            sage: P(0).gcd(Q(0))
            0
            sage: x.gcd(1)
            1

            sage: k.<a> = GF(9)
            sage: R.<x,y> = PolynomialRing(k)
            sage: f = R.change_ring(GF(3)).gen()
            sage: g = x+y
            sage: g.gcd(f)
            1
            sage: x.gcd(R.change_ring(GF(3)).gen())
            x

            sage: Pol.<x,y,z> = ZZ[]
            sage: p = x*y - 5*y^2 + x*z - z^2 + z
            sage: q = -3*x^2*y^7*z + 2*x*y^6*z^3 + 2*x^2*y^3*z^4 + x^2*y^5 - 7*x*y^5*z
            sage: (21^3*p^2*q).gcd(35^2*p*q^2) == -49*p*q
            True
        """
#        cdef poly *_res
#        cdef ring *_ring = self._parent_ring
#        cdef MPolynomial_libsingular _right = <MPolynomial_libsingular>right
#
#        if _right._poly == NULL:
#            return self
#        elif self._poly == NULL:
#            return right
#        elif p_IsOne(self._poly, _ring):
#            return self
#        elif p_IsOne(_right._poly, _ring):
#            return right
#
#        res = new_MP(self._parent, _res)
#        return res

        assert (<MPolynomial_flint>left)._parent == (<MPolynomial_flint>right)._parent

        cdef MPolynomial_flint p = MPolynomial_flint.__new__(MPolynomial_flint)
        p._parent = (<MPolynomial_flint>left)._parent

        if not fmpz_mpoly_gcd(p._poly, (<MPolynomial_flint>left)._poly, (<MPolynomial_flint>right)._poly, (<MPolynomialRing_flint>left._parent)._ctx):
            raise RuntimeError("GCD failed")

        return p

    def factor(self, proof=None):
        cdef fmpz_mpoly_factor_t f
        fmpz_mpoly_factor_init(f, (<MPolynomialRing_flint>self._parent)._ctx)
        if not fmpz_mpoly_factor(f, self._poly, (<MPolynomialRing_flint>self._parent)._ctx):
            raise RuntimeError("multivariate factorization failed")
        factors = []
        cdef MPolynomial_flint p
        for i in range(f.num):
            p = MPolynomial_flint.__new__(MPolynomial_flint)
            p._parent = self._parent
            fmpz_mpoly_init(p._poly, (<MPolynomialRing_flint>self._parent)._ctx)
            # would like to fmpz_mpoly_swap, but this is an inline function in FLINT
            # and I'm not sure how to import that into Cython
            fmpz_mpoly_set(p._poly, &f.poly[i], (<MPolynomialRing_flint>self._parent)._ctx)
            factors.append((p, f.exp[i]))
        return Factorization(factors, unit=fmpz_get_si(f.constant))

    @coerce_binop
    def lcm(left, right, algorithm=None, **kwds):
        """
        Return the least common multiple of left and right.
        """
#        cdef poly *_res
#        cdef ring *_ring = self._parent_ring
#        cdef MPolynomial_libsingular _right = <MPolynomial_libsingular>right
#
#        if _right._poly == NULL:
#            return self
#        elif self._poly == NULL:
#            return right
#        elif p_IsOne(self._poly, _ring):
#            return self
#        elif p_IsOne(_right._poly, _ring):
#            return right
#
#        res = new_MP(self._parent, _res)
#        return res

        return left*right/left.gcd(right)

    def __reduce__(self):
        """
        Serialize this polynomial.

        EXAMPLES::

            sage: P.<x,y,z> = PolynomialRing(QQ,3, order='degrevlex')
            sage: f = 27/113 * x^2 + y*z + 1/2
            sage: f == loads(dumps(f))
            True

            sage: P = PolynomialRing(GF(127),3,names='abc')
            sage: a,b,c = P.gens()
            sage: f = 57 * a^2*b + 43 * c + 1
            sage: f == loads(dumps(f))
            True

        TESTS:

        Verify that :trac:`9220` is fixed.

            sage: R=QQ['x']
            sage: S=QQ['x','y']
            sage: h=S.0^2
            sage: parent(h(R.0,0))
            Univariate Polynomial Ring in x over Rational Field
        """
        return unpickle_MPolynomial_flint, (self._parent, self.dict())

def unpickle_MPolynomial_flint(MPolynomialRing_flint R, d):
    """
    Deserialize an ``MPolynomial_flint`` object

    INPUT:

    - ``R`` - the base ring
    - ``d`` - a Python dictionary as returned by :meth:`MPolynomial_flint.dict`

    EXAMPLES::

        sage: P.<x,y> = PolynomialRing(QQ)
        sage: loads(dumps(x)) == x # indirect doctest
        True
    """

    cdef MPolynomial_flint p = MPolynomial_flint.__new__(MPolynomial_flint)
    p._parent = R
    cdef int N = R.ngens()
    cdef ulong *exp = <ulong *>malloc(sizeof(ulong) * N)
    cdef int _i, _e

    try:
        for mon,c in d.iteritems():
            for i in range(N):
                exp[i] = 0
            for i,e in mon.sparse_iter():
                _i = i
                if _i >= N:
                    raise TypeError("variable index too big")
                _e = e
                if _e <= 0:
                    raise TypeError("exponent too small")
                exp[_i] = _e

            fmpz_mpoly_push_term_si_ui(p._poly, mpz_get_si(Integer(c).value), exp, R._ctx)

    except:
        free(exp)
        raise

    free(exp)

    fmpz_mpoly_sort_terms(p._poly, R._ctx)
    fmpz_mpoly_combine_like_terms(p._poly, R._ctx)
    fmpz_mpoly_assert_canonical(p._poly, R._ctx)

    return p
