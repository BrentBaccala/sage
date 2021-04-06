# distutils: depends = flint/flint.h flint/fmpz.h flint/fmpz_poly.h flint/fmpz_mat.h flint/fmpq.h flint/fmpq_poly.h flint/fmpq_mat.h flint/fmpz_mod_poly.h flint/nmod_poly.h flint/fq.h flint/fq_nmod.h flint/ulong_extras.h flint/padic.h flint/padic_poly.h flint/qadic.h flint/fmpz_poly_q.h

"""
Declarations for FLINT types
"""

#*****************************************************************************
#       Copyright (C) 2014 Jeroen Demeyer
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.libs.gmp.types cimport *

# Use these typedefs in lieu of flint's ulong and slong macros
ctypedef mp_limb_t ulong
ctypedef mp_limb_signed_t slong


# flint/flint.h:
cdef extern from "flint_wrap.h":
    ctypedef void* flint_rand_t
    cdef long FLINT_BITS
    cdef long FLINT_D_BITS
    ctypedef ulong flint_bitcnt_t

# flint/fmpz.h:
cdef extern from "flint_wrap.h":
    ctypedef slong fmpz
    ctypedef fmpz fmpz_t[1]

    bint COEFF_IS_MPZ(fmpz)
    mpz_ptr COEFF_TO_PTR(fmpz)

    ctypedef struct fmpz_preinvn_struct:
        mp_ptr dinv
        long n
        mp_bitcnt_t norm

    ctypedef fmpz_preinvn_struct[1] fmpz_preinvn_t

# flint/fmpz_poly.h:
cdef extern from "flint_wrap.h":
    ctypedef struct fmpz_poly_struct:
        fmpz* coeffs
        long alloc
        long length

    ctypedef fmpz_poly_struct fmpz_poly_t[1]

# flint/fmpz_mat.h:
cdef extern from "flint_wrap.h":
    ctypedef struct fmpz_mat_struct:
        pass

    ctypedef fmpz_mat_struct fmpz_mat_t[1]

# flint/fmpq.h:
cdef extern from "flint_wrap.h":
    ctypedef struct fmpq:
        pass

    ctypedef fmpq fmpq_t[1]

# flint/fmpq_poly.h:
cdef extern from "flint_wrap.h":
    ctypedef struct fmpq_poly_struct:
        pass

    ctypedef fmpq_poly_struct fmpq_poly_t[1]

# flint/fmpq_mat.h:
cdef extern from "flint_wrap.h":
    ctypedef struct fmpq_mat_struct:
        pass

    ctypedef fmpq_mat_struct fmpq_mat_t[1]

# flint/fmpz_poly_mat.h:
cdef extern from "flint_wrap.h":
    ctypedef struct fmpz_poly_mat_struct:
        pass

    ctypedef fmpz_poly_mat_struct fmpz_poly_mat_t[1]

# flint/fmpz_mod_poly.h:
cdef extern from "flint_wrap.h":
    ctypedef struct fmpz_mod_poly_struct:
        pass

    ctypedef fmpz_mod_poly_struct fmpz_mod_poly_t[1]

# flint/mpoly.h:
cdef extern from "flint_wrap.h":
    ctypedef struct mpoly_ctx_struct:
        pass

    ctypedef mpoly_ctx_struct mpoly_ctx_t[1]

# flint/fmpz_mpoly.h:
cdef extern from "flint_wrap.h":
    ctypedef struct fmpz_mpoly_struct:
        fmpz * coeffs
        ulong * exps
        slong alloc
        slong length
        flint_bitcnt_t bits

    ctypedef fmpz_mpoly_struct fmpz_mpoly_t[1]

    ctypedef struct fmpz_mpoly_ctx_struct:
        mpoly_ctx_t minfo

    ctypedef fmpz_mpoly_ctx_struct fmpz_mpoly_ctx_t[1]

# flint/nmod_poly.h:
cdef extern from "flint_wrap.h":
    ctypedef struct nmod_t:
        mp_limb_t n
        mp_limb_t ninv
        mp_bitcnt_t norm

    ctypedef struct nmod_poly_struct:
        mp_limb_t *coeffs
        long alloc
        long length
        nmod_t mod

    ctypedef nmod_poly_struct nmod_poly_t[1]

    ctypedef struct nmod_poly_factor_struct:
        nmod_poly_t p
        long *exp
        long num
        long alloc

    ctypedef nmod_poly_factor_struct nmod_poly_factor_t[1]

# flint/fq.h:
cdef extern from "flint_wrap.h":
    ctypedef struct fq_ctx_struct:
        fmpz_mod_poly_t modulus

    ctypedef fq_ctx_struct fq_ctx_t[1]

    ctypedef fmpz_poly_struct fq_struct
    ctypedef fmpz_poly_t fq_t

# flint/fq_nmod.h:
cdef extern from "flint_wrap.h":
    ctypedef struct fq_nmod_ctx_struct:
        nmod_poly_t modulus

    ctypedef fq_nmod_ctx_struct fq_nmod_ctx_t[1]

    ctypedef nmod_poly_struct fq_nmod_struct
    ctypedef nmod_poly_t fq_nmod_t

# flint/ulong_extras.h:
cdef extern from "flint_wrap.h":
    ctypedef struct n_factor_t:
        int num
        unsigned long exp[15]
        unsigned long p[15]

# flint/padic.h:
cdef extern from "flint_wrap.h":
    ctypedef struct padic_struct:
        fmpz u
        long v

    ctypedef padic_struct padic_t[1]

    cdef enum padic_print_mode:
        PADIC_TERSE
        PADIC_SERIES
        PADIC_VAL_UNIT

    ctypedef struct padic_ctx_struct:
        fmpz_t p
        long N
        double pinv
        fmpz* pow
        long min
        long max

    ctypedef padic_ctx_struct padic_ctx_t[1]

    ctypedef struct padic_inv_struct:
        long n
        fmpz *pow
        fmpz *u

    ctypedef padic_inv_struct padic_inv_t[1]

# flint/padic_poly.h:
cdef extern from "flint_wrap.h":
    ctypedef struct padic_poly_struct:
        fmpz *coeffs
        long alloc
        long length
        long val
        long N

    ctypedef padic_poly_struct padic_poly_t[1]

# flint/qadic.h:
cdef extern from "flint_wrap.h":
    ctypedef struct qadic_ctx_struct:
        padic_ctx_struct pctx
        fmpz *a
        long *j
        long len
        char *var

    ctypedef qadic_ctx_struct qadic_ctx_t[1]

    ctypedef padic_poly_struct qadic_struct
    ctypedef padic_poly_t qadic_t

# flint/fmpz_poly_q.h:
cdef extern from "flint_wrap.h":
    ctypedef struct fmpz_poly_q_struct:
        fmpz_poly_struct *num
        fmpz_poly_struct *den

    ctypedef fmpz_poly_q_struct fmpz_poly_q_t[1]
