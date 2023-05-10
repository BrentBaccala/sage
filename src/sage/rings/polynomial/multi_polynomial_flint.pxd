from sage.libs.flint.types cimport fmpz_mpoly_ctx_struct, fmpz_mpoly_ctx_t, fmpz_mpoly_t

from sage.rings.polynomial.multi_polynomial_ring_base cimport MPolynomialRing_base
from sage.rings.polynomial.multi_polynomial cimport MPolynomial
from sage.rings.integer cimport Integer
from sage.structure.parent cimport Parent

ctypedef struct encoding_format:
    # number of 64-bits words used to encode (and length of variables array)
    int words
    # array of integers, either non-zero to encode that many vars as deglex64,
    # or zero to encode the coefficient as sint64
    int * variables

cdef class MPolynomialRing_flint(MPolynomialRing_base):
    cdef fmpz_mpoly_ctx_t _ctx
    cdef object _bnames
    cdef const char ** _cnames
    cdef encoding_format _encoding_format

cdef class MPolynomial_flint(MPolynomial):
    cdef fmpz_mpoly_t _poly
    cdef fmpz_mpoly_ctx_struct * _ctx
    cpdef quo_rem(left, right)
