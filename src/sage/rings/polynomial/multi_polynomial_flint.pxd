from sage.libs.flint.types cimport fmpz_mpoly_ctx_t, fmpz_mpoly_t

from sage.rings.polynomial.multi_polynomial_ring_base cimport MPolynomialRing_base
from sage.rings.polynomial.multi_polynomial cimport MPolynomial
from sage.rings.integer cimport Integer
from sage.structure.parent cimport Parent

cdef class MPolynomialRing_flint(MPolynomialRing_base):
    cdef fmpz_mpoly_ctx_t _ctx
    cdef object _bnames
    cdef const char ** _cnames

cdef class MPolynomial_flint(MPolynomial):
    cdef fmpz_mpoly_t _poly
