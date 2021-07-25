# distutils: libraries = flint
# distutils: depends = flint/fmpz_mpoly.h

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport mpz_t
from sage.libs.flint.types cimport *
from sage.libs.flint.mpoly cimport *

# flint/fmpz_mpoly.h
cdef extern from "flint_wrap.h":

    void fmpz_mpoly_factor_init(fmpz_mpoly_factor_t f, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_factor_clear(fmpz_mpoly_factor_t f, const fmpz_mpoly_ctx_t ctx)
    int fmpz_mpoly_factor(fmpz_mpoly_factor_t f, const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
    int fmpz_mpoly_factor_expand(fmpz_mpoly_t A, const fmpz_mpoly_factor_t f, const fmpz_mpoly_ctx_t ctx)
