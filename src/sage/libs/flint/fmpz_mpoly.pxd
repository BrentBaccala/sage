# distutils: libraries = flint
# distutils: depends = flint/fmpz_mpoly.h

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport mpz_t
from sage.libs.flint.types cimport *
from sage.libs.flint.mpoly cimport *

# flint/fmpz_mpoly.h
cdef extern from "flint_wrap.h":

# Context object

    void fmpz_mpoly_ctx_init(fmpz_mpoly_ctx_t, slong, const ordering_t)
    void fmpz_mpoly_ctx_init_rand(fmpz_mpoly_ctx_t, flint_rand_t, slong);
    void fmpz_mpoly_ctx_clear(fmpz_mpoly_ctx_t)

#    void fmpz_mpoly_init(fmpz_mpoly_t, const fmpz_mpoly_ctx_t)
#    void fmpz_mpoly_init2(fmpz_mpoly_t, slong, const fmpz_mpoly_ctx_t)
#    void fmpz_mpoly_init3(fmpz_mpoly_t, slong, flint_bitcnt_t, const fmpz_mpoly_ctx_t)

#    char * fmpz_mpoly_get_str_pretty(const fmpz_mpoly_t A, const char **x, const fmpz_mpoly_ctx_t ctx)

#    void fmpz_mpoly_gen(fmpz_mpoly_t, slong, const fmpz_mpoly_ctx_t)

#    void fmpz_mpoly_add(fmpz_mpoly_t, const fmpz_mpoly_t, const fmpz_mpoly_t, const fmpz_mpoly_ctx_t)


# Memory management

    void fmpz_mpoly_init(fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_init2(fmpz_mpoly_t A, slong alloc, const fmpz_mpoly_ctx_t ctx);

    #void fmpz_mpoly_init3(fmpz_mpoly_t A, slong alloc, flint_bitcnt_t bits, const fmpz_mpoly_ctx_t ctx);

    void _fmpz_mpoly_realloc(fmpz ** Acoeff, ulong ** Aexp, slong * Aalloc, slong len, slong N);

    void fmpz_mpoly_realloc(fmpz_mpoly_t A, slong alloc, const fmpz_mpoly_ctx_t ctx);

    void _fmpz_mpoly_fit_length(fmpz ** Acoeff, ulong ** Aexp, slong * Aalloc, slong len, slong N);

    void fmpz_mpoly_fit_length(fmpz_mpoly_t A, slong len, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_clear(fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx);

# Input/output

    int fmpz_mpoly_set_str_pretty(fmpz_mpoly_t A, const char * str, const char ** x, const fmpz_mpoly_ctx_t ctx);

    #char * _fmpz_mpoly_get_str_pretty(const fmpz * poly, const ulong * exps, slong len, const char ** x, slong bits, const mpoly_ctx_t mctx);

    char * fmpz_mpoly_get_str_pretty(const fmpz_mpoly_t A, const char ** x, const fmpz_mpoly_ctx_t ctx);

    #int _fmpz_mpoly_fprint_pretty(FILE * file, const fmpz * poly, const ulong * exps, slong len, const char ** x_in, flint_bitcnt_t bits, const mpoly_ctx_t mctx);

    int fmpz_mpoly_fprint_pretty(FILE * file, const fmpz_mpoly_t A, const char ** x, const fmpz_mpoly_ctx_t ctx);

#  Basic manipulation

    void fmpz_mpoly_gen(fmpz_mpoly_t poly, slong i, const fmpz_mpoly_ctx_t ctx);

    int fmpz_mpoly_is_gen(const fmpz_mpoly_t poly, slong k, const fmpz_mpoly_ctx_t ctx);

    void _fmpz_mpoly_set(fmpz * poly1, ulong * exps1, const fmpz * poly2, const ulong * exps2, slong n, slong N);

    void fmpz_mpoly_set(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx);

    int _fmpz_mpoly_equal(fmpz * poly1, ulong * exps1, const fmpz * poly2, const ulong * exps2, slong n, slong N);

    int fmpz_mpoly_equal(const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx);

# Constants

    int fmpz_mpoly_is_fmpz(const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_get_fmpz(fmpz_t c, const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_set_fmpz(fmpz_mpoly_t A, const fmpz_t c, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_set_ui(fmpz_mpoly_t A, ulong c, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_set_si(fmpz_mpoly_t A, slong c, const fmpz_mpoly_ctx_t ctx);

    int fmpz_mpoly_equal_fmpz(const fmpz_mpoly_t A, const fmpz_t c, const fmpz_mpoly_ctx_t ctx);

    int fmpz_mpoly_equal_ui(const fmpz_mpoly_t A, ulong c, const fmpz_mpoly_ctx_t ctx);

    int fmpz_mpoly_equal_si(const fmpz_mpoly_t A, slong c, const fmpz_mpoly_ctx_t ctx);

# Degrees

    slong fmpz_mpoly_degree_si(const fmpz_mpoly_t A, slong var, const fmpz_mpoly_ctx_t ctx);

    slong fmpz_mpoly_total_degree_si(const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx);


# Coefficients

    void fmpz_mpoly_get_coeff_fmpz_monomial(fmpz_t c, const fmpz_mpoly_t A, const fmpz_mpoly_t M, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_set_coeff_fmpz_monomial(fmpz_mpoly_t A, const fmpz_t c, const fmpz_mpoly_t M, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_get_coeff_fmpz_fmpz(fmpz_t c, const fmpz_mpoly_t A, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx);

    ulong fmpz_mpoly_get_coeff_ui_fmpz(           const fmpz_mpoly_t A, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx);

    slong fmpz_mpoly_get_coeff_si_fmpz(           const fmpz_mpoly_t A, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_get_coeff_fmpz_ui(fmpz_t c, const fmpz_mpoly_t A, const ulong * exp, const fmpz_mpoly_ctx_t ctx);

    ulong fmpz_mpoly_get_coeff_ui_ui(           const fmpz_mpoly_t A, const ulong * exp, const fmpz_mpoly_ctx_t ctx);

    slong fmpz_mpoly_get_coeff_si_ui(           const fmpz_mpoly_t A, const ulong * exp, const fmpz_mpoly_ctx_t ctx);

    void _fmpz_mpoly_set_coeff_fmpz_fmpz(fmpz_mpoly_t A, const fmpz_t c, const fmpz * exp, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_set_coeff_fmpz_fmpz(fmpz_mpoly_t A, const fmpz_t c, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_set_coeff_ui_fmpz(fmpz_mpoly_t A, const ulong c, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_set_coeff_si_fmpz(fmpz_mpoly_t A, const slong c, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_set_coeff_fmpz_ui(fmpz_mpoly_t A, const fmpz_t c, const ulong * exp, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_set_coeff_ui_ui(fmpz_mpoly_t A, const ulong c, const ulong * exp, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_set_coeff_si_ui(fmpz_mpoly_t A, const slong c, const ulong * exp, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_get_coeff_vars_ui(fmpz_mpoly_t C, const fmpz_mpoly_t A, const slong * vars, const ulong * exps, slong length, const fmpz_mpoly_ctx_t ctx);

# comparison

    int fmpz_mpoly_cmp(const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx);

# container operations

    int fmpz_mpoly_is_canonical(const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_resize(fmpz_mpoly_t A, slong new_length, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_get_term_coeff_fmpz(fmpz_t c, const fmpz_mpoly_t A, slong i, const fmpz_mpoly_ctx_t ctx);

    ulong fmpz_mpoly_get_term_coeff_ui(           const fmpz_mpoly_t A, slong i, const fmpz_mpoly_ctx_t ctx);

    slong fmpz_mpoly_get_term_coeff_si(           const fmpz_mpoly_t A, slong i, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_set_term_coeff_fmpz(fmpz_mpoly_t A, slong i, const fmpz_t c, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_set_term_coeff_ui(fmpz_mpoly_t A, slong i, ulong c, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_set_term_coeff_si(fmpz_mpoly_t A, slong i, slong c, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_get_term_exp_fmpz(fmpz ** exp, const fmpz_mpoly_t A, slong i, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_get_term_exp_ui(ulong * exp, const fmpz_mpoly_t A, slong i, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_get_term_exp_si(slong * exp, const fmpz_mpoly_t A, slong i, const fmpz_mpoly_ctx_t ctx);

    ulong fmpz_mpoly_get_term_var_exp_ui(const fmpz_mpoly_t A, slong i, slong var, const fmpz_mpoly_ctx_t ctx);

    slong fmpz_mpoly_get_term_var_exp_si(const fmpz_mpoly_t A, slong i, slong var, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_set_term_exp_fmpz(fmpz_mpoly_t A, slong i, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_set_term_exp_ui(fmpz_mpoly_t A, slong i, const ulong * exp, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_get_term(fmpz_mpoly_t M, const fmpz_mpoly_t A, slong i, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_get_term_monomial(fmpz_mpoly_t M, const fmpz_mpoly_t A, slong i, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_push_term_fmpz_fmpz(fmpz_mpoly_t A, const fmpz_t c, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_push_term_ui_fmpz(fmpz_mpoly_t A, ulong c, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_push_term_si_fmpz(fmpz_mpoly_t A, slong c, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_push_term_fmpz_ui(fmpz_mpoly_t A, const fmpz_t c, const ulong * exp, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_push_term_ui_ui(fmpz_mpoly_t A, ulong c, const ulong * exp, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_push_term_si_ui(fmpz_mpoly_t A, slong c, const ulong * exp, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_sort_terms(fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_combine_like_terms(fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_reverse(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_assert_canonical(const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx);

    #void _fmpz_mpoly_radix_sort1(fmpz_mpoly_t A, slong left, slong right, flint_bitcnt_t pos, ulong cmpmask, ulong totalmask);

    #void _fmpz_mpoly_radix_sort(fmpz_mpoly_t A, slong left, slong right, flint_bitcnt_t pos, slong N, ulong * cmpmask);

    void _fmpz_mpoly_push_exp_ffmpz(fmpz_mpoly_t A, const fmpz * exp, const fmpz_mpoly_ctx_t ctx);

    void _fmpz_mpoly_push_exp_pfmpz(fmpz_mpoly_t A, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx);

    void _fmpz_mpoly_push_exp_ui(fmpz_mpoly_t A, const ulong * exp, const fmpz_mpoly_ctx_t ctx);


# Random generation

    #void fmpz_mpoly_randtest_bound(fmpz_mpoly_t A, flint_rand_t state, slong length, flint_bitcnt_t coeff_bits, ulong exp_bound, const fmpz_mpoly_ctx_t ctx);

    #void fmpz_mpoly_randtest_bounds(fmpz_mpoly_t A, flint_rand_t state, slong length, flint_bitcnt_t coeff_bits, ulong * exp_bounds, const fmpz_mpoly_ctx_t ctx);

    #void fmpz_mpoly_randtest_bits(fmpz_mpoly_t A, flint_rand_t state, slong length, flint_bitcnt_t coeff_bits, flint_bitcnt_t exp_bits, const fmpz_mpoly_ctx_t ctx);


# Addition/Subtraction

    void fmpz_mpoly_add_fmpz(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_t c, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_add_ui(fmpz_mpoly_t A, const fmpz_mpoly_t B, ulong c, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_add_si(fmpz_mpoly_t A, const fmpz_mpoly_t B, slong c, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_sub_fmpz(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_t c, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_sub_ui(fmpz_mpoly_t A, const fmpz_mpoly_t B, ulong c, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_sub_si(fmpz_mpoly_t A, const fmpz_mpoly_t B, slong c, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_add(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_t C, const fmpz_mpoly_ctx_t ctx);

    slong _fmpz_mpoly_add(fmpz * poly1, ulong * exps1, const fmpz * poly2, const ulong * exps2, slong len2, const fmpz * poly3, const ulong * exps3, slong len3, slong N, const ulong * cmpmask);

    void fmpz_mpoly_sub(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_t C, const fmpz_mpoly_ctx_t ctx);

    slong _fmpz_mpoly_sub(fmpz * poly1, ulong * exps1, const fmpz * poly2, const ulong * exps2, slong len2, const fmpz * poly3, const ulong * exps3, slong len3, slong N, const ulong * cmpmask);


# Scalar operations

    void fmpz_mpoly_neg(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_scalar_mul_fmpz(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_t c, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_scalar_mul_si(fmpz_mpoly_t A, const fmpz_mpoly_t B, slong c, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_scalar_mul_ui(fmpz_mpoly_t A, const fmpz_mpoly_t B, ulong c, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_scalar_divexact_fmpz(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_t c, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_scalar_divexact_si(fmpz_mpoly_t A, const fmpz_mpoly_t B, slong c, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_scalar_divexact_ui(fmpz_mpoly_t A, const fmpz_mpoly_t B, ulong c, const fmpz_mpoly_ctx_t ctx);

    int fmpz_mpoly_scalar_divides_fmpz(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_t c, const fmpz_mpoly_ctx_t ctx);

    int fmpz_mpoly_scalar_divides_si(fmpz_mpoly_t A, const fmpz_mpoly_t B, slong c, const fmpz_mpoly_ctx_t ctx);

    int fmpz_mpoly_scalar_divides_ui(fmpz_mpoly_t A, const fmpz_mpoly_t B, ulong c, const fmpz_mpoly_ctx_t ctx);


# Differentiation/Integration

    void fmpz_mpoly_derivative(fmpz_mpoly_t A, const fmpz_mpoly_t B, slong var, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_integral(fmpz_mpoly_t A, fmpz_t scale, const fmpz_mpoly_t B, slong var, const fmpz_mpoly_ctx_t ctx);


# Evaluation

    #int _fmpz_pow_ui_is_not_feasible(flint_bitcnt_t bbits, ulong e);

    #int _fmpz_pow_fmpz_is_not_feasible(flint_bitcnt_t bbits, const fmpz_t e);

    int fmpz_mpoly_evaluate_all_fmpz(fmpz_t ev, const fmpz_mpoly_t A, fmpz * const * vals, const fmpz_mpoly_ctx_t ctx);

    int fmpz_mpoly_evaluate_one_fmpz(fmpz_mpoly_t A, const fmpz_mpoly_t B, slong var, const fmpz_t val, const fmpz_mpoly_ctx_t ctx);

    int fmpz_mpoly_compose_fmpz_poly(fmpz_poly_t A, const fmpz_mpoly_t B, fmpz_poly_struct * const * C, const fmpz_mpoly_ctx_t ctxB);

    void _fmpz_mpoly_compose_mat(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mat_t M, const fmpz_mpoly_ctx_t ctxB, const fmpz_mpoly_ctx_t ctxAC);

    int fmpz_mpoly_compose_fmpz_mpoly_geobucket(fmpz_mpoly_t A, const fmpz_mpoly_t B, fmpz_mpoly_struct * const * C, const fmpz_mpoly_ctx_t ctxB, const fmpz_mpoly_ctx_t ctxAC);

    int fmpz_mpoly_compose_fmpz_mpoly_horner(fmpz_mpoly_t A, const fmpz_mpoly_t B, fmpz_mpoly_struct * const * C, const fmpz_mpoly_ctx_t ctxB, const fmpz_mpoly_ctx_t ctxAC);

    int fmpz_mpoly_compose_fmpz_mpoly(fmpz_mpoly_t A, const fmpz_mpoly_t B, fmpz_mpoly_struct * const * C, const fmpz_mpoly_ctx_t ctxB, const fmpz_mpoly_ctx_t ctxAC);

    void fmpz_mpoly_compose_fmpz_mpoly_gen(fmpz_mpoly_t A, const fmpz_mpoly_t B, const slong * c, const fmpz_mpoly_ctx_t ctxB, const fmpz_mpoly_ctx_t ctxAC);


# Multiplication

    void fmpz_mpoly_mul(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_t C, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_mul_johnson(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_t C, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_mul_heap_threaded(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_t C, const fmpz_mpoly_ctx_t ctx);

    int fmpz_mpoly_mul_array(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_t C, const fmpz_mpoly_ctx_t ctx);

    int fmpz_mpoly_mul_array_threaded(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_t C, const fmpz_mpoly_ctx_t ctx);

    int fmpz_mpoly_mul_dense(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_t C, const fmpz_mpoly_ctx_t ctx);

    #slong _fmpz_mpoly_mul_johnson(fmpz ** poly1, ulong ** exp1, slong * alloc, const fmpz * poly2, const ulong * exp2, slong len2, const fmpz * poly3, const ulong * exp3, slong len3, flint_bitcnt_t bits, slong N, const ulong * cmpmask);

    void _fmpz_mpoly_mul_johnson_maxfields(fmpz_mpoly_t A, const fmpz_mpoly_t B, fmpz * maxBfields, const fmpz_mpoly_t C, fmpz * maxCfields, const fmpz_mpoly_ctx_t ctx);

    #void _fmpz_mpoly_mul_heap_threaded_pool_maxfields(fmpz_mpoly_t A, const fmpz_mpoly_t B, fmpz * maxBfields, const fmpz_mpoly_t C, fmpz * maxCfields, const fmpz_mpoly_ctx_t ctx, const thread_pool_handle * handles, slong num_handles);

    int _fmpz_mpoly_mul_array_DEG(fmpz_mpoly_t A, const fmpz_mpoly_t B, fmpz * maxBfields, const fmpz_mpoly_t C, fmpz * maxCfields, const fmpz_mpoly_ctx_t ctx);

    int _fmpz_mpoly_mul_array_LEX(fmpz_mpoly_t A, const fmpz_mpoly_t B, fmpz * maxBfields, const fmpz_mpoly_t C, fmpz * maxCfields, const fmpz_mpoly_ctx_t ctx);

    #int _fmpz_mpoly_mul_array_threaded_pool_DEG(fmpz_mpoly_t A, const fmpz_mpoly_t B, fmpz * maxBfields, const fmpz_mpoly_t C, fmpz * maxCfields, const fmpz_mpoly_ctx_t ctx, const thread_pool_handle * handles, slong num_handles);

    #int _fmpz_mpoly_mul_array_threaded_pool_LEX(fmpz_mpoly_t A, const fmpz_mpoly_t B, fmpz * maxBfields, const fmpz_mpoly_t C, fmpz * maxCfields, const fmpz_mpoly_ctx_t ctx, const thread_pool_handle * handles, slong num_handles);

    int _fmpz_mpoly_mul_dense(fmpz_mpoly_t P, const fmpz_mpoly_t A, fmpz * maxAfields, const fmpz_mpoly_t B, fmpz * maxBfields, const fmpz_mpoly_ctx_t ctx);

# Combined Addition and Multiplication

    void fmpz_mpoly_addmul_multi(fmpz_mpoly_t A, const fmpz_mpoly_struct ** Blist, const slong * Blengths, const slong Bnumseq, const fmpz_mpoly_ctx_t ctx);

# Powering

    int fmpz_mpoly_pow_fmpz(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_t k, const fmpz_mpoly_ctx_t ctx);

    int fmpz_mpoly_pow_ui(fmpz_mpoly_t A, const fmpz_mpoly_t B, ulong k, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_pow_fps(fmpz_mpoly_t A, const fmpz_mpoly_t B, ulong k, const fmpz_mpoly_ctx_t ctx);

    #slong _fmpz_mpoly_pow_fps(fmpz ** poly1, ulong ** exp1, slong * alloc, const fmpz * poly2, const ulong * exp2, slong len2, ulong k, flint_bitcnt_t bits, slong N, const ulong * cmpmask);


# Division

    int fmpz_mpoly_divides(fmpz_mpoly_t Q, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx);

    int fmpz_mpoly_divides_monagan_pearce(fmpz_mpoly_t Q, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx);

    int fmpz_mpoly_divides_heap_threaded(fmpz_mpoly_t Q, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx);

    #int _fmpz_mpoly_divides_heap_threaded_pool(fmpz_mpoly_t Q, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx, const thread_pool_handle * handles, slong num_handles);

    slong _fmpz_mpoly_divides_array(fmpz ** poly1, ulong ** exp1, slong * alloc, const fmpz * poly2, const ulong * exp2, slong len2, const fmpz * poly3, const ulong * exp3, slong len3, slong * mults, slong num, slong bits);

    int fmpz_mpoly_divides_array(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3, const fmpz_mpoly_ctx_t ctx);


    #int mpoly_divides_select_exps(fmpz_mpoly_t S, fmpz_mpoly_ctx_t zctx, slong nworkers, ulong * Aexp, slong Alen, ulong * Bexp, slong Blen, flint_bitcnt_t bits);

    #slong _fmpz_mpoly_divides_monagan_pearce(fmpz ** poly1, ulong ** exp1, slong * alloc, const fmpz * poly2, const ulong * exp2, slong len2, const fmpz * poly3, const ulong * exp3, slong len3, flint_bitcnt_t bits, slong N, const ulong * cmpmask);

    void fmpz_mpoly_divrem(fmpz_mpoly_t Q, fmpz_mpoly_t R, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_quasidivrem(fmpz_t scale, fmpz_mpoly_t Q, fmpz_mpoly_t R, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_div(fmpz_mpoly_t Q, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_quasidiv(fmpz_t scale, fmpz_mpoly_t Q, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_divrem_ideal(fmpz_mpoly_struct ** Q, fmpz_mpoly_t R, const fmpz_mpoly_t A, fmpz_mpoly_struct * const * B, slong len, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_quasidivrem_ideal(fmpz_t scale, fmpz_mpoly_struct ** Q, fmpz_mpoly_t R, const fmpz_mpoly_t A, fmpz_mpoly_struct * const * B, slong len, const fmpz_mpoly_ctx_t ctx);


    slong _fmpz_mpoly_div_monagan_pearce(fmpz ** polyq, ulong ** expq, slong * allocq, const fmpz * poly2, const ulong * exp2, slong len2, const fmpz * poly3, const ulong * exp3, slong len3, slong bits, slong N, const ulong * cmpmask);

    void fmpz_mpoly_div_monagan_pearce(fmpz_mpoly_t q, const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3, const fmpz_mpoly_ctx_t ctx);

    slong _fmpz_mpoly_divrem_monagan_pearce(slong * lenr, fmpz ** polyq, ulong ** expq, slong * allocq, fmpz ** polyr, ulong ** expr, slong * allocr, const fmpz * poly2, const ulong * exp2, slong len2, const fmpz * poly3, const ulong * exp3, slong len3, slong bits, slong N, const ulong * cmpmask);

    void fmpz_mpoly_divrem_monagan_pearce(fmpz_mpoly_t q, fmpz_mpoly_t r, const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3, const fmpz_mpoly_ctx_t ctx);

    slong _fmpz_mpoly_divrem_array(slong * lenr, fmpz ** polyq, ulong ** expq, slong * allocq, fmpz ** polyr, ulong ** expr, slong * allocr, const fmpz * poly2, const ulong * exp2, slong len2, const fmpz * poly3, const ulong * exp3, slong len3, slong * mults, slong num, slong bits);

    int fmpz_mpoly_divrem_array(fmpz_mpoly_t q, fmpz_mpoly_t r, const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_quasidivrem_heap(fmpz_t scale, fmpz_mpoly_t q, fmpz_mpoly_t r, const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_quasidiv_heap(fmpz_t scale, fmpz_mpoly_t q, const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3, const fmpz_mpoly_ctx_t ctx);

    slong_fmpz_mpoly_divrem_ideal_monagan_pearce(fmpz_mpoly_struct ** polyq, fmpz ** polyr, ulong ** expr, slong * allocr, const fmpz * poly2, const ulong * exp2, slong len2, fmpz_mpoly_struct * const * poly3, ulong * const * exp3, slong len, slong N, slong bits, const fmpz_mpoly_ctx_t ctx, const ulong * cmpmask);

    void fmpz_mpoly_divrem_ideal_monagan_pearce(fmpz_mpoly_struct ** q, fmpz_mpoly_t r, const fmpz_mpoly_t poly2, fmpz_mpoly_struct * const * poly3, slong len, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_quasidivrem_ideal_heap(fmpz_t scale, fmpz_mpoly_struct ** q, fmpz_mpoly_t r, const fmpz_mpoly_t poly2, fmpz_mpoly_struct * const * poly3, slong len, const fmpz_mpoly_ctx_t ctx);

# GCD

    void fmpz_mpoly_term_content(fmpz_mpoly_t M, const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx);

    int fmpz_mpoly_gcd(fmpz_mpoly_t G, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx);

    #int _fmpz_mpoly_gcd_threaded_pool(fmpz_mpoly_t G, flint_bitcnt_t Gbits, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx, const thread_pool_handle * handles, slong num_handles);

    int fmpz_mpoly_gcd_cofactors(fmpz_mpoly_t G, fmpz_mpoly_t Abar, fmpz_mpoly_t Bbar, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx);

    #int _fmpz_mpoly_gcd_cofactors_threaded_pool(fmpz_mpoly_t G, flint_bitcnt_t Gbits, fmpz_mpoly_t Abar, flint_bitcnt_t Abarbits, fmpz_mpoly_t Bbar, flint_bitcnt_t Bbarbits, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx, const thread_pool_handle * handles, slong num_handles);

    void fmpz_mpoly_deflation(fmpz * shift, fmpz * stride, const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_deflate(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz * shift, const fmpz * stride, const fmpz_mpoly_ctx_t ctx);

    void fmpz_mpoly_inflate(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz * shift, const fmpz * stride, const fmpz_mpoly_ctx_t ctx);

    int fmpz_mpoly_gcd_prs(fmpz_mpoly_t G, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx);

    int fmpz_mpoly_gcd_brown(fmpz_mpoly_t G, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx);

    int fmpz_mpoly_gcd_zippel(fmpz_mpoly_t G, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx);

    int fmpz_mpoly_gcd_berlekamp_massey(fmpz_mpoly_t G, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx);

    int fmpz_mpoly_gcd_brown_threaded(fmpz_mpoly_t G, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx);

    int fmpz_mpoly_gcd_berlekamp_massey_threaded(fmpz_mpoly_t G, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx);

# Univariates

    #void fmpz_mpoly_univar_init(fmpz_mpoly_univar_t A, const fmpz_mpoly_ctx_t ctx);

    #void fmpz_mpoly_univar_clear(fmpz_mpoly_univar_t A, const fmpz_mpoly_ctx_t ctx);

    #void fmpz_mpoly_univar_fit_length(fmpz_mpoly_univar_t A, slong length, const fmpz_mpoly_ctx_t ctx);

    #void fmpz_mpoly_univar_print_pretty(const fmpz_mpoly_univar_t A, const char ** x, const fmpz_mpoly_ctx_t ctx);

    #void fmpz_mpoly_univar_assert_canonical(fmpz_mpoly_univar_t A, const fmpz_mpoly_ctx_t ctx);

    #void fmpz_mpoly_to_univar(fmpz_mpoly_univar_t A, const fmpz_mpoly_t B, slong var, const fmpz_mpoly_ctx_t ctx);

    #void fmpz_mpoly_from_univar_bits(fmpz_mpoly_t A, flint_bitcnt_t Abits, const fmpz_mpoly_univar_t B, slong var, const fmpz_mpoly_ctx_t ctx);

    #void fmpz_mpoly_from_univar(fmpz_mpoly_t A, const fmpz_mpoly_univar_t B, slong var, const fmpz_mpoly_ctx_t ctx);

# mpolyd

    #void fmpz_mpolyd_init(fmpz_mpolyd_t poly, slong nvars);

    #void fmpz_mpolyd_fit_length(fmpz_mpolyd_t poly, slong len);

    #void fmpz_mpolyd_clear(fmpz_mpolyd_t poly);

# mpolyu

    #int fmpz_mpolyu_is_canonical(const fmpz_mpolyu_t A, const fmpz_mpoly_ctx_t ctx);

    #void fmpz_mpolyu_init(fmpz_mpolyu_t A, flint_bitcnt_t bits, const fmpz_mpoly_ctx_t ctx);

    #void fmpz_mpolyu_clear(fmpz_mpolyu_t A, const fmpz_mpoly_ctx_t uctx);

    #void fmpz_mpolyu_swap(fmpz_mpolyu_t A, fmpz_mpolyu_t B, const fmpz_mpoly_ctx_t uctx);

    #void fmpz_mpolyu_zero(fmpz_mpolyu_t A, const fmpz_mpoly_ctx_t uctx);

    #void fmpz_mpolyu_print_pretty(const fmpz_mpolyu_t poly, const char ** x, const fmpz_mpoly_ctx_t ctx);

    #void fmpz_mpolyuu_print_pretty(const fmpz_mpolyu_t poly, const char ** x, slong nmainvars, const fmpz_mpoly_ctx_t ctx);

    #void fmpz_mpolyu_fit_length(fmpz_mpolyu_t A, slong length, const fmpz_mpoly_ctx_t uctx);

    #void fmpz_mpolyu_one(fmpz_mpolyu_t A, const fmpz_mpoly_ctx_t uctx);

    #void fmpz_mpolyu_set(fmpz_mpolyu_t A, const fmpz_mpolyu_t B, const fmpz_mpoly_ctx_t uctx);

    #void fmpz_mpoly_to_mpoly_perm_deflate_threaded_pool(fmpz_mpoly_t A, const fmpz_mpoly_ctx_t lctx, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx, const slong * perm, const ulong * shift, const ulong * stride, const thread_pool_handle * handles, slong num_handles);

    #void fmpz_mpoly_from_mpoly_perm_inflate(fmpz_mpoly_t A, flint_bitcnt_t Abits, const fmpz_mpoly_ctx_t ctx, const fmpz_mpoly_t B,  const fmpz_mpoly_ctx_t lctx, const slong * perm, const ulong * shift, const ulong * stride);

    #void fmpz_mpoly_to_mpolyu_perm_deflate_threaded_pool(fmpz_mpolyu_t A, const fmpz_mpoly_ctx_t uctx, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx, const slong * perm, const ulong * shift, const ulong * stride, const ulong * maxexps, const thread_pool_handle * handles, slong num_handles);

    #void fmpz_mpoly_from_mpolyu_perm_inflate(fmpz_mpoly_t A, flint_bitcnt_t Abits, const fmpz_mpoly_ctx_t ctx, const fmpz_mpolyu_t B,  const fmpz_mpoly_ctx_t uctx, const slong * perm, const ulong * shift, const ulong * stride);

    #void fmpz_mpoly_to_mpolyuu_perm_deflate_threaded_pool(fmpz_mpolyu_t A, const fmpz_mpoly_ctx_t uctx, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx, const slong * perm, const ulong * shift, const ulong * stride, const ulong * maxexps, const thread_pool_handle * handles, slong num_handles);

    #void fmpz_mpoly_from_mpolyuu_perm_inflate(fmpz_mpoly_t A, flint_bitcnt_t Abits, const fmpz_mpoly_ctx_t ctx, const fmpz_mpolyu_t B, const fmpz_mpoly_ctx_t uctx, const slong * perm, const ulong * shift, const ulong * stride);

    #int fmpz_mpolyuu_divides(fmpz_mpolyu_t Q, const fmpz_mpolyu_t A, const fmpz_mpolyu_t B, slong main_nvars, const fmpz_mpoly_ctx_t ctx);

    #int fmpz_mpolyuu_divides_threaded_pool(fmpz_mpolyu_t Q, const fmpz_mpolyu_t A, const fmpz_mpolyu_t B, slong main_nvars, const fmpz_mpoly_ctx_t ctx, const thread_pool_handle * handles, slong num_handles);

    #void fmpz_mpolyu_fmpz_content(fmpz_t c, fmpz_mpolyu_t A, const fmpz_mpoly_ctx_t ctx);

    #void fmpz_mpolyu_mul_fmpz(fmpz_mpolyu_t A, fmpz_mpolyu_t B, fmpz_t c, const fmpz_mpoly_ctx_t ctx);

    #void fmpz_mpolyu_divexact_fmpz(fmpz_mpolyu_t A, fmpz_mpolyu_t B, fmpz_t c, const fmpz_mpoly_ctx_t ctx);

    #void fmpz_mpolyu_divexact_mpoly(fmpz_mpolyu_t A, fmpz_mpolyu_t B, int saveB, fmpz_mpoly_t c, const fmpz_mpoly_ctx_t ctx);

    #void fmpz_mpolyu_divexact_mpoly_inplace(fmpz_mpolyu_t A, fmpz_mpoly_t c, const fmpz_mpoly_ctx_t ctx);

    #void fmpz_mpolyu_mul_mpoly(fmpz_mpolyu_t A, fmpz_mpolyu_t B, fmpz_mpoly_t c, const fmpz_mpoly_ctx_t ctx);

    #void fmpz_mpolyu_mul_mpoly_inplace(fmpz_mpolyu_t A, fmpz_mpoly_t c, const fmpz_mpoly_ctx_t ctx);

    #void fmpz_mpolyu_shift_right(fmpz_mpolyu_t A, ulong s);

    #void fmpz_mpolyu_shift_left(fmpz_mpolyu_t A, ulong s);

    #void fmpz_mpolyu_content_fmpz(fmpz_t g, const fmpz_mpolyu_t A, const fmpz_mpoly_ctx_t ctx);

    #int fmpz_mpolyu_content_mpoly_threaded_pool(fmpz_mpoly_t g, const fmpz_mpolyu_t A, const fmpz_mpoly_ctx_t ctx, const thread_pool_handle * handles, slong num_handles);

    #void fmpz_mpoly_height(fmpz_t max, const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx);

    #void fmpz_mpoly_heights(fmpz_t max, fmpz_t sum, const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx);

    #int fmpz_mpolyl_gcd_brown(fmpz_mpoly_t G, fmpz_mpoly_t Abar, fmpz_mpoly_t Bbar, fmpz_mpoly_t A, fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx, const mpoly_gcd_info_t I);

    #int fmpz_mpolyl_gcd_brown_threaded_pool(fmpz_mpoly_t G, fmpz_mpoly_t Abar, fmpz_mpoly_t Bbar, fmpz_mpoly_t A, fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx, const mpoly_gcd_info_t I, const thread_pool_handle * handles, slong num_handles);

    #int fmpz_mpolyu_gcdm_zippel(fmpz_mpolyu_t G, fmpz_mpolyu_t Abar, fmpz_mpolyu_t Bbar, fmpz_mpolyu_t A, fmpz_mpolyu_t B, const fmpz_mpoly_ctx_t ctx, mpoly_zipinfo_t zinfo, flint_rand_t randstate);

    #int fmpz_mpolyuu_gcd_berlekamp_massey(fmpz_mpolyu_t G, fmpz_mpolyu_t Abar, fmpz_mpolyu_t Bbar, fmpz_mpolyu_t A, fmpz_mpolyu_t B, const fmpz_mpoly_t Gamma, const fmpz_mpoly_ctx_t ctx);

    #int fmpz_mpolyuu_gcd_berlekamp_massey_threaded_pool(fmpz_mpolyu_t G, fmpz_mpolyu_t Abar, fmpz_mpolyu_t Bbar, fmpz_mpolyu_t A, fmpz_mpolyu_t B, const fmpz_mpoly_t Gamma, const fmpz_mpoly_ctx_t ctx, const thread_pool_handle * handles, slong num_handles);

# gcd_helper_eval_interp

    #void fmpz_mpolyu_interp_reduce_p(nmod_mpolyu_t Ap, const nmod_mpoly_ctx_t ctxp, fmpz_mpolyu_t A, const fmpz_mpoly_ctx_t ctx);

    #void fmpz_mpolyu_interp_lift_p(fmpz_mpolyu_t A, const fmpz_mpoly_ctx_t ctx, nmod_mpolyu_t Ap, const nmod_mpoly_ctx_t ctxp);

    #int fmpz_mpolyu_interp_mcrt_p(flint_bitcnt_t * coeffbits, fmpz_mpolyu_t H, const fmpz_mpoly_ctx_t ctx, const fmpz_t m, const nmod_mpolyu_t A, const nmod_mpoly_ctx_t ctxp);

    #void fmpz_mpoly_interp_reduce_p_mpolyn(nmod_mpolyn_t E, const nmod_mpoly_ctx_t pctx, const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx);

    #void fmpz_mpoly_interp_lift_p_mpolyn(fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx, const nmod_mpolyn_t B, const nmod_mpoly_ctx_t pctx);

    #int fmpz_mpoly_interp_crt_p_mpolyn(fmpz_mpoly_t F, fmpz_mpoly_t T, const fmpz_mpoly_ctx_t ctx, fmpz_t modulus, const nmod_mpolyn_t A, const nmod_mpoly_ctx_t pctx);

# geobuckets

    #void fmpz_mpoly_geobucket_init(fmpz_mpoly_geobucket_t B, const fmpz_mpoly_ctx_t ctx);

    #void fmpz_mpoly_geobucket_clear(fmpz_mpoly_geobucket_t B, const fmpz_mpoly_ctx_t ctx);

    #void fmpz_mpoly_geobucket_empty(fmpz_mpoly_t p, fmpz_mpoly_geobucket_t B, const fmpz_mpoly_ctx_t ctx);

    #void fmpz_mpoly_geobucket_print(fmpz_mpoly_geobucket_t B, const char ** x, const fmpz_mpoly_ctx_t ctx);

    #void fmpz_mpoly_geobucket_fit_length(fmpz_mpoly_geobucket_t B, slong i, const fmpz_mpoly_ctx_t ctx);

    #void _fmpz_mpoly_geobucket_fix(fmpz_mpoly_geobucket_t B, slong i, const fmpz_mpoly_ctx_t ctx);

    #void fmpz_mpoly_geobucket_set(fmpz_mpoly_geobucket_t B, fmpz_mpoly_t p, const fmpz_mpoly_ctx_t ctx);

    #void fmpz_mpoly_geobucket_add(fmpz_mpoly_geobucket_t B, fmpz_mpoly_t p, const fmpz_mpoly_ctx_t ctx);

    #void fmpz_mpoly_geobucket_sub(fmpz_mpoly_geobucket_t B, fmpz_mpoly_t p, const fmpz_mpoly_ctx_t ctx);

    #void fmpz_mpoly_geobucket_set_fmpz(fmpz_mpoly_geobucket_t B, fmpz_t c, const fmpz_mpoly_ctx_t ctx);

    #void fmpz_mpoly_geobucket_gen(fmpz_mpoly_geobucket_t B, slong var, const fmpz_mpoly_ctx_t ctx);

    #void fmpz_mpoly_geobucket_add_inplace(fmpz_mpoly_geobucket_t B1, fmpz_mpoly_geobucket_t B2, const fmpz_mpoly_ctx_t ctx);

    #void fmpz_mpoly_geobucket_sub_inplace(fmpz_mpoly_geobucket_t B1, fmpz_mpoly_geobucket_t B2, const fmpz_mpoly_ctx_t ctx);

    #void fmpz_mpoly_geobucket_neg_inplace(fmpz_mpoly_geobucket_t B1, const fmpz_mpoly_ctx_t ctx);

    #void fmpz_mpoly_geobucket_mul_inplace(fmpz_mpoly_geobucket_t B1, fmpz_mpoly_geobucket_t B2, const fmpz_mpoly_ctx_t ctx);

    #void fmpz_mpoly_geobucket_pow_ui_inplace(fmpz_mpoly_geobucket_t B1, ulong k, const fmpz_mpoly_ctx_t ctx);

    #void fmpz_mpoly_geobucket_pow_fmpz_inplace(fmpz_mpoly_geobucket_t B1, const fmpz_t k, const fmpz_mpoly_ctx_t ctx);

    #int fmpz_mpoly_geobucket_divides_inplace(fmpz_mpoly_geobucket_t B1, fmpz_mpoly_geobucket_t B2, const fmpz_mpoly_ctx_t ctx);


# Helpers for gcd_berlekamp_massey

#   nmod_mpoly "skeletons" - just the coefficients

    #void nmod_mpolyc_init(nmod_mpolyc_t A);

    #void nmod_mpolyc_clear(nmod_mpolyc_t A);

    #void nmod_mpolyc_fit_length(nmod_mpolyc_t A, slong length);

    #void nmod_mpolycu_init(nmod_mpolycu_t A);

    #void nmod_mpolycu_clear(nmod_mpolycu_t A);

    #void nmod_mpolycu_fit_length(nmod_mpolycu_t A, slong length);

#   fmpz_mpoly "skeletons" - just the coefficients


    #void fmpz_mpolyc_init(fmpz_mpolyc_t A);

    #void fmpz_mpolyc_clear(fmpz_mpolyc_t A);

    #void fmpz_mpolyc_fit_length(fmpz_mpolyc_t A, slong length);

    #void fmpz_mpolycu_init(fmpz_mpolycu_t A);

    #void fmpz_mpolycu_clear(fmpz_mpolycu_t A);

    #void fmpz_mpolycu_fit_length(fmpz_mpolycu_t A, slong length);

    #void nmod_zip_init(nmod_zip_t Z);

    #void nmod_zip_clear(nmod_zip_t Z);

    #void nmod_zip_set_lengths(nmod_zip_t A, slong mlength, slong elength);

    #void nmod_zip_print(const nmod_zip_t Z, slong elength);

    #void nmod_zip_mpolyu_init(nmod_zip_mpolyu_t Z);

    #void nmod_zip_mpolyu_clear(nmod_zip_mpolyu_t Z);

    #void nmod_zip_mpolyu_fit_length(nmod_zip_mpolyu_t A, slong length);

    #void nmod_zip_mpolyu_fit_poly(nmod_zip_mpolyu_t Z, fmpz_mpolyu_t H, slong eval_length);

    #void nmod_zip_mpolyu_set_skel(nmod_zip_mpolyu_t Z, const nmod_mpoly_ctx_t ctx_sp, const fmpz_mpolyu_t A, const mp_limb_t * alpha, const fmpz_mpoly_ctx_t ctx);

    #void nmod_zip_mpolyuu_print(const nmod_zip_mpolyu_t A);

    #int nmod_zip_mpolyuu_add_point(nmod_zip_mpolyu_t L, const nmod_mpolyn_t A);

    #nmod_zip_find_coeffs_ret_t nmod_zip_find_coeffs(nmod_zip_t Z, nmod_poly_t master, slong pointcount, const nmodf_ctx_t ffinfo);

    #nmod_zip_find_coeffs_ret_t nmod_mpolyu_zip_find_coeffs(nmod_zip_mpolyu_t Z, const nmod_mpoly_ctx_t ctx_sp);

    #int fmpz_mpolyu_addinterp_zip(fmpz_mpolyu_t H, const fmpz_t Hmodulus, const nmod_zip_mpolyu_t Z, const nmodf_ctx_t ffinfo);

    #int fmpz_mpoly_repack_bits_inplace(fmpz_mpoly_t A, flint_bitcnt_t Abits, const fmpz_mpoly_ctx_t ctx);

    #int fmpz_mpolyu_repack_bits(fmpz_mpolyu_t A, flint_bitcnt_t Abits, const fmpz_mpoly_ctx_t ctx);

    #void fmpz_mpoly_set_bits(fmpz_mpoly_t A, flint_bitcnt_t Abits, const fmpz_mpoly_ctx_t ctx);

    #void fmpz_mpolyu_set_bits(fmpz_mpolyu_t A, flint_bitcnt_t Abits, const fmpz_mpoly_ctx_t ctx);

    #int fmpz_mpolyuu_eval_all_but_one_nmod(slong * out_deg, nmod_poly_t out, const fmpz_mpolyu_t A, slong var, mp_limb_t * values, const fmpz_mpoly_ctx_t ctx);

    #slong fmpz_mpolyuu_gcd_degree_bound_minor(slong * Adeg, slong * Bdeg, const fmpz_mpolyu_t A, const fmpz_mpolyu_t B, slong var, const fmpz_mpoly_ctx_t ctx, flint_rand_t state);

    #void fmpz_mpoly_ksub_content(fmpz_t content, const fmpz_mpoly_t A, const ulong * subdegs, const fmpz_mpoly_ctx_t ctx);

# Helpers for array methods

    #void _fmpz_mpoly_mul_array_chunked_DEG(fmpz_mpoly_t P, const fmpz_mpoly_t A, const fmpz_mpoly_t B, ulong degb, const fmpz_mpoly_ctx_t ctx);

    #void _fmpz_mpoly_mul_array_chunked_LEX(fmpz_mpoly_t P, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const ulong * mults, const fmpz_mpoly_ctx_t ctx);

    #void _fmpz_mpoly_addmul_array1_slong1(ulong * poly1, const slong * poly2, const ulong * exp2, slong len2, const slong * poly3, const ulong * exp3, slong len3);

    #void _fmpz_mpoly_addmul_array1_slong(ulong * poly1, const slong * poly2, const ulong * exp2, slong len2, const slong * poly3, const ulong * exp3, slong len3);

    #void _fmpz_mpoly_addmul_array1_slong2(ulong * poly1, const slong * poly2, const ulong * exp2, slong len2, const slong * poly3, const ulong * exp3, slong len3);

    #void _fmpz_mpoly_addmul_array1_fmpz(fmpz * poly1, const fmpz * poly2, const ulong * exp2, slong len2, const fmpz * poly3, const ulong * exp3, slong len3);

    #void _fmpz_mpoly_submul_array1_slong(ulong * poly1, const slong * poly2, const ulong * exp2, slong len2, const slong * poly3, const ulong * exp3, slong len3);

    #void _fmpz_mpoly_submul_array1_slong2(ulong * poly1, const slong * poly2, const ulong * exp2, slong len2, const slong * poly3, const ulong * exp3, slong len3);

    #void _fmpz_mpoly_submul_array1_slong1(ulong * poly1, const slong * poly2, const ulong * exp2, slong len2, const slong * poly3, const ulong * exp3, slong len3);

    #void _fmpz_mpoly_submul_array1_fmpz(fmpz * poly1, const fmpz * poly2, const ulong * exp2, slong len2, const fmpz * poly3, const ulong * exp3, slong len3);

    #void _fmpz_mpoly_submul_array1_slong_1(ulong * poly1, slong d, const ulong exp2, const slong * poly3, const ulong * exp3, slong len3);

    #void _fmpz_mpoly_submul_array1_slong2_1(ulong * poly1, slong d, const ulong exp2, const slong * poly3, const ulong * exp3, slong len3);

    #void _fmpz_mpoly_submul_array1_fmpz_1(fmpz * poly1, const fmpz_t d, ulong exp2, const fmpz * poly3, const ulong * exp3, slong len3);

    #void mpoly_main_variable_split_LEX(slong * ind, ulong * pexp, const ulong * Aexp, slong l1, slong Alen, const ulong * mults, slong num, slong Abits);

    #void mpoly_main_variable_split_DEG(slong * ind, ulong * pexp, const ulong * Aexp, slong l1, slong Alen, ulong deg, slong num, slong Abits);

    #slong fmpz_mpoly_append_array_sm1_LEX(fmpz_mpoly_t P, slong Plen, ulong * coeff_array, const ulong * mults, slong num, slong array_size, slong top);
    #slong fmpz_mpoly_append_array_sm2_LEX(fmpz_mpoly_t P, slong Plen, ulong * coeff_array, const ulong * mults, slong num, slong array_size, slong top);
    #slong fmpz_mpoly_append_array_sm3_LEX(fmpz_mpoly_t P, slong Plen, ulong * coeff_array, const ulong * mults, slong num, slong array_size, slong top);
    #slong fmpz_mpoly_append_array_fmpz_LEX(fmpz_mpoly_t P, slong Plen, fmpz * coeff_array, const ulong * mults, slong num, slong array_size, slong top);

    #slong fmpz_mpoly_append_array_sm1_DEGLEX(fmpz_mpoly_t P, slong Plen, ulong * coeff_array, slong top, slong nvars, slong degb);
    #slong fmpz_mpoly_append_array_sm2_DEGLEX(fmpz_mpoly_t P, slong Plen, ulong * coeff_array, slong top, slong nvars, slong degb);
    #slong fmpz_mpoly_append_array_sm3_DEGLEX(fmpz_mpoly_t P, slong Plen, ulong * coeff_array, slong top, slong nvars, slong degb);
    #slong fmpz_mpoly_append_array_fmpz_DEGLEX(fmpz_mpoly_t P, slong Plen, fmpz * coeff_array, slong top, slong nvars, slong degb);

    #slong fmpz_mpoly_append_array_sm1_DEGREVLEX(fmpz_mpoly_t P, slong Plen, ulong * coeff_array, slong top, slong nvars, slong degb);
    #slong fmpz_mpoly_append_array_sm2_DEGREVLEX(fmpz_mpoly_t P, slong Plen, ulong * coeff_array, slong top, slong nvars, slong degb);
    #slong fmpz_mpoly_append_array_sm3_DEGREVLEX(fmpz_mpoly_t P, slong Plen, ulong * coeff_array, slong top, slong nvars, slong degb);
    #slong fmpz_mpoly_append_array_fmpz_DEGREVLEX(fmpz_mpoly_t P, slong Plen, fmpz * coeff_array, slong top, slong nvars, slong degb);

    #slong _fmpz_mpoly_from_ulong_array(fmpz ** poly1, ulong ** exp1, slong * alloc, ulong * poly2, const slong * mults, slong num, slong bits, slong k);

    #slong _fmpz_mpoly_from_ulong_array2(fmpz ** poly1, ulong ** exp1, slong * alloc, ulong * poly2, const slong * mults, slong num, slong bits, slong k);

    #slong _fmpz_mpoly_from_ulong_array1(fmpz ** poly1, ulong ** exp1, slong * alloc, ulong * poly2, const slong * mults, slong num, slong bits, slong k);

    #slong _fmpz_mpoly_from_fmpz_array(fmpz ** poly1, ulong ** exp1, slong * alloc, fmpz * poly2, const slong * mults, slong num, slong bits, slong k);

    #void _fmpz_mpoly_to_ulong_array2(ulong * p, const fmpz * coeffs, const ulong * exps, slong len);

    #void _fmpz_mpoly_to_ulong_array1(ulong * p, const fmpz * coeffs, const ulong * exps, slong len);

    #void _fmpz_mpoly_to_ulong_array(ulong * p, const fmpz * coeffs, const ulong * exps, slong len);

    #void _fmpz_mpoly_to_fmpz_array(fmpz * p, const fmpz * coeffs, const ulong * exps, slong len);
