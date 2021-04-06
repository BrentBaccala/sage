
from sage.libs.flint.types cimport *

# flint/mpoly.h
cdef extern from "flint_wrap.h":

    ctypedef enum ordering_t: ORD_LEX, ORD_DEGLEX, ORD_DEGREVLEX

# context

    # Now in types.pxd

    #ctypedef struct mpoly_ctx_struct:
    #    pass

    #ctypedef mpoly_ctx_struct mpoly_ctx_t[1]

    void mpoly_ctx_init(mpoly_ctx_t ctx, slong nvars, const ordering_t ord)

    void mpoly_ctx_init_rand(mpoly_ctx_t mctx, flint_rand_t state, slong max_nvars)

    #void mpoly_monomial_randbits_fmpz(fmpz * exp, flint_rand_t state, flint_bitcnt_t exp_bits, const mpoly_ctx_t mctx)

    void mpoly_ctx_clear(mpoly_ctx_t mctx)


# heaps

# trees

    ctypedef struct mpoly_rbnode_struct:
        pass

    ctypedef mpoly_rbnode_struct mpoly_rbnode_t[1]

    ctypedef struct mpoly_rbtree_struct:
        pass

    ctypedef mpoly_rbtree_struct mpoly_rbtree_t[1]

    void mpoly_rbtree_init(mpoly_rbtree_t tree)

    void mpoly_rbnode_clear(mpoly_rbtree_t tree, mpoly_rbnode_t node, void ** dataout, slong * keysout, slong * idx)

    void mpoly_rbtree_clear(mpoly_rbtree_t tree, void ** dataout, slong * keysout)

    mpoly_rbnode_struct * mpoly_rbtree_get(int * new_node, mpoly_rbtree_struct *tree, slong rcx)

    mpoly_rbnode_struct * mpoly_rbtree_get_fmpz(int * new_node, mpoly_rbtree_struct *tree, fmpz_t rcx)

# Orderings


#  Monomials

    void mpoly_monomial_mul_fmpz(ulong * exp2, const ulong * exp3, slong N, fmpz_t c)



# generators

    void mpoly_gen_fields_ui(ulong * exp, slong var, const mpoly_ctx_t mctx)

    void mpoly_gen_fields_fmpz(fmpz * exp, slong var, const mpoly_ctx_t mctx)

    #flint_bitcnt_t mpoly_gen_bits_required(slong var, const mpoly_ctx_t mctx)

    #void mpoly_gen_offset_shift_sp(slong * offset, slong * shift, slong var, flint_bitcnt_t bits, const mpoly_ctx_t mctx)

    #void mpoly_gen_monomial_offset_shift_sp(ulong * mexp, slong * offset, slong * shift, slong var, flint_bitcnt_t bits, const mpoly_ctx_t mctx)

    #void mpoly_gen_monomial_sp(ulong * oneexp, slong var, flint_bitcnt_t bits, const mpoly_ctx_t mctx)

    #slong mpoly_gen_offset_mp(slong var,  flint_bitcnt_t bits, const mpoly_ctx_t mctx)

    #slong mpoly_gen_monomial_offset_mp(ulong * mexp, slong var, flint_bitcnt_t bits, const mpoly_ctx_t mctx)

    void fmpz_mat_mul_vec(fmpz * v, const fmpz_mat_t M, fmpz * u)

    void mpoly_compose_mat_gen(fmpz_mat_t M, const slong * c, const mpoly_ctx_t mctxB, const mpoly_ctx_t mctxAC)

    #void mpoly_compose_mat_fill_column(fmpz_mat_t M, const ulong * Cexp, flint_bitcnt_t Cbits, slong Bvar, const mpoly_ctx_t mctxB, const mpoly_ctx_t mctxAC)

# Monomial arrays

    #void mpoly_get_cmpmask(ulong * cmpmask, slong N, flint_bitcnt_t bits, const mpoly_ctx_t mctx)

    #void mpoly_get_ovfmask(ulong * ovfmask, slong N, flint_bitcnt_t bits, const mpoly_ctx_t mctx)

    #flint_bitcnt_t mpoly_exp_bits_required_ui(const ulong * user_exp, const mpoly_ctx_t mctx)
    #flint_bitcnt_t mpoly_exp_bits_required_ffmpz(const fmpz * user_exp, const mpoly_ctx_t mctx)
    #flint_bitcnt_t mpoly_exp_bits_required_pfmpz(fmpz * const * user_exp, const mpoly_ctx_t mctx)

    #void mpoly_pack_vec_ui(ulong * exp1, const ulong * exp2, flint_bitcnt_t bits, slong nfields, slong len)
    #void mpoly_pack_vec_fmpz(ulong * exp1, const fmpz * exp2, flint_bitcnt_t bits, slong nfields, slong len)

    #void mpoly_unpack_vec_ui(ulong * exp1, const ulong * exp2, flint_bitcnt_t bits, slong nfields, slong len)

    #void mpoly_unpack_vec_fmpz(fmpz * exp1, const ulong * exp2, flint_bitcnt_t bits, slong nfields, slong len)

    void mpoly_get_monomial_ui_unpacked_ffmpz(ulong * user_exps, const fmpz * poly_exps, const mpoly_ctx_t mctx)

    void mpoly_get_monomial_ffmpz_unpacked_ffmpz(fmpz * user_exps, const fmpz * poly_exps, const mpoly_ctx_t mctx)

    void mpoly_get_monomial_pfmpz_unpacked_ffmpz(fmpz ** user_exps, const fmpz * poly_exps, const mpoly_ctx_t mctx)

    void mpoly_get_monomial_ui_unpacked_ui(ulong * user_exps, const ulong * poly_exps, const mpoly_ctx_t mctx)

    #void mpoly_get_monomial_ui_sp(ulong * user_exps, const ulong * poly_exps, flint_bitcnt_t bits, const mpoly_ctx_t mctx)

    #void mpoly_get_monomial_ui_mp(ulong * user_exps, const ulong * poly_exps, flint_bitcnt_t bits, const mpoly_ctx_t mctx)

    #void mpoly_get_monomial_si_mp(slong * user_exps, const ulong * poly_exps, flint_bitcnt_t bits, const mpoly_ctx_t mctx)

    #ulong mpoly_get_monomial_var_exp_ui_sp(const ulong * poly_exps, slong var, flint_bitcnt_t bits, const mpoly_ctx_t mctx)

    #ulong mpoly_get_monomial_var_exp_ui_mp(const ulong * poly_exps, slong var, flint_bitcnt_t bits, const mpoly_ctx_t mctx)

    #slong mpoly_get_monomial_var_exp_si_mp(const ulong * poly_exps, slong var, flint_bitcnt_t bits, const mpoly_ctx_t mctx)

    #void mpoly_get_monomial_ffmpz(fmpz * exps, const ulong * poly_exps, flint_bitcnt_t bits, const mpoly_ctx_t mctx)

    #void mpoly_get_monomial_pfmpz(fmpz ** exps, const ulong * poly_exps, flint_bitcnt_t bits, const mpoly_ctx_t mctx)

    #void mpoly_set_monomial_ui(ulong * exp1, const ulong * exp2, flint_bitcnt_t bits, const mpoly_ctx_t mctx)
    #void mpoly_set_monomial_ffmpz(ulong * exp1, const fmpz * exp2, flint_bitcnt_t bits, const mpoly_ctx_t mctx)
    #void mpoly_set_monomial_pfmpz(ulong * exp1, fmpz * const * exp2, flint_bitcnt_t bits, const mpoly_ctx_t mctx)

    #int mpoly_repack_monomials(ulong * exps1, flint_bitcnt_t bits1, const ulong * exps2, flint_bitcnt_t bits2, slong len, const mpoly_ctx_t mctx)

    void mpoly_pack_monomials_tight(ulong * exp1, const ulong * exp2, slong len, const slong * mults,  slong num, slong bits)

    void mpoly_unpack_monomials_tight(ulong * e1, ulong * e2, slong len, slong * mults, slong num, slong bits)

    int mpoly_monomial_exists(slong * index, const ulong * poly_exps, const ulong * exp, slong len, slong N, const ulong * cmpmask)

    #slong mpoly_monomial_index_ui(const ulong * Aexp, flint_bitcnt_t Abits, slong Alength, const ulong * exp, const mpoly_ctx_t mctx)

    #slong mpoly_monomial_index_pfmpz(const ulong * Aexp, flint_bitcnt_t Abits, slong Alength, fmpz * const * exp, const mpoly_ctx_t mctx)

    #slong mpoly_monomial_index_monomial(const ulong * Aexp, flint_bitcnt_t Abits, slong Alength, const ulong * Mexp, flint_bitcnt_t Mbits, const mpoly_ctx_t mctx)

    #void mpoly_min_fields_ui_sp(ulong * min_fields, const ulong * poly_exps, slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx)

    #void mpoly_min_fields_fmpz(fmpz * min_fields, const ulong * poly_exps, slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx)

    #void mpoly_max_fields_ui_sp(ulong * max_fields, const ulong * poly_exps, slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx)

    #void mpoly_max_fields_fmpz(fmpz * max_fields, const ulong * poly_exps, slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx)

    #int mpoly_degrees_fit_si(const ulong * poly_exps, slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx)

    #void mpoly_degrees_si(slong * user_degs, const ulong * poly_exps, slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx)

    #void mpoly_degrees_si_threaded(slong * user_degs, const ulong * poly_exps, slong len,  flint_bitcnt_t bits, const mpoly_ctx_t mctx, const thread_pool_handle * handles, slong num_handles)

    #void mpoly_degrees_ffmpz(fmpz * user_degs, const ulong * poly_exps, slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx)

    #void mpoly_degrees_pfmpz(fmpz ** user_degs, const ulong * poly_exps, slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx)

    #slong mpoly_degree_si(const ulong * poly_exps, slong len, flint_bitcnt_t bits, slong var, const mpoly_ctx_t mctx)

    #void mpoly_degree_fmpz(fmpz_t deg, const ulong * poly_exps, slong len, flint_bitcnt_t bits, slong var, const mpoly_ctx_t mctx)

    #int  mpoly_total_degree_fits_si(const ulong * exps, slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx)

    slong mpoly_total_degree_si(const ulong * exps, slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx)

    void mpoly_total_degree_fmpz(fmpz_t totdeg, const ulong * exps, slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx)

    #void mpoly_total_degree_fmpz_ref(fmpz_t totdeg, const ulong * exps, slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx)

    void mpoly_search_monomials( slong ** e_ind, ulong * e, slong * e_score, slong * t1, slong * t2, slong *t3, slong lower, slong upper, const ulong * a, slong a_len, const ulong * b, slong b_len, slong N, const ulong * cmpmask)

    void mpoly_main_variable_split_LEX(slong * ind, ulong * pexp, const ulong * Aexp, slong l1, slong Alen, const ulong * mults, slong num, slong Abits)

    void mpoly_main_variable_split_DEG(slong * ind, ulong * pexp, const ulong * Aexp,  slong l1, slong Alen, ulong deg, slong num, slong Abits)

    #int mpoly_term_exp_fits_si(ulong * exps, flint_bitcnt_t bits, slong n, const mpoly_ctx_t mctx)

    #int mpoly_term_exp_fits_ui(ulong * exps, flint_bitcnt_t bits, slong n, const mpoly_ctx_t mctx)

    #int mpoly_is_gen(ulong * exps, slong var, flint_bitcnt_t bits, const mpoly_ctx_t mctx)

    #int mpoly_monomials_valid_test(ulong * exps, slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx)

    #int mpoly_monomials_overflow_test(ulong * exps, slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx)

    #int mpoly_monomials_inorder_test(ulong * exps, slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx)

    void mpoly_reverse(ulong * Aexp, const ulong * Bexp, slong len, slong N)

    #void mpoly_monomials_deflation(fmpz * shift, fmpz * stride, const ulong * Aexps, flint_bitcnt_t Abits, slong Alength, const mpoly_ctx_t mctx)

    #void mpoly_monomials_deflate(ulong * Aexps, flint_bitcnt_t Abits, const ulong * Bexps, flint_bitcnt_t Bbits, slong Blength, const fmpz * shift, const fmpz * stride, const mpoly_ctx_t mctx)

    #void mpoly_monomials_inflate(ulong * Aexps, flint_bitcnt_t Abits, const ulong * Bexps, flint_bitcnt_t Bbits, slong Blength, const fmpz * shift, const fmpz * stride, const mpoly_ctx_t mctx)

    #void _mpoly_gen_shift_right(ulong * Aexp, flint_bitcnt_t Abits, slong Alength, slong var, ulong amount, const mpoly_ctx_t mctx)

    #void _mpoly_gen_shift_left(ulong * Aexp, flint_bitcnt_t Abits, slong Alength, slong var, ulong amount, const mpoly_ctx_t mctx)

    #int mpoly_monomial_cmp_general(ulong * Aexp, flint_bitcnt_t Abits, ulong * Bexp, flint_bitcnt_t Bbits, const mpoly_ctx_t mctx)

    #void mpoly_monomials_shift_right_ui(ulong * Aexps, flint_bitcnt_t Abits, slong Alength, const ulong * user_exps, const mpoly_ctx_t mctx)

    #int mpoly_monomial_cofactors(fmpz * Abarexps, fmpz * Bbarexps, const ulong * Aexps, flint_bitcnt_t Abits, const ulong * Bexps, flint_bitcnt_t Bbits, slong length,  const mpoly_ctx_t mctx)

# info related to gcd calculation

    ctypedef struct mpoly_gcd_info_struct:
        pass

    ctypedef mpoly_gcd_info_struct mpoly_gcd_info_t[1]

    void mpoly_gcd_info_init(mpoly_gcd_info_t I, slong nvars)

    void mpoly_gcd_info_clear(mpoly_gcd_info_t I)

    #void mpoly_gcd_info_limits(ulong * Amax_exp, ulong * Amin_exp, slong * Amax_exp_count, slong * Amin_exp_count, const ulong * Aexps, flint_bitcnt_t Abits, slong Alength, const mpoly_ctx_t mctx)

    #void mpoly_gcd_info_stride(ulong * strides, const ulong * Aexps, flint_bitcnt_t Abits, slong Alength, const ulong * Amax_exp, const ulong * Amin_exp, const ulong * Bexps, flint_bitcnt_t Bbits, slong Blength, const ulong * Bmax_exp, const ulong * Bmin_exp, const mpoly_ctx_t mctx)

    void mpoly_gcd_info_set_perm(mpoly_gcd_info_t I, slong Alength, slong Blength, const mpoly_ctx_t mctx)

    slong mpoly_gcd_info_get_brown_upper_limit(const mpoly_gcd_info_t I, slong var, slong bound)

    void mpoly_gcd_info_measure_brown(mpoly_gcd_info_t I, slong Alength, slong Blength, const mpoly_ctx_t mctx)

    void mpoly_gcd_info_measure_bma(mpoly_gcd_info_t I, slong Alength, slong Blength, const mpoly_ctx_t mctx)

    void mpoly_gcd_info_measure_zippel(mpoly_gcd_info_t I, slong Alength, slong Blength, const mpoly_ctx_t mctx)

    ctypedef struct mpoly_zipinfo_struct:
        pass

    ctypedef mpoly_zipinfo_struct mpoly_zipinfo_t[1]

    void mpoly_zipinfo_init(mpoly_zipinfo_t zinfo, slong nvars)

    void mpoly_zipinfo_clear(mpoly_zipinfo_t zinfo)

    void _fmpz_vec_content_chained(fmpz_t res, const fmpz * vec, slong len)
