#ifndef M_ELLCURV_H_	/* Include guard */
#define M_ELLCURV_H_

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#include "m_ellcurv_struct.h"
#include "mpn_l.h"

#define n_temp_addh2 3
#define n_temp_duph2 3
#define n_temp_mul n_temp_duph2
#define n_temp_diff (n_temp_duph2 + 1)
#define n_temp_fase2 (n_temp_mul + 3)

#define n_p_temp_mul 2
#define n_p_temp_fase2 (n_p_temp_mul + 2)

#define n_temp_fact n_temp_fase2 + 1	//gcd
#define n_p_temp_fact n_p_temp_fase2 + 2	//p q

void m_ell_addh_t(const m_ellc * e, m_ellp * r, const m_ellp * p, const m_ellp * q, const m_ellp * diff);	// calculate p + q
void m_ell_duph_t(const m_ellc * e, m_ellp * r, const m_ellp * p);	// calculate 2p
void m_ell_mul_t(const mpz_t k, const m_ellc * e, m_ellp * r, const m_ellp * p);
void m_ell_diff_t(m_ellp_rep * rep, mpz_rep * beta, const m_ellc * e, const m_ellp * p);

void m_ell_mul(const mpz_t k, const mpz_t n, const mpz_t e_C2, m_ellp * r, const m_ellp * p, m_ellp_temp * p_temp, mpz_temp * temp);
void m_ell_diff(m_ellp_rep * rep, mpz_rep * beta, const mpz_t n, const mpz_t e_C2, const m_ellp * p, mpz_temp * temp);
void m_ell_fase2(mpz_t g, unsigned long b1, unsigned long b2, const mpz_t n,
		 const mpz_t e_C2, const m_ellp * p, const m_ellp_rep rep,
		 const mpz_rep beta, const unsigned char vdiff[], m_ellp_temp * p_temp, mpz_temp * temp);

void check_mul_ui(unsigned long k, const m_ellc * e, m_ellp * r, const m_ellp * p, m_ellp_temp * p_temp, mpz_temp * temp);

m_ellfact_res *m_ell_fact(gmp_randstate_t state, const mpz_t n, unsigned long b1, unsigned long b2, unsigned long max_iter);

#define n_temp_maddh (3 + 1)
#define n_temp_mduph (3 + 1)
#define n_temp_mmul n_temp_mduph
#define n_temp_mdiff n_temp_mduph + 1
#define n_temp_mfase2 n_temp_mmul + 4

#define n_p_temp_mmul 2
#define n_p_temp_mfase2 n_p_temp_mmul + 2

#define n_temp_mfact n_temp_mfase2 + 2	//g g_r
#define n_p_temp_mfact n_p_temp_mfase2 + 2	//p q

void mm_ell_mul_t(const mpz_t k, const mpz_t n, mpz_t e_C2, m_ellp * r, m_ellp * p);

void mm_ell_mul(const mpz_t k, mpz_t e_C2, m_ellp * r, m_ellp * p, const mform_data * mdata, m_ellp_temp * p_temp, mpz_temp * temp);

void mm_ell_diff(m_ellp * rep, mpz_t * beta, const unsigned long d, const mpz_t e_C2, const m_ellp * p, const mform_data * mdata, mpz_temp * temp);

void mm_ell_fase2(mpz_t g, unsigned long b1, unsigned long b2, const mpz_t e_C2, const m_ellp * p, m_ellp * rep, mpz_t * beta,
		  const unsigned long d, const unsigned char vdiff[], const mform_data * mdata, m_ellp_temp * p_temp, mpz_temp * temp);

m_ellfact_res *mm_ell_fact(gmp_randstate_t state, const mpz_t n, unsigned long b1, unsigned long b2, unsigned long max_iter);

#endif //M_ELLCURV_H
