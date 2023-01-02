#ifndef M_ELLCURV_LIB_H_	/* Include guard */
#define M_ELLCURV_LIB_H_

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#include "matbase_lib.h"
#include "base_lib.h"
#include "m_ellcurv_struct.h"

#define n_temp_addh2 3
#define n_temp_duph2 3
#define n_temp_mul n_temp_duph2
#define n_temp_diff (n_temp_duph2 + 1)
#define n_temp_fase2 (n_temp_mul + 3)

#define n_p_temp_mul 2
#define n_p_temp_fase2 (n_p_temp_mul + 2)

#define n_temp_fact n_temp_fase2 + 1	//gcd
#define n_p_temp_fact n_p_temp_fase2 + 2	//p q

typedef struct m_fact_param {
	unsigned long b1;
	unsigned long b2;
	unsigned long max_iter;
	mpz_t n;
	mpz_t k;
	unsigned char *vdiff;
} m_fact_param;

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

void m_ell_fact(mpz_t fact, gmp_randstate_t state, const m_fact_param * param, unsigned int *iter, int *fase_found);
int m_ell_fact_param_init(m_fact_param * param, const mpz_t n, unsigned long b1, unsigned long b2, unsigned long max_iter);
void m_ell_fact_param_clear(m_fact_param * param);

#endif //M_ELLCURV_LIB_H
