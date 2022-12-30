#ifndef MM_ELLCURV_LIB_H_   /* Include guard */
#define MM_ELLCURV_LIB_H_

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h> 

#include "matbase_lib.h"
#include "base_lib.h"
#include "m_ellcurv_struct.h"
#include "mpn_l.h"

#define n_temp_maddh (3 + 1)
#define n_temp_mduph (3 + 1)
#define n_temp_mmul n_temp_mduph
#define n_temp_mdiff n_temp_mduph + 1
#define n_temp_mfase2 n_temp_mmul + 4

#define n_p_temp_mmul 2
#define n_p_temp_mfase2 n_p_temp_mmul + 2

#define n_temp_mfact n_temp_mfase2 + 2 //g g_r
#define n_p_temp_mfact n_p_temp_mfase2 + 2 //p q

typedef struct mm_fact_param
{
	unsigned long b1;
	unsigned long b2;
	unsigned int max_iter;
	mform_data mdata;
	mpz_t k;
	unsigned char *vdiff;
}mm_fact_param;

void mm_ell_mul_t(const mpz_t k, const mpz_t n, mpz_t e_C2, m_ellp *r, m_ellp *p);

void mm_ell_mul(const mpz_t k, mpz_t e_C2, m_ellp *r, m_ellp *p, const mform_data *mdata, m_ellp_temp *p_temp, mpz_temp *temp);				
void mm_ell_diff(m_ellp *rep, mpz_t *beta, const unsigned long d, mpz_t e_C2, const m_ellp *p, const mform_data *mdata, mpz_temp *temp);
void mm_ell_fase2(mpz_t g, unsigned long b1, unsigned long b2, const mpz_t e_C2, const m_ellp *p, m_ellp *rep, mpz_t *beta, 
					const unsigned long d, const unsigned char vdiff[], const mform_data *mdata, m_ellp_temp *p_temp, mpz_temp *temp);					

#endif //MM_ELLCURV_LIB_H
