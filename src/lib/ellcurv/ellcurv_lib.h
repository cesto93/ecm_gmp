#ifndef ELLCURV_LIB_H_		/* Include guard */
#define ELLCURV_LIB_H_

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#include "matbase_lib.h"
#include "base_lib.h"
#include "ellcurv_struct.h"

#define MAX_ITER 750
#define ELL_FACT_NOT_FOUND -1
#define ell_fact_FASE(n) (n / MAX_ITER)
#define ell_fact_ITER(n) (n % MAX_ITER)

int ell_c_delta_t(const ell_c e, mpz_t ret);	//return 1 on error

void ell_p_setrandom(ell_p * p, const mpz_t n, gmp_randstate_t state);
void ell_c_setrandom(ell_c * e, const ell_p p, gmp_randstate_t state,
		     mpz_temp * temp);

//return 1 if factor is found
int ell_add_t(const ell_c e, ell_p * r, const ell_p p, const ell_p q);	// calculate p + q
int ell_dup_t(const ell_c e, ell_p * r, const ell_p p);	// calculate 2p

//return ELL_FACT_NOT_FOUND if fact not found
int ell_c_fact(mpz_t factors[], const mpz_t n, const mpz_t b1, const mpz_t b2,
	       gmp_randstate_t state);
//return 1 if factor is found
int ell_p_mul(mpz_t factor, const ell_c e, ell_p * r, const ell_p p,
	      const mpz_t k, mpz_temp * temp, ell_rep * p_temp);
int ell_p_fact_fase_2(mpz_t factor, const ell_c e, const ell_p r,
		      const mpz_t b1, const mpz_t b2, mpz_temp * temp,
		      ell_rep * p_temp, ell_rep * rep);
int create_r_difference(ell_rep * rep, mpz_t factor, const ell_c e,
			const ell_p p, unsigned long start,
			unsigned long end, mpz_temp * temp);

int get_rep_size(const mpz_t b2);
#endif //ELLCURV_LIB_H
