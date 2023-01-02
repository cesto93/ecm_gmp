#include "m_ellcurv_fact.h"

long factorize(mpz_t factors[], const mpz_t n, unsigned long b1, unsigned long b2, unsigned long max_iter)
{
	const unsigned long n_size = mpz_size(n) * mp_bits_per_limb;
	int fase_found = -1;
	unsigned long iter_done = 0;
	m_fact_param param;
	unsigned long seed;
	gmp_randstate_t state;
	mpz_t fact;

	if (syscall(SYS_getrandom, &seed, sizeof(unsigned long), 0) == -1)
		error_msg("error in getrandom at m_ell_fact\n");

	mpz_init2(fact, n_size);
	gmp_randinit_default(state);
	gmp_randseed_ui(state, seed);
	m_ell_fact_param_init(&param, n, b1, b2, max_iter);

	m_ell_fact(factors[0], state, &param, &iter_done, &fase_found);

	gmp_randclear(state);
	m_ell_fact_param_clear(&param);

	if (iter_done == max_iter) {
		return ELL_FACT_NOT_FOUND;
	} else {		//FACT FOUND
		mpz_divexact(factors[1], n, factors[0]);
		return iter_done + max_iter * fase_found;
	}
}

long mfactorize(mpz_t factors[], const mpz_t n, unsigned long b1, unsigned long b2, unsigned long max_iter)
{
	const unsigned long n_size = mpz_size(n) * mp_bits_per_limb;
	int fase_found = -1;
	unsigned long iter_done = 0;
	mm_fact_param param;
	unsigned long seed;
	gmp_randstate_t state;
	mpz_t fact;

	if (syscall(SYS_getrandom, &seed, sizeof(unsigned long), 0) == -1)
		error_msg("error in getrandom at m_ell_fact\n");

	mpz_init2(fact, n_size);
	gmp_randinit_default(state);
	gmp_randseed_ui(state, seed);

	if (mm_ell_fact_param_init(&param, n, b1, b2, max_iter)) {	// FOUND IN INVERTION OF R
		mpz_set(factors[0], param.mdata.R2);
		mpz_divexact(factors[1], n, factors[0]);
		return iter_done;	//FASE 0
	}

	mm_ell_fact(factors[0], state, &param, &iter_done, &fase_found);

	gmp_randclear(state);
	mm_ell_fact_param_clear(&param);

	if (iter_done == max_iter) {
		return ELL_FACT_NOT_FOUND;
	} else {		//FACT FOUND
		mpz_divexact(factors[1], n, factors[0]);
		return iter_done + max_iter * fase_found;
	}
}