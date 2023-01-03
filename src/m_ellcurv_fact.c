#include "m_ellcurv_fact.h"

#include <linux/random.h>
#include <sys/syscall.h>
#include <unistd.h>

#include "m_ellcurv_lib.h"
#include "mm_ellcurv_lib.h"

m_ellfact_res *factorize(const mpz_t n, unsigned long b1, unsigned long b2, unsigned long max_iter)
{
	unsigned long seed;
	gmp_randstate_t state;
	m_ellfact_res *res;

	if (syscall(SYS_getrandom, &seed, sizeof(unsigned long), 0) == -1)
		error_msg("error in getrandom at m_ell_fact\n");

	gmp_randinit_default(state);
	gmp_randseed_ui(state, seed);
	res = m_ell_fact(state, n, b1, b2, max_iter);
	gmp_randclear(state);

	return res;
}

long mfactorize(mpz_t factors[], const mpz_t n, unsigned long b1, unsigned long b2, unsigned long max_iter)
{
	int fase_found = -1;
	unsigned long iter_done = 0;
	unsigned long seed;
	gmp_randstate_t state;
	mpz_t *fact;

	if (syscall(SYS_getrandom, &seed, sizeof(unsigned long), 0) == -1)
		error_msg("error in getrandom at m_ell_fact\n");

	gmp_randinit_default(state);
	gmp_randseed_ui(state, seed);

	fact = mm_ell_fact(state, n, b1, b2, max_iter, &iter_done, &fase_found);

	gmp_randclear(state);

	if (iter_done == max_iter) {
		return ELL_FACT_NOT_FOUND;
	} else {		//FACT FOUND
		mpz_set(factors[0], *fact);
		mpz_divexact(factors[1], n, factors[0]);
		return iter_done + max_iter * fase_found;
	}
}
