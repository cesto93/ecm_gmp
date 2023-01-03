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

m_ellfact_res *mfactorize(const mpz_t n, unsigned long b1, unsigned long b2, unsigned long max_iter)
{
	unsigned long seed;
	gmp_randstate_t state;
	m_ellfact_res *res;

	if (syscall(SYS_getrandom, &seed, sizeof(unsigned long), 0) == -1)
		error_msg("error in getrandom at m_ell_fact\n");

	gmp_randinit_default(state);
	gmp_randseed_ui(state, seed);
	res = mm_ell_fact(state, n, b1, b2, max_iter);
	gmp_randclear(state);

	return res;
}
