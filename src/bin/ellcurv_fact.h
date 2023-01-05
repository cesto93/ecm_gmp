#ifndef ELLCURV_FACT_H_	/* Include guard */
#define ELLCURV_FACT_H_

#define NDEBUG

#include "../lib/m_ellcurv.h"

#include <linux/random.h>
#include <sys/syscall.h>
#include <unistd.h>

#define _GNU_SOURCE
#define FACT_REP_SIZE 950

#ifndef MM_ENABLE
#define factorize(n, b1, b2, max_iter) mfactorize(n, b1, b2, max_iter)
#else
#define factorize(n, b1, b2, max_iter) mmfactorize(n, b1, b2, max_iter)
#endif

static inline m_ellfact_res *mfactorize(const mpz_t n, unsigned long b1, unsigned long b2, unsigned long max_iter)
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

static inline m_ellfact_res *mmfactorize(const mpz_t n, unsigned long b1, unsigned long b2, unsigned long max_iter)
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

m_ellfact_res *tfactorize(const mpz_t n, unsigned long b1, unsigned long b2, unsigned long max_iter, unsigned long thread_num);

#endif //ELLCURV_FACT_H
