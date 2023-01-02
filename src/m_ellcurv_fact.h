#ifndef M_ELLCURV_FACT_H_	/* Include guard */
#define M_ELLCURV_FACT_H_

#define NO_LINUX
#define MM_ENABLE
#define NO_THREAD
#define THREAD_NUM 4
#define NDEBUG

#include <pthread.h>
#include <unistd.h>
#include <sys/syscall.h>
#include <linux/random.h>

#include "m_ellcurv_struct.h"
#include "m_ellcurv_lib.h"
#include "mm_ellcurv_lib.h"

#define _GNU_SOURCE

#define FACT_REP_SIZE 950

#ifdef NO_THREAD
#undef THREAD_NUM
#define THREAD_NUM 1
#endif

#define ELL_FACT_NOT_FOUND -1
#define ell_fact_FASE(n, max_iter) (n / max_iter)
#define ell_fact_ITER(n, max_iter) (n % max_iter)

#ifdef MM_ENABLE
	typedef mm_fact_param fact_param;
	#define ell_fact(fact, state, param, res, fase_found) mm_ell_fact(fact, state, param, res, fase_found)
	#define fact_param_init(param, n, b1, b2, max_iter) mm_ell_fact_param_init(param, n, b1, b2, max_iter)
#else
	typedef m_fact_param fact_param;
	#define ell_fact(fact, state, param, res, fase_found) m_ell_fact(fact, state, param, res, fase_found)
	#define fact_param_init(param, n, b1, b2, max_iter) m_ell_fact_param_init(param, n, b1, b2, max_iter)
#endif

typedef struct fact_tddata {
	unsigned int index;
	pthread_t *tids;
	pthread_mutex_t *mtx;
	mpz_t fact;
	unsigned long iter_done;
	int fase_found;
	gmp_randstate_t state;
	const fact_param *const_param;
} fact_tddata;

//return ELL_FACT_NOT_FOUND if fact not found
long factorize(mpz_t factors[], const mpz_t n, unsigned long b1,
	      unsigned long b2, unsigned long max_iter);

#endif //M_ELLCURV_FACT_H
