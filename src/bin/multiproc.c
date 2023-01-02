#include <pthread.h>
#include <unistd.h>
#include <sys/syscall.h>
#include <linux/random.h>

#include "m_ellcurv_struct.h"
#include "m_ellcurv_lib.h"
#include "mm_ellcurv_lib.h"
#include "m_ellcurv_fact.h"

#define _GNU_SOURCE

#define THREAD_NUM 4

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

static void td_data_init(fact_tddata *td, const unsigned long n_size, const unsigned long seed, pthread_mutex_t *mtx, unsigned int index)
{
	gmp_randinit_default(td->state);
	mpz_init2(td->fact, n_size);

	gmp_randseed_ui(td->state, seed);
	
	td->index = index;
	td->fase_found = -1;
	td->mtx = mtx;

	td->tids = malloc(sizeof(pthread_t));
	if (td->tids == NULL)
		error_msg("error in td_data_init()");
}

static void td_data_clear(fact_tddata * td)
{
	gmp_randclear(td->state);
}

void *fact_tjob(void *arg)
{
	fact_tddata *const td = (fact_tddata *) arg;

	ell_fact(td->fact, td->state, td->const_param, &(td->iter_done), &(td->fase_found));

	if (pthread_mutex_lock(td->mtx))
		error_msg("errore nella lock");

	if (td->fase_found != -1) {
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, NULL);
		for (unsigned int i = 0; i < THREAD_NUM; i++)
			if (i != td->index)
				pthread_cancel(td->tids[i]);
	}

	if (pthread_mutex_unlock(td->mtx))
		error_msg("errore nella unlock");
	
	pthread_exit(0);
}

long factorize(mpz_t factors[], const mpz_t n, unsigned long b1, unsigned long b2, unsigned long max_iter)
{
	const unsigned long n_size = mpz_size(n) * mp_bits_per_limb;
	int fase_found = -1;
	unsigned long iter_done = 0;
	fact_param param;
	fact_tddata td[THREAD_NUM];
	unsigned long seeds[THREAD_NUM];
	pthread_mutex_t td_mutex;

	if (pthread_mutex_init(&td_mutex, NULL))
		error_msg("error in mutex init\n");

	if (syscall(SYS_getrandom, seeds, sizeof(unsigned long) * THREAD_NUM, 0) == -1)
		error_msg("error in getrandom at m_ell_fact\n");

	for (int i = 0; i < THREAD_NUM; i++) {
		td_data_init(td + i, n_size, seeds[i], &td_mutex, i);
	}

	//PRECALCULATE PARAM
#ifdef MM_ENABLE
	if (fact_param_init(&param, n, b1, b2, max_iter / THREAD_NUM)) { //FOUND IN INVERTION OF R
		mpz_set(factors[0], param.mdata.R2);
		mpz_divexact(factors[1], n, factors[0]);
		return iter_done;	//FASE 0
	}
#else
	fact_param_init(&param, n, b1, b2, max_iter / THREAD_NUM);
#endif

	for (int i = 0; i < THREAD_NUM; i++) {
		td[i].const_param = &param;
		if (pthread_create(td[i].tids, NULL, fact_tjob, (void *)&td[i]))	//THE JOB
			error_msg("error in pthread_create");
	}

	for (int i = 0; i < THREAD_NUM; i++) {
		if (pthread_join(*(td[i].tids), NULL))
			error_msg("error in pthread_join");
		iter_done += td[i].iter_done;
		if (td[i].fase_found != -1) {
			mpz_set(factors[0], td[i].fact);
			fase_found = td[i].fase_found;
		}
	}

	for (int i = 0; i < THREAD_NUM; i++)
		td_data_clear(&td[i]);
	fact_param_clear(&param);

	if (iter_done == max_iter) {
		return ELL_FACT_NOT_FOUND;
	} else { //FACT FOUND
		mpz_divexact(factors[1], n, factors[0]);
		return iter_done + max_iter * fase_found;
	}
}