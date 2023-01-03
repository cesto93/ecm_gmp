#include <pthread.h>
#include <unistd.h>
#include <sys/syscall.h>
#include <linux/random.h>

#include "ellcurv_fact.h"

#define _GNU_SOURCE

#define THREAD_NUM 4

typedef struct fact_tddata {
	unsigned int index;
	pthread_t *tids;
	pthread_mutex_t *mtx;
	m_ellfact_res *res;
	const mpz_t n;
	unsigned long b1;
	unsigned long b2;
	unsigned long max_iter;
} fact_tddata;

static void td_data_init(fact_tddata *td, pthread_mutex_t *mtx, unsigned int index, const mpz_t n, unsigned long b1, unsigned long b2, unsigned long max_iter)
{	
	td->index = index;
	td->mtx = mtx;
	td->b1 = b1;
	td->b2 = b2;
	td->max_iter = max_iter;

	td->tids = malloc(sizeof(pthread_t));
	if (td->tids == NULL)
		error_msg("error in td_data_init()");
}

static void td_data_clear(fact_tddata * td)
{
	mpz_clears(res->fact[0], res->fact[1], NULL);
	free(td->res);
	free(td->tids);
}

static void *fact_tjob(void *arg)
{
	fact_tddata *const td = (fact_tddata *) arg;

	td->res = factorize(td->n, td->b1, td->b2, td->max_iter);

	if (pthread_mutex_lock(td->mtx))
		error_msg("errore nella lock");

	if (td->res->fase_found != ELL_FACT_NOT_FOUND) {
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, NULL);
		for (unsigned int i = 0; i < THREAD_NUM; i++)
			if (i != td->index)
				pthread_cancel(td->tids[i]);
	}

	if (pthread_mutex_unlock(td->mtx))
		error_msg("errore nella unlock");
	
	pthread_exit(0);
}

long tfactorize(const mpz_t n, unsigned long b1, unsigned long b2, unsigned long max_iter)
{
	fact_tddata td[THREAD_NUM];
	pthread_mutex_t td_mutex;
	m_ellfact_res *res;

	if (pthread_mutex_init(&td_mutex, NULL))
		error_msg("error in mutex init\n");

	for (int i = 0; i < THREAD_NUM; i++) {
		td_data_init(td + i, &td_mutex, i, n, b1, b2, max_iter / THREAD_NUM);
		if (pthread_create(td[i].tids, NULL, fact_tjob, (void *)&td[i]))	//THE JOB
			error_msg("error in pthread_create");
	}

	for (int i = 0; i < THREAD_NUM; i++) {
		if (pthread_join(*(td[i].tids), NULL))
			error_msg("error in pthread_join");
		res->iter += td[i].res->iter;
		
		if (td[i].res->fase_found != ELL_FACT_NOT_FOUND) {
			mpz_set(res->fact[0], td[i].res->fact[0]);
			mpz_set(res->fact[1], td[i].res->fact[1]);
			res->fase_found = td[i].res->fase_found;
		}
	}

	for (int i = 0; i < THREAD_NUM; i++)
		td_data_clear(&td[i]);

	return res;
}