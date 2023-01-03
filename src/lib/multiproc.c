#include <pthread.h>
#include <unistd.h>
#include <sys/syscall.h>
#include <linux/random.h>

#include "../bin/ellcurv_fact.h"

#define _GNU_SOURCE

typedef struct fact_tddata {
	unsigned int index;
	pthread_t *tids;
	pthread_mutex_t *mtx;
	m_ellfact_res *res;
	mpz_t n;
	unsigned long b1;
	unsigned long b2;
	unsigned long max_iter;
	unsigned long thread_num;
} fact_tddata;

static void td_data_init(fact_tddata *td, pthread_mutex_t *mtx, unsigned int index, const mpz_t n, unsigned long b1, unsigned long b2, unsigned long max_iter)
{	
	td->index = index;
	td->mtx = mtx;
	td->b1 = b1;
	td->b2 = b2;
	td->max_iter = max_iter;
	mpz_set(td->n, n);

	td->tids = malloc(sizeof(pthread_t));
	if (td->tids == NULL)
		error_msg("error in td_data_init()");
}

static void td_data_clear(fact_tddata * td)
{
	mpz_clears(td->res->fact[0], td->res->fact[1], NULL);
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
		for (unsigned int i = 0; i < td->thread_num; i++)
			if (i != td->index)
				pthread_cancel(td->tids[i]);
	}

	if (pthread_mutex_unlock(td->mtx))
		error_msg("errore nella unlock");
	
	pthread_exit(0);
}

m_ellfact_res *tfactorize(const mpz_t n, unsigned long b1, unsigned long b2, unsigned long max_iter, unsigned long thread_num)
{
	fact_tddata td[thread_num];
	pthread_mutex_t td_mutex;
	m_ellfact_res *res;
	
	res = malloc(sizeof(*res));
	if (res == NULL)
		return NULL;

	if (pthread_mutex_init(&td_mutex, NULL))
		error_msg("error in mutex init\n");

	for (unsigned long i = 0; i < thread_num; i++) {
		td_data_init(td + i, &td_mutex, i, n, b1, b2, max_iter / thread_num);
		td->thread_num = thread_num;
		if (pthread_create(td[i].tids, NULL, fact_tjob, (void *)&td[i]))	//THE JOB
			error_msg("error in pthread_create");
	}

	for (unsigned long i = 0; i < thread_num; i++) {
		if (pthread_join(*(td[i].tids), NULL))
			error_msg("error in pthread_join");
		res->iter += td[i].res->iter;
		
		if (td[i].res->fase_found != ELL_FACT_NOT_FOUND) {
			mpz_set(res->fact[0], td[i].res->fact[0]);
			mpz_set(res->fact[1], td[i].res->fact[1]);
			res->fase_found = td[i].res->fase_found;
		}
	}

	for (unsigned long i = 0; i < thread_num; i++)
		td_data_clear(&td[i]);

	return res;
}