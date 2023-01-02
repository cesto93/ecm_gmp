#include "m_ellcurv_fact.h"

static void lock(pthread_mutex_t * mtx)
{
	if (pthread_mutex_lock(mtx))
		error_msg("errore nella lock");
}

static void unlock(pthread_mutex_t * mtx)
{
	if (pthread_mutex_unlock(mtx))
		error_msg("errore nella unlock");
}

static void td_data_init(fact_tddata * td, const unsigned long n_size)
{
	gmp_randinit_default(td->state);
	mpz_init2(td->fact, n_size);
	mpz_init2(td->e_C2, n_size);
	m_ellp_rep_init2(&(td->rep), FACT_REP_SIZE, n_size);
	mpz_rep_init2(&(td->beta), FACT_REP_SIZE, n_size);
}

static void td_data_clear(fact_tddata * td)
{
	gmp_randclear(td->state);
	mpz_clears(td->fact, td->e_C2, NULL);
	m_ellp_rep_clear(&(td->rep));
	mpz_rep_clear(&(td->beta));
	mpz_temp_clear(&(td->temp));
	m_ellp_temp_clear(&(td->p_temp));
}

static int fact_param_init(fact_param * param, const mpz_t n, unsigned long b1, unsigned long b2, unsigned long max_iter, mpz_temp * temp)
{
	int vdiff_size = get_vdiff_size(b2);
	if (vdiff_size == -1)
		error_msg("b2 to big\n");
	if ((param->vdiff = malloc(vdiff_size)) == NULL)
		error_msg("error in malloc at m_ell_fact\n");
	mpz_init2(param->k, bigk_size_bits(b1));

	create_bigk(param->k, b1, temp);
	get_prime_diff(b1, 1, b2, param->vdiff, temp);
	mpz_realloc2(param->k, mpz_size(param->k) * mp_bits_per_limb);
	param->b1 = b1;
	param->b2 = b2;
	param->max_iter = max_iter / THREAD_NUM;

#ifdef MM_ENABLE
	if (mform_data_init(&(param->mdata), n, temp))	//CAN'T INVERT R TO MOD N
		return 1;
#else
	mpz_init(param->n);
	mpz_set(param->n, n);
#endif

	return 0;
}

static void fact_param_clear(fact_param * param)
{
	mpz_clears(param->k, NULL);
	free(param->vdiff);

#ifdef MM_ENABLE
	mform_data_clear(&(param->mdata));
#else
	mpz_clear(param->n);
#endif
}

void *fact_tjob(void *arg)
{
	fact_tddata *const td = (fact_tddata *) arg;

	ell_fact(td->fact, td->state, td->e_C2, td->const_param, &(td->rep),
		 &(td->beta), &(td->temp), &(td->p_temp), &(td->iter_done),
		 &(td->fase_found));
	lock(td->mtx);
	if (td->fase_found != -1) {
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, NULL);
		for (unsigned int i = 0; i < THREAD_NUM; i++)
			if (i != td->index)
				pthread_cancel(td->tids[i]);
	}
	unlock(td->mtx);
	pthread_exit(0);
}

//HIGH_LEVEL FUNCT
long factorize(mpz_t factors[], const mpz_t n, unsigned long b1, unsigned long b2, unsigned long max_iter)
{
	const unsigned long n_size = mpz_size(n) * mp_bits_per_limb;
	int fase_found = -1;
	unsigned long iter_done = 0;
	fact_param param;
	fact_tddata td[THREAD_NUM];
	pthread_t tids[THREAD_NUM];
	unsigned long seeds[THREAD_NUM];
	pthread_mutex_t td_mutex;

	if (pthread_mutex_init(&td_mutex, NULL))
		error_msg("error in mutex init\n");

#ifndef NO_LINUX
	if (syscall(SYS_getrandom, seeds, sizeof(unsigned long) * THREAD_NUM, 0)
	    == -1)
		error_msg("error in getrandom at m_ell_fact\n");
#else
	seeds[0] = time(NULL);
	for (int i = 1; i < THREAD_NUM; i++)
		seeds[i] = seeds[0] + i;
#endif

	for (int i = 0; i < THREAD_NUM; i++) {
		td_data_init(td + i, n_size);
		gmp_randseed_ui(td[i].state, seeds[i]);
		td[i].tids = tids;
		td[i].index = i;
		td[i].fase_found = -1;
		td[i].mtx = &td_mutex;
	}

	//PRECALCULATE PARAM
#ifdef MM_ENABLE
	if (fact_param_init(&param, n, b1, b2, max_iter, &(td[0].temp)))	//FOUND IN INVERTION OF R
	{
		mpz_set(factors[0], param.mdata.R2);
		mpz_divexact(factors[1], n, factors[0]);
		return iter_done;	//FASE 0
	}
#else
	fact_param_init(&param, n, b1, b2, max_iter, &(td[0].temp));
#endif

#ifndef NO_THREAD
	for (int i = 0; i < THREAD_NUM; i++) {
		td[i].const_param = &param;
		if (pthread_create(&tids[i], NULL, fact_tjob, (void *)&td[i]))	//THE JOB
			error_msg("error in pthread_create");
	}

	for (int i = 0; i < THREAD_NUM; i++) {
		if (pthread_join(tids[i], NULL))
			error_msg("error in pthread_join");
		iter_done += td[i].iter_done;
		if (td[i].fase_found != -1) {
			mpz_set(factors[0], td[i].fact);
			fase_found = td[i].fase_found;
		}
	}
#else
	ell_fact(factors[0], td[0].state, td[0].e_C2, &param, &(td[0].rep),
		 &(td[0].beta), &(td[0].temp), &(td[0].p_temp), &iter_done,
		 &fase_found);
#endif

	for (int i = 0; i < THREAD_NUM; i++)
		td_data_clear(&td[i]);
	fact_param_clear(&param);

	if (iter_done == max_iter)
		return ELL_FACT_NOT_FOUND;
	else			//FACT_FOUND
	{
		mpz_divexact(factors[1], n, factors[0]);
		return iter_done + max_iter * fase_found;
	}
}
