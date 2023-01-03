#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#include "ellcurv_fact.h"

#define DEFAULT_COEFF 100	// b2 = b1 * coeff
#define DEFAULT_MAX_ITER 5000

int main(int argc, char *argv[])
{
	unsigned long b1, b2, max_iter;
	struct timespec start, end, diff;
	m_ellfact_res *res;
	mpz_t n;

	if (argc < 3) {
		perror("usage: ./start.o n b1 [b2] [max_iter]");
		exit(EXIT_FAILURE);
	}

	mpz_inits(n, NULL);

	mpz_set_str(n, argv[1], 10);
	b1 = str2l(argv[2]);	//SMOOTHNESS bound
	if (argc > 3)
		b2 = str2l(argv[3]);
	else
		b2 = b1 * DEFAULT_COEFF;
	if (argc > 4)
		max_iter = str2l(argv[4]);
	else
		max_iter = DEFAULT_MAX_ITER;

	gmp_printf("n = %Zd b1 = %lu b2 = %lu\n", n, b1, b2);

	get_current_time(&start);
	res = factorize(n, b1, b2, max_iter);
	get_current_time(&end);
	timespec_diff(&start, &end, &diff);

	if (res->fase_found != ELL_FACT_NOT_FOUND)
		gmp_printf("factors are: \n%Zd\n%Zd\n", res->fact[0], res->fact[1]);
	else
		gmp_printf("factor not found\n");
	printf("iter %lu\t fase %u\t time[s] %lu.%lu\n", res->iter, res->fase_found, diff.tv_sec, diff.tv_nsec);
	exit(EXIT_SUCCESS);
}
