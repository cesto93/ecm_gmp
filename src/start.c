#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#include "base_lib.h"
#include "m_ellcurv_fact.h"

#define DEFAULT_COEFF 100 // b2 = b1 * coeff
#define DEFAULT_MAX_ITER 2000

int main(int argc, char *argv[])
{
	unsigned long b1, b2, max_iter;
	struct timespec start, end, diff;
	mpz_t n;
	mpz_t factors[2];
	
	if (argc < 3)
	{
		perror("usage: ./start.o n b1 [b2] [max_iter]");
		exit(EXIT_FAILURE);
	}
	
	mpz_inits(factors[0], factors[1], n, NULL);
	
	mpz_set_str(n, argv[1], 10);
	b1 = str2l(argv[2]); //SMOOTHNESS bound
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
	int res = factorize(factors, n, b1, b2, max_iter); 
	get_current_time(&end);
	timespec_diff(&start, &end, &diff);
	
	if (res != ELL_FACT_NOT_FOUND)
		gmp_printf("factors are: \n%Zd\n%Zd\n", factors[0], factors[1]);
	else
		gmp_printf("factor not found\n");
	printf("iter %lu\t fase %lu\t time[s] %lu.%lu\n", ell_fact_ITER(res, max_iter), ell_fact_FASE(res, max_iter), diff.tv_sec, diff.tv_nsec);	
	exit(EXIT_SUCCESS);
}


