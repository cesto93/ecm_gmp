#ifndef MATBASE_LIB_H_		/* Include guard */
#define MATBASE_LIB_H_

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <errno.h>
#include <gmp.h>

#include "base_lib.h"

#define mpz_is_zero(n) (mpz_sgn(n) == 0)
#define mpz_equal(a, b) (mpz_cmp(a, b) == 0)

#define n_temp_get_prime_diff 3
#define n_temp_create_bigk 2
#define n_temp_mop2 1

#define bigk_size_bits(b) (b * 1.44)	//1.442695040888963407359924681001892137426645954152985934135

typedef struct mpz_temp		//temp variable container
{
	mpz_t *t;
	unsigned int lenght;
	unsigned int index;
} mpz_temp;

typedef struct mpz_rep		//rep variable container
{
	mpz_t *v;
	unsigned long lenght;
} mpz_rep;

void mpz_temp_init(mpz_temp * temp, unsigned int lenght);
void mpz_temp_init2(mpz_temp * temp, unsigned int lenght,
		    unsigned long size);
void mpz_temp_clear(mpz_temp * temp);
#define mpz_temp_free(temp,n) (temp)->index = ((temp)->index - n)
#define mpz_temp_space(temp) (temp->lenght - (temp->index))

#define mpz_temp_get(rop, temp)			\
do {									\
	rop = ((temp)->t + (temp)->index);	\
	((temp)->index)++;					\
} while(0)

void mpz_rep_init(mpz_rep * rep, unsigned long lenght);
void mpz_rep_init2(mpz_rep * rep, unsigned long lenght,
		   unsigned long size);
void mpz_rep_clear(mpz_rep * rep);

void get_randprime(mpz_t prime, const mpz_t offset, const mpz_t range,
		   gmp_randstate_t state);

void create_bigk(mpz_t k, const unsigned long b, mpz_temp * temp);	//create k the number product of all b-smooth numbers
void get_prime_diff(const unsigned long start, int sub_offs,
		    const unsigned long end, unsigned char v[],
		    mpz_temp * temp);

int get_vdiff_size(const unsigned long b2);

static inline unsigned long str2l(const char *s)
{
	errno = 0;
	char *endptr;
	unsigned long n = strtol(s, &endptr, 10);
	if (*endptr != '\0' || errno != 0)
		error_msg("error in str2l conversion\n");
	return n;
}

#define abs_modn(op, n)			\
do {							\
	while (mpz_sgn(op) == -1)	\
		mpz_add(op, op, n);		\
} while (0)

#define lazy_high_bit(res, op)							\
do {													\
	res = 0;											\
	for (unsigned long n = op; n != 0; n >>= 1)			\
		res++;											\
	res--;												\
} while (0)

#define stupid_log10_inf(res, op)						\
do {													\
	res = 0;											\
	for (unsigned long n = op; n > 9; n = n / 10)		\
		res++;											\
} while (0)

#endif //MATBASE_LIB_H_
