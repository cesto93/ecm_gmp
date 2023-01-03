#include "matbase_lib.h"

void get_randprime(mpz_t prime, const mpz_t offset, const mpz_t range, gmp_randstate_t state)
{
	mpz_urandomm(prime, state, range);
	mpz_add(prime, prime, offset);
	mpz_nextprime(prime, prime);
}

void create_bigk(mpz_t k, const unsigned long b, mpz_temp * temp)
{
	mpz_t *i;		// iterator of primes
	mpz_t *pow;		// pow = (i ^ k) <= B
	if (mpz_temp_space(temp) < n_temp_create_bigk)
		error_msg("error in temp_get at create_bigk\n");
	mpz_temp_get(i, temp);
	mpz_temp_get(pow, temp);
	mpz_set_ui(*i, 2);
	mpz_set_ui(k, 1);
	while (mpz_cmp_ui(*i, b) != 1)	// while (i <= b)
	{
		mpz_set(*pow, *i);
		while (mpz_cmp_ui(*pow, b) != 1)	// while pow <= B
		{
			mpz_mul(*pow, *pow, *i);
			mpz_mul(k, k, *i);
		}
		mpz_nextprime(*i, *i);
	}
	mpz_temp_free(temp, n_temp_create_bigk);
}

void get_prime_diff(const unsigned long start, int sub_offs, const unsigned long end, unsigned char v[], mpz_temp * temp)
{
	unsigned long i = 0;
	mpz_t *prev, *next, *diff;	// iterator of primes
	if (mpz_temp_space(temp) < n_temp_get_prime_diff)
		error_msg("error in temp_get at ell_curv_delta\n");
	mpz_temp_get(prev, temp);
	mpz_temp_get(next, temp);
	mpz_temp_get(diff, temp);

	mpz_set_ui(*prev, start - sub_offs);

	for (mpz_nextprime(*next, *prev); (mpz_cmp_ui(*next, end) <= 0); mpz_nextprime(*next, *next)) {
		mpz_sub(*diff, *next, *prev);
		v[i] = (char)mpz_get_ui(*diff);
		i++;
		mpz_set(*prev, *next);
	}
	v[i] = 0;
	mpz_temp_free(temp, n_temp_get_prime_diff);
}

int get_vdiff_size(const unsigned long b2)
{
	const unsigned int v_size[] = { 164, 1204, 9424, 77269, 654987, 5682957 };
	unsigned int size;
	stupid_log10_inf(size, b2);
	size++;
	size = size - 3;
	STAMP("digits %d vdiff_size = %d\n", size, v_size[size]);
	if (size < 6)
		return v_size[size];
	return -1;
}
