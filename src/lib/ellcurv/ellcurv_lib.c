#include "ellcurv_lib.h"

static int ell_add(const ell_c e, ell_p * r, ell_p p, ell_p q, mpz_temp * temp);	// calculate p + q
static int ell_dup(const ell_c e, ell_p * r, ell_p p, mpz_temp * temp);
static int ell_c_delta(const ell_c e, mpz_t ret, mpz_temp * temp);

int ell_c_delta_t(const ell_c e, mpz_t ret)
{
	int r;
	mpz_temp t;
	mpz_temp_init(&t, 5);
	r = ell_c_delta(e, ret, &t);
	mpz_temp_clear(&t);
	return r;
}

void ell_p_setrandom(ell_p * p, const mpz_t n, gmp_randstate_t state)
{
	mpz_urandomm(p->x, state, n);
	mpz_urandomm(p->y, state, n);
}

void ell_c_setrandom(ell_c * e, const ell_p p, gmp_randstate_t state, mpz_temp * temp)	//must grant delta != 0 & point exist
{
	mpz_t *t1;
	if (mpz_temp_space(temp) < 1)
		error_msg("error in temp_get at ell_curv_delta\n");
	mpz_temp_get(t1, temp);

	mpz_urandomm(e->A, state, e->n);
	mpz_mul(e->B, p.y, p.y);
	mpz_mod(e->B, e->B, e->n);	// b= y^2
	mpz_mul(*t1, p.x, p.x);
	mpz_add(*t1, *t1, e->A);	// t1 = x^2 + a
	mpz_submul(e->B, *t1, p.x);	// b= y^2 -x(x^2 + a)
	mpz_mod(e->B, e->B, e->n);

	mpz_temp_free(temp, 1);
}

int ell_add_t(ell_c e, ell_p * r, ell_p p, ell_p q)
{
	int i;
	mpz_temp t;
	mpz_temp_init(&t, 5);
	i = ell_add(e, r, p, q, &t);
	mpz_temp_clear(&t);
	return i;
}

int ell_dup_t(ell_c e, ell_p * r, ell_p p)
{
	int i;
	mpz_temp t;
	mpz_temp_init(&t, 5);
	i = ell_dup(e, r, p, &t);
	mpz_temp_clear(&t);
	return i;
}

static int ell_c_delta(const ell_c e, mpz_t ret, mpz_temp * temp)
{
	mpz_t *t1, *t2;
	if (mpz_temp_space(temp) < 2)
		error_msg("error in temp_get at ell_curv_delta\n");
	mpz_temp_get(t1, temp);
	mpz_temp_get(t2, temp);

	mpz_mul(*t1, e.A, e.A);
	mpz_mod(*t1, *t1, e.n);
	mpz_mul(*t1, *t1, e.A);
	mpz_mod(*t1, *t1, e.n);	// t1 = A^3
	mpz_mul(*t2, e.B, e.B);
	mpz_mod(*t2, *t2, e.n);	// t2 = B^2
	mpz_mul_si(ret, *t1, 4);	// ret = 4 * A^3
	//mpz_mul_2exp(ret, *t1, 2); 
	mpz_addmul_ui(ret, *t2, 27);	// ret 4 * A^3 + 27 * B^2
	mpz_mod(ret, ret, e.n);

	mpz_temp_free(temp, 2);
	return mpz_is_zero(ret);	//check delta!=0 mod n
}

static int ell_add(const ell_c e, ell_p * r, const ell_p p, const ell_p q, mpz_temp * temp)	// calculate p + q
{
	mpz_t *t1, *t2, *t3;
	if (mpz_temp_space(temp) < 3)
		error_msg("error in temp_get at ell_add\n");
	mpz_temp_get(t1, temp);
	mpz_temp_get(t2, temp);
	mpz_temp_get(t3, temp);

	mpz_sub(*t1, p.x, q.x);	// t1 = (xp -xq)
	mpz_sub(*t2, p.y, q.y);	// t2 = (yp - yq)
	mpz_mod(*t1, *t1, e.n);

	if (mpz_invert(*t3, *t1, e.n) == 0)	// t3 = 1 / (xp-xq)
	{
		mpz_set(r->x, *t1);
		r->is_inf = 1;

		mpz_temp_free(temp, 3);
		return 1;
	}
	mpz_mul(*t2, *t2, *t3);	// t2 = delta = (yp - yq) / (xp - xq)
	mpz_mod(*t2, *t2, e.n);	// mod

	mpz_mul(*t3, *t2, *t2);	// t3 = delta^2 
	mpz_sub(*t3, *t3, p.x);
	mpz_sub(*t3, *t3, q.x);	// t3 = delta^2 -xp -xq
	mpz_mod(*t3, *t3, e.n);	// t3 = xr

	mpz_sub(*t1, p.x, *t3);	// t1 = (xp- xr)
	mpz_mul(*t1, *t1, *t2);	// t1 = (xp - xr) * delta
	mpz_sub(*t1, *t1, p.y);	// t1 = (xp- xr) * delta - yp

	mpz_mod(r->y, *t1, e.n);	// y = (xp- xr) * delta - yp
	mpz_set(r->x, *t3);	// x = delta^2 -xp -xq

	mpz_temp_free(temp, 3);
	return 0;
}

static int ell_dup(const ell_c e, ell_p * r, const ell_p p, mpz_temp * temp)
{
	mpz_t *t1, *t2, *t3;
	if (mpz_temp_space(temp) < 3)
		error_msg("error in temp_get at ell_add\n");
	mpz_temp_get(t1, temp);
	mpz_temp_get(t2, temp);
	mpz_temp_get(t3, temp);

	mpz_mul(*t1, p.x, p.x);	//t1 = xp^2
	mpz_mod(*t1, *t1, e.n);
	mpz_mul_ui(*t1, *t1, 3);	//t1 = 3xp^2
	mpz_add(*t1, *t1, e.A);	//t1 = 3xp^2 + A
	mpz_mul_ui(*t2, p.y, 2);	//t2 = 2yp
	mpz_mod(*t2, *t2, e.n);

	if (mpz_invert(*t2, *t2, e.n) == 0)	// t2 = 1/ 2yp
	{
		r->is_inf = 1;
		mpz_set(r->x, *t2);

		mpz_temp_free(temp, 3);
		return 1;
	}
	mpz_mul(*t1, *t1, *t2);	//t1 = delta = (3xp^2 + A)/2y
	mpz_mod(*t1, *t1, e.n);
	mpz_mul(*t3, *t1, *t1);	// t3 = delta ^2
	mpz_submul_ui(*t3, p.x, 2);	// t3 = delta^2 - 2xp
	mpz_mod(*t3, *t3, e.n);	// t3 = xr

	mpz_sub(*t2, p.x, *t3);	// t2 = xp - xr
	mpz_mul(*t2, *t1, *t2);	// t2 = delta(xp-xr)
	mpz_sub(*t2, *t2, p.y);	// t2 = delta(xp-xr) -yp

	mpz_mod(r->y, *t2, e.n);	// ry = delta(xp-xr) -yp
	mpz_set(r->x, *t3);	// rx = delta ^ 2 -2xp

	mpz_temp_free(temp, 3);
	return 0;
}

int ell_c_fact(mpz_t factors[], const mpz_t n, const mpz_t b1, const mpz_t b2,
	       gmp_randstate_t state)
{
#define p_temp_size 2
#define temp_size 6
	unsigned long b2_size, n_size;
	int rep_size;

	mpz_t k, delta;
	ell_c e;
	ell_p p, r;
	mpz_temp temp;
	ell_rep rep, p_temp;

	if ((rep_size = get_rep_size(b2)) == -1)
		error_msg("b2 too big\n");

	b2_size = mpz_size(b2) * mp_bits_per_limb;
	n_size = mpz_size(n) * mp_bits_per_limb;

	mpz_init2(k, b2_size * 2);
	mpz_init2(delta, n_size);
	ell_c_init2(&e, n_size);
	ell_p_init2(&p, n_size);
	ell_p_init2(&r, n_size);
	mpz_temp_init2(&temp, temp_size, n_size);
	ell_rep_init(&rep, rep_size);
	ell_rep_init(&p_temp, p_temp_size);

	mpz_set_ui(factors[0], 0);
	mpz_set_ui(factors[1], 0);
	mpz_set(e.n, n);
	create_bigk(k, b1, &temp);

	int i, fase_found = -1;
	for (i = 0; i < MAX_ITER; i++) {
		ell_p_setrandom(&p, e.n, state);
		ell_c_setrandom(&e, p, state, &temp);
		GMP_DEBUG_STAMP
		    ("testing p: (%Zd, %Zd)\non curv: Y^2 = X^3 +%ZdX + %Zd\n",
		     p.x, p.y, e.A, e.B);
		STAMP("testing curve %d\n", i);

		if (ell_c_delta(e, delta, &temp))	//check delta!=0 mod n
			printf("curve not valid becouse delta must be != 0\n");
		else {
			mpz_gcd(delta, delta, e.n);
			if (mpz_cmp_ui(delta, 1) != 0) {
				GMP_DEBUG_STAMP("gcd(delta, %Zd ) is %Zd\n",
						e.n, delta);
				fase_found = 0;
				mpz_set(factors[0], delta);
				break;
			}
			if (ell_p_mul(factors[0], e, &r, p, k, &temp, &p_temp)) {
				fase_found = 1;
				break;
			}
			if (ell_p_fact_fase_2
			    (factors[0], e, r, b1, b2, &temp, &p_temp, &rep)) {
				fase_found = 2;
				break;
			}
		}
	}
	if (!mpz_is_zero(factors[0]))
		mpz_divexact(factors[1], e.n, factors[0]);

	mpz_clears(k, delta, NULL);
	ell_c_clear(&e);
	ell_p_clear(&p);
	ell_p_clear(&r);
	mpz_temp_clear(&temp);
	ell_rep_clear(&rep);
	ell_rep_clear(&p_temp);

	if (i == MAX_ITER)
		return ELL_FACT_NOT_FOUND;
	return i + MAX_ITER * fase_found;
}

int ell_p_mul(mpz_t factor, const ell_c e, ell_p * r, const ell_p p,
	      const mpz_t k, mpz_temp * temp, ell_rep * p_temp)
{
	ell_p *pow = (p_temp->p);
	unsigned long exp = 0, bit_count = 0, bit_k;
	mpz_set_ui(factor, 0);
	bit_k = mpz_popcount(k);
	ell_p_set(pow, p);

	while (mpz_tstbit(k, exp) != 1) {
		if (ell_dup(e, pow, *pow, temp))	// p = 2 * p
		{
			mpz_set(factor, pow->x);
			goto end;
		}
		exp++;
	}			//when we exit we have p = p0 * 2 ^ exp   where exp is the lowest bit of k
	ell_p_set(r, *pow);	// r = p0 * 2 ^ exp
	bit_count++;

	while (bit_count < bit_k) {
		exp++;
		if (ell_dup(e, pow, *pow, temp))	// p = p0 * 2 ^ exp
		{
			mpz_set(factor, pow->x);
			goto end;
		}
		if (mpz_tstbit(k, exp) == 1) {
			if (ell_add(e, r, *r, *pow, temp))	// r = r + p
			{
				mpz_set(factor, r->x);
				goto end;
			}
			bit_count++;
		}
	}
end:
	if (!mpz_is_zero(factor) && (mpz_cmp(factor, e.n) != 0)) {
		mpz_gcd(factor, factor, e.n);
		return 1;
	}
	return 0;
}

int create_r_difference(ell_rep * rep, mpz_t factor, const ell_c e,
			const ell_p p, unsigned long start,
			unsigned long end, mpz_temp * temp)
{
	unsigned long empty_i;
	if (start <= 4) {
		if (ell_dup(e, &(rep->p[0]), p, temp)) {
			mpz_set(factor, rep->p[0].x);
			goto end;
		}
		if (ell_dup(e, &(rep->p[1]), rep->p[0], temp)) {
			mpz_set(factor, rep->p[1].x);
			goto end;
		}
		empty_i = 2;
	} else
		empty_i = start / 2;

	for (; empty_i < end / 2; empty_i++)	//end va calcolato
	{
		if (ell_add
		    (e, &(rep->p[empty_i]), rep->p[empty_i - 1], rep->p[0],
		     temp)) {
			mpz_set(factor, rep->p[empty_i].x);
			goto end;
		}
	}
end:
	if (!mpz_is_zero(factor) && (mpz_cmp(factor, e.n) != 0)) {
		mpz_gcd(factor, factor, e.n);
		return 1;
	}
	return 0;
}

int ell_p_fact_fase_2(mpz_t factor, const ell_c e, const ell_p r,
		      const mpz_t b1, const mpz_t b2, mpz_temp * temp,
		      ell_rep * p_temp, ell_rep * rep)
{
	ell_p *q = &p_temp->p[1];	// q = qi * r 
	mpz_t *prime, *old, *diff;
	unsigned long diff_ui, max_diff = 0;

	if (mpz_temp_space(temp) < 3)
		error_msg("error in temp_get at ell_point_fact2\n");
	mpz_temp_get(prime, temp);
	mpz_temp_get(old, temp);
	mpz_temp_get(diff, temp);

	mpz_set_ui(factor, 0);
	mpz_sub_ui(*prime, b1, 1);
	mpz_nextprime(*prime, *prime);	//prime = q0

	if (!ell_p_mul(factor, e, q, r, *prime, temp, p_temp))	// q = r * qo
		while (mpz_cmp(*prime, b2) != 1) {
			mpz_set(*old, *prime);
			mpz_nextprime(*prime, *prime);
			mpz_sub(*diff, *prime, *old);
			diff_ui = mpz_get_ui(*diff);

			if (diff_ui > max_diff) {
				if (create_r_difference
				    (rep, factor, e, r, max_diff, diff_ui + 4,
				     temp))
					goto end;
				max_diff = diff_ui + 4;
			}

			if (ell_add(e, q, rep->p[(diff_ui / 2) - 1], *q, temp))	// q = r * q1
			{
				mpz_set(factor, q->x);
				break;
			}
		}
end:
	mpz_temp_free(temp, 3);

	if (!mpz_is_zero(factor) && (mpz_cmp(factor, e.n) != 0)) {
		mpz_gcd(factor, factor, e.n);
		return 1;
	}
	return 0;
}

int get_rep_size(const mpz_t b2)
{
	size_t size = mpz_sizeinbase(b2, 10) - 3;
	const unsigned int array_size[] = { 10, 18, 36, 57, 77, 111 };
	DEBUG_STAMP("rep_size = %d\n", array_size[size]);
	if (size < 6)
		return array_size[size] + 4;
	else
		return -1;
}
