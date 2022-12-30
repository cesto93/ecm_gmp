#include "m_ellcurv_struct.h"

static inline void m_ellp_init_base(m_ellp * p)
{
	mpz_init(p->X);
	mpz_init(p->Z);
}

static inline void m_ellp_init2_base(m_ellp * p, mp_bitcnt_t size)
{
	mpz_init2(p->X, size);
	mpz_init2(p->Z, size);
}

static inline void m_ellp_clear_base(m_ellp * p)
{
	mpz_clear(p->X);
	mpz_clear(p->Z);
}

void m_ellc_init(m_ellc * e)
{
	mpz_init(e->C2);
	mpz_init(e->n);
}

void m_ellc_init2(m_ellc * e, mp_bitcnt_t size)
{
	mpz_init2(e->C2, size);
	mpz_init2(e->n, size);
}

void m_ellc_clear(m_ellc * e)
{
	mpz_clear(e->C2);
	mpz_clear(e->n);
}

void m_ellp_init(m_ellp * p)
{
	m_ellp_init_base(p);
}

void m_ellp_init2(m_ellp * p, mp_bitcnt_t size)
{
	m_ellp_init2_base(p, size);
}

void m_ellp_clear(m_ellp * p)
{
	m_ellp_clear_base(p);
}

void m_ellp_temp_init(m_ellp_temp * temp, unsigned long int lenght)
{
	temp->p = allocate(sizeof(m_ellp) * lenght);
	temp->index = 0;
	temp->lenght = lenght;
	unsigned int i;
	for (i = 0; i < temp->lenght; i++)
		m_ellp_init_base(&(temp->p[i]));
}

void m_ellp_temp_init2(m_ellp_temp * temp, unsigned long int lenght,
		       mp_bitcnt_t size)
{
	temp->index = 0;
	unsigned int i;
	temp->p = allocate(sizeof(m_ellp) * lenght);
	temp->lenght = lenght;
	for (i = 0; i < temp->lenght; i++)
		m_ellp_init2_base(&(temp->p[i]), size);
}

void m_ellp_temp_clear(m_ellp_temp * temp)
{
	unsigned int i;
	for (i = 0; i < temp->lenght; i++)
		m_ellp_clear_base(&(temp->p[i]));
	temp->lenght = 0;
	free(temp->p);
}

void m_ellp_rep_init(m_ellp_rep * rep, unsigned long int lenght)
{
	unsigned int i;
	rep->p = allocate(sizeof(m_ellp) * lenght);
	rep->lenght = lenght;
	for (i = 0; i < rep->lenght; i++)
		m_ellp_init_base(&(rep->p[i]));
}

void m_ellp_rep_init2(m_ellp_rep * rep, unsigned long int lenght,
		      mp_bitcnt_t size)
{
	unsigned int i;
	rep->p = allocate(sizeof(m_ellp) * lenght);
	rep->lenght = lenght;
	for (i = 0; i < rep->lenght; i++)
		m_ellp_init2_base(&(rep->p[i]), size);
}

void m_ellp_rep_clear(m_ellp_rep * rep)
{
	unsigned int i;
	for (i = 0; i < rep->lenght; i++)
		m_ellp_clear_base(&(rep->p[i]));
	rep->lenght = 0;
	free(rep->p);
}

static inline void m_ell_setrand_uv(gmp_randstate_t state, mpz_t sigma, mpz_t u,
				    mpz_t v, const mpz_t n)
{
	mpz_sub_ui(u, n, 6);	//t_u = n-6
	mpz_urandomm(sigma, state, u);	// σ [6, n-1]
	mpz_add_ui(sigma, sigma, 6);	//t3 = σ

	mpz_mul(u, sigma, sigma);
	mpz_sub_ui(u, u, 5);
	mpz_mod(u, u, n);	//u = (σ^2 − 5) modn
	mpz_mul_ui(v, sigma, 4);
	mpz_mod(v, v, n);	//v = 4σ modn
	GMP_DEBUG_STAMP("sigma = %Zd\nu = %Zd\nv = %Zd\n", sigma, u, v);
}

int m_ell_setrand2(const mpz_t n, mpz_t e_C2, m_ellp * p, gmp_randstate_t state,
		   mpz_temp * temp)
{
	mpz_t *t1, *t2, *t3, *t4;
	mpz_temp_get(t1, temp);
	mpz_temp_get(t2, temp);
	mpz_temp_get(t3, temp);
	mpz_temp_get(t4, temp);

	m_ell_setrand_uv(state, *t3, *t1, *t2, n);	// t1 = u    t2 = v t3 = sigma

	mpz_mul(*t3, *t1, *t1);
	mpz_mul(*t3, *t3, *t1);
	mpz_mod(p->X, *t3, n);	//px = u^3 FINAL

	mpz_mul(*t3, *t2, *t2);
	mpz_mul(*t3, *t3, *t2);
	mpz_mod(p->Z, *t3, n);	//py = v^3 FINAL 

	mpz_sub(*t3, *t2, *t1);	//t3 = v-u
	mpz_mul(*t4, *t3, *t3);
	mpz_mul(*t4, *t4, *t3);	//t4 = [v-u] ^3
	mpz_mod(*t3, *t4, n);	//t3 = [v-u]^3 mod

	mpz_mul_ui(*t1, *t1, 3);	//t1 = 3u   //LAST U
	mpz_add(*t1, *t1, *t2);	//t1 = [3u] + v
	mpz_mul(*t3, *t3, *t1);	//t3 = [(v-u)^3] * [3u + v]

	mpz_mul(*t1, p->X, *t2);	//t1 = [u^3] *v //LAST V
	mpz_mul_ui(*t2, *t1, 4);
	mpz_mod(*t2, *t2, n);	//t2 = 4 * [u^3 *v]

	if (mpz_invert(*t2, *t2, n) == 0)	//t2 = 1/ [4u^3 * v]
	{
		mpz_set(p->X, *t2);
		mpz_temp_free(temp, n_temp_setrand2);
		return 1;
	}
	mpz_mul(*t3, *t3, *t2);	//t3 = c+2 = [(v − u)^3 * (3u + v)] * [1/(4u^3 * v)]
	/*mpz_sub_ui(*t4, *t3, 2); //t4 = c
	   mpz_mod(e->C, *t4, e->n); //C = [(v − u)^3 * (3u + v)/(4u^3 * v)] − 2 mod n */

	mpz_mul(*t1, *t1, *t2);	//t1 = 1/4
	mpz_mul(*t2, *t3, *t1);	//t2 = [c+2] *  [1/ 4]
	mpz_mod(e_C2, *t2, n);	//c2 = t2

	mpz_temp_free(temp, 4);
	return 0;
}

int m_ell_setrand2_t(const mpz_t n, mpz_t e_C2, m_ellp * p, gmp_randstate_t state)	//return 1 on error and not invertible in p->X
{
	int i;
	mpz_temp temp;
	mpz_temp_init(&temp, n_temp_setrand2);
	i = m_ell_setrand2(n, e_C2, p, state, &temp);
	mpz_temp_clear(&temp);
	return i;
}
