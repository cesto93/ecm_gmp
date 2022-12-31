#include "ellcurv_struct.h"

void ell_c_init(ell_c * e)
{
	mpz_init(e->A);
	mpz_init(e->B);
	mpz_init(e->n);
}

void ell_c_init2(ell_c * e, unsigned long size)
{
	mpz_init2(e->A, size);
	mpz_init2(e->B, size);
	mpz_init2(e->n, size);
}

void ell_c_clear(ell_c * e)
{
	mpz_clear(e->A);
	mpz_clear(e->B);
	mpz_clear(e->n);
}

void ell_p_init(ell_p * p)
{
	mpz_init(p->x);
	mpz_init(p->y);
	p->is_inf = 0;
}

void ell_p_init2(ell_p * p, unsigned long size)
{
	mpz_init2(p->x, size);
	mpz_init2(p->y, size);
	p->is_inf = 0;
}

void ell_p_inits(ell_p * p, unsigned int size)
{
	for (unsigned int i = 0; i < size; i++) {
		mpz_init(p[i].x);
		mpz_init(p[i].y);
		p->is_inf = 0;
	}
}

void ell_p_clear(ell_p * p)
{
	mpz_clear(p->x);
	mpz_clear(p->y);
}

void ell_p_clears(ell_p * p, unsigned int size)
{
	for (unsigned int i = 0; i < size; i++) {
		mpz_clear(p[i].x);
		mpz_clear(p[i].y);
	}
}

void ell_rep_init(ell_rep * rep, unsigned int size)
{
	rep->p = allocate(sizeof(ell_p) * size);
	rep->size = size;
	for (unsigned int i = 0; i < size; i++) {
		mpz_init(rep->p[i].x);
		mpz_init(rep->p[i].y);
		rep->p[i].is_inf = 0;
	}
}

void ell_rep_clear(ell_rep * rep)
{
	for (unsigned int i = 0; i < rep->size; i++) {
		mpz_clear(rep->p[i].x);
		mpz_clear(rep->p[i].y);
	}
	rep->size = 0;
	free(rep->p);
}
