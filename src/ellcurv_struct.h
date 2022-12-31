#ifndef ELLCURV_STRUCT_H_	/* Include guard */
#define ELLCURV_STRUCT_H_

#include <gmp.h>

#include "matbase_lib.h"

#define ell_p_is_inf(p) (p.is_inf == 1)

#define ell_p_set(rop, op) \
		{ \
			mpz_set(rop->x, (op).x); \
			mpz_set(rop->y, (op).y); \
		}

typedef struct ell_c {
	mpz_t A;		//the x coefficient 
	mpz_t B;		//the note term
	mpz_t n;		//the operation are made modn
} ell_c;

typedef struct ell_p		// point must respect equation y^2 = x^3 + Ax + b
{
	mpz_t x;		//coordinate x
	mpz_t y;		//coordinate y
	short int is_inf;
} ell_p;

typedef struct ell_rep		// p = n * r
{
	ell_p *p;
	unsigned int size;
} ell_rep;

void ell_c_init(ell_c * e);
void ell_c_init2(ell_c * e, unsigned long size);
void ell_c_clear(ell_c * e);

void ell_p_init(ell_p * p);
void ell_p_init2(ell_p * p, unsigned long size);
void ell_p_inits(ell_p * p, unsigned int size);
void ell_p_clear(ell_p * p);
void ell_p_clears(ell_p * p, unsigned int size);

void ell_rep_init(ell_rep * rep, unsigned int size);
void ell_rep_clear(ell_rep * rep);

#endif //ELLCURV_STRUCT_H
