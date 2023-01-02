#ifndef M_ELLCURV_STRUCT_H_	/* Include guard */
#define M_ELLCURV_STRUCT_H_

#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>

#include "matbase_lib.h"

typedef struct m_ellc		//y^2 = x^3 + C(σ)x^2 + x
{
	mpz_t C2;		// (c+2) / 4
	mpz_t n;		//the operation are made modn
} m_ellc;

typedef struct m_ellp		// point must respect equation y^2 = x^3 + C(σ)x^2 + x
{
	mpz_t X;		//coordinate x
	mpz_t Z;		//coordinate y
} m_ellp;

typedef struct m_ellp_temp	//temp variable container
{
	m_ellp *p;
	unsigned long lenght;
	unsigned long index;
} m_ellp_temp;

typedef struct m_ellp_rep {
	m_ellp *p;
	unsigned long lenght;
} m_ellp_rep;

#define n_temp_setrand2 4
//return 1 on error and not invertible in p->X
int m_ell_setrand2_t(const mpz_t n, mpz_t e_C2, m_ellp * p, gmp_randstate_t state);
int m_ell_setrand2(const mpz_t n, mpz_t e_C2, m_ellp * p, gmp_randstate_t state, mpz_temp * temp);

static inline void m_ellc_init(m_ellc * e)
{
	mpz_init(e->C2);
	mpz_init(e->n);
}

static inline void m_ellc_init2(m_ellc * e, mp_bitcnt_t size)
{
	mpz_init2(e->C2, size);
	mpz_init2(e->n, size);
}

static inline void m_ellc_clear(m_ellc * e)
{
	mpz_clear(e->C2);
	mpz_clear(e->n);
}

#define m_ellc_to_mform(rop, op, mdata, temp) to_mform((rop)->C2, (op)->C2, mdata, temp)
#define m_ellc_from_mform(rop, op, mdata, temp) from_mform((rop)->C2, (op)->C2, mdata, temp)

static inline void m_ellp_init(m_ellp * p)
{
	mpz_init(p->X);
	mpz_init(p->Z);
}

static inline void m_ellp_init2(m_ellp * p, mp_bitcnt_t size)
{
	mpz_init2(p->X, size);
	mpz_init2(p->Z, size);
}

static inline void m_ellp_clear(m_ellp * p)
{
	mpz_clear(p->X);
	mpz_clear(p->Z);
}

#define m_ellp_set(rop, op) \
do { \
	mpz_set((rop)->X, (op)->X); \
	mpz_set((rop)->Z, (op)->Z); \
} while (0)

#define m_ellp_to_mform(rop, op, mdata, temp) \
do { \
	to_mform((rop)->X, (op)->X, mdata, temp); \
	to_mform((rop)->Z, (op)->Z, mdata, temp); \
} while (0)

#define m_ellp_from_mform(rop, op, mdata, temp) \
do { \
	from_mform((rop)->X, (op)->X, mdata, temp); \
	from_mform((rop)->Z, (op)->Z, mdata, temp); \
} while (0)

static inline void m_ellp_temp_init(m_ellp_temp * temp, unsigned long lenght)
{
	temp->p = allocate(sizeof(m_ellp) * lenght);
	temp->index = 0;
	temp->lenght = lenght;
	for (unsigned int i = 0; i < temp->lenght; i++)
		m_ellp_init(&(temp->p[i]));
}

static inline void m_ellp_temp_init2(m_ellp_temp * temp, unsigned long lenght, mp_bitcnt_t size)
{
	temp->index = 0;
	temp->p = allocate(sizeof(m_ellp) * lenght);
	temp->lenght = lenght;
	for (unsigned int i = 0; i < temp->lenght; i++)
		m_ellp_init2(&(temp->p[i]), size);
}

static inline void m_ellp_temp_clear(m_ellp_temp * temp)
{
	for (unsigned int i = 0; i < temp->lenght; i++)
		m_ellp_clear(&(temp->p[i]));
	temp->lenght = 0;
	free(temp->p);
}

#define m_ellp_temp_free(temp,n) (temp)->index = ((temp)->index - n)

#define m_ellp_temp_space(temp) (temp->lenght - (temp->index))

#define m_ellp_temp_get(rop, temp) \
do { \
	rop = ((temp)->p + (temp)->index); \
	((temp)->index)++; \
} while (0)

static inline void m_ellp_rep_init(m_ellp_rep * rep, unsigned long lenght)
{
	rep->p = allocate(sizeof(m_ellp) * lenght);
	rep->lenght = lenght;
	for (unsigned int i = 0; i < rep->lenght; i++)
		m_ellp_init(&(rep->p[i]));
}

static inline void m_ellp_rep_init2(m_ellp_rep * rep, unsigned long lenght, mp_bitcnt_t size)
{
	rep->p = allocate(sizeof(m_ellp) * lenght);
	rep->lenght = lenght;
	for (unsigned int i = 0; i < rep->lenght; i++)
		m_ellp_init2(&(rep->p[i]), size);
}

static inline void m_ellp_rep_clear(m_ellp_rep * rep)
{
	for (unsigned int i = 0; i < rep->lenght; i++)
		m_ellp_clear(&(rep->p[i]));
	rep->lenght = 0;
	free(rep->p);
}

#endif //M_ELLCURV_STRUCT_H
