#ifndef M_ELLCURV_FACT_H_	/* Include guard */
#define M_ELLCURV_FACT_H_

#define MM_ENABLE

#define NDEBUG

#include "m_ellcurv_struct.h"

#define _GNU_SOURCE

#define FACT_REP_SIZE 950

#define ELL_FACT_NOT_FOUND -1
#define ell_fact_FASE(n, max_iter) (n / max_iter)
#define ell_fact_ITER(n, max_iter) (n % max_iter)

//return ELL_FACT_NOT_FOUND if fact not found
long factorize(mpz_t factors[], const mpz_t n, unsigned long b1, unsigned long b2, unsigned long max_iter);

long mfactorize(mpz_t factors[], const mpz_t n, unsigned long b1, unsigned long b2, unsigned long max_iter);

#endif //M_ELLCURV_FACT_H
