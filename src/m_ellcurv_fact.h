#ifndef M_ELLCURV_FACT_H_	/* Include guard */
#define M_ELLCURV_FACT_H_

#define MM_ENABLE

#define NDEBUG

#include "m_ellcurv_struct.h"

#define _GNU_SOURCE

#define FACT_REP_SIZE 950

//return ELL_FACT_NOT_FOUND if fact not found
m_ellfact_res *factorize(const mpz_t n, unsigned long b1, unsigned long b2, unsigned long max_iter);

long mfactorize(mpz_t factors[], const mpz_t n, unsigned long b1, unsigned long b2, unsigned long max_iter);

#endif //M_ELLCURV_FACT_H
