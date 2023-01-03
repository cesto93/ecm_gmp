#include "m_ellcurv_lib.h"

#include "matbase_lib.h"
#include "base_lib.h"

#define FACT_REP_SIZE 950

typedef struct m_fact_param {
	unsigned long b1;
	unsigned long b2;
	unsigned long max_iter;
	mpz_t n;
	mpz_t k;
	mpz_t e_C2;
	unsigned char *vdiff;
} m_fact_param;

//when used first check temp space
static inline void m_ell_addh2(const mpz_t n, m_ellp * r, const m_ellp * p, const m_ellp * q, const m_ellp * diff, mpz_temp * temp);
static inline void m_ell_duph2(const mpz_t n, const mpz_t e_C2, m_ellp * r, const m_ellp * p, mpz_temp * temp);
static void m_ell_mul_ui(unsigned long k, const mpz_t n, const mpz_t e_C2, m_ellp * r, const m_ellp * p, m_ellp_temp * p_temp, mpz_temp * temp);

void m_ell_addh_t(const m_ellc * e, m_ellp * r, const m_ellp * p, const m_ellp * q, const m_ellp * diff)	// calculate p + q
{
	mpz_temp temp;
	mpz_temp_init(&temp, n_temp_addh2);
	m_ell_addh2(e->n, r, p, q, diff, &temp);
	mpz_temp_clear(&temp);
}

void m_ell_duph_t(const m_ellc * e, m_ellp * r, const m_ellp * p)	// calculate 2p
{
	mpz_temp temp;
	mpz_temp_init(&temp, n_temp_duph2);
	m_ell_duph2(e->n, e->C2, r, p, &temp);
	mpz_temp_clear(&temp);
}

void m_ell_mul_t(const mpz_t k, const m_ellc * e, m_ellp * r, const m_ellp * p)
{
	mpz_temp temp;
	m_ellp_temp p_temp;
	mpz_temp_init(&temp, n_temp_mul + 1);	//max(n_temp_addh, n_temp_duph)
	m_ellp_temp_init(&p_temp, n_p_temp_mul);
	m_ell_mul(k, e->n, e->C2, r, p, &p_temp, &temp);
	mpz_temp_clear(&temp);
	m_ellp_temp_clear(&p_temp);
}

void m_ell_diff_t(m_ellp_rep * rep, mpz_rep * beta, const m_ellc * e, const m_ellp * p)
{
	mpz_temp temp;
	mpz_temp_init(&temp, n_temp_diff);	//max(n_temp_addh, n_temp_duph)
	m_ell_diff(rep, beta, e->n, e->C2, p, &temp);
	mpz_temp_clear(&temp);
}

static inline void m_ell_addh2(const mpz_t n, m_ellp * r, const m_ellp * p, const m_ellp * q, const m_ellp * diff, mpz_temp * temp)
{
	mpz_t *t1, *t2, *t3;
	mpz_temp_get(t1, temp);
	mpz_temp_get(t2, temp);
	mpz_temp_get(t3, temp);

	mpz_sub(*t1, p->X, p->Z);
	mpz_add(*t2, q->X, q->Z);
	mpz_mul(*t3, *t1, *t2);	//t3 = (x1-z1)(x2+z2)
	mpz_add(*t1, p->X, p->Z);
	mpz_sub(*t2, q->X, q->Z);
	mpz_mul(*t1, *t1, *t2);	//t1 = (x1+z1)(x2-z2)

	mpz_add(*t2, *t3, *t1);	//t2 = (x1-z1)(x2+z2) + (x1+z1)(x2-z2)
	mpz_sub(*t1, *t3, *t1);	//t1 = (x1-z1)(x2+z2) - (x1+z1)(x2-z2)

	mpz_mul(*t3, *t2, *t2);	//t3 = t2^2 
	mpz_mul(*t3, *t3, diff->Z);	// t3 := ((x1-z1)(x2+z2) + (x1+z1)(x2-z2))^2 * z-
	mpz_mul(*t2, *t1, *t1);	//t2 = t1^2
	mpz_mul(*t2, *t2, diff->X);	//t2 := ((x1-z1)(x2+z2) - (x1+z1)(x2-z2))^2 * x-

	mpz_mod(r->X, *t3, n);	//rx = t3
	mpz_mod(r->Z, *t2, n);	//rz = t3

	mpz_temp_free(temp, n_temp_addh2);
}

static inline void m_ell_duph2(const mpz_t n, const mpz_t e_C2, m_ellp * r, const m_ellp * p, mpz_temp * temp)
{
	mpz_t *t1, *t2, *t3;
	mpz_temp_get(t1, temp);
	mpz_temp_get(t2, temp);
	mpz_temp_get(t3, temp);

	mpz_add(*t1, p->X, p->Z);
	mpz_mul(*t2, *t1, *t1);	//t2 = (x1+z1)^2
	mpz_sub(*t1, p->X, p->Z);
	mpz_mul(*t3, *t1, *t1);	//t3 = (x1-z1)^2

	mpz_mul(*t1, *t2, *t3);	// t1 = (x1+z1)^2 * (x1-z1)^2
	mpz_mod(r->X, *t1, n);	//rx = t1

	mpz_sub(*t1, *t2, *t3);	//t1 = 4x1z1 = (x1+z1)^2 - (x1-z1)^2
	mpz_addmul(*t3, e_C2, *t1);	//t3 = (x1-z1)^2 + ((c+2)/4) * 4x1z1
	mpz_mul(*t2, *t3, *t1);	//t2 = 4x1z1 * ((x1-z1)^2 + ((c+2)/4) * 4x1z1)
	mpz_mod(r->Z, *t2, n);	//rz = t2

	mpz_temp_free(temp, n_temp_duph2);
}

void m_ell_mul(const mpz_t k, const mpz_t n, const mpz_t e_C2, m_ellp * r, const m_ellp * p, m_ellp_temp * p_temp, mpz_temp * temp)
{
	m_ellp *u, *t;
	m_ellp_temp_get(u, p_temp);
	m_ellp_temp_get(t, p_temp);

	m_ellp_set(u, p);
	m_ell_duph2(n, e_C2, t, p, temp);
	unsigned long high;
	high = mpz_sizeinbase(k, 2) - 1;
	DEBUG_STAMP("higher bit is in pos %lu\n", high);

	if (mpz_cmp_ui(k, 2) > 0) {
		for (unsigned long i = high - 1; i > 0; i--)
			if (mpz_tstbit(k, i)) {
				m_ell_addh2(n, u, t, u, p, temp);
				m_ell_duph2(n, e_C2, t, t, temp);
			} else {
				m_ell_addh2(n, t, u, t, p, temp);
				m_ell_duph2(n, e_C2, u, u, temp);
			}
		if (mpz_tstbit(k, 0))
			m_ell_addh2(n, r, u, t, p, temp);
		else
			m_ell_duph2(n, e_C2, r, u, temp);
	} else {
		gmp_printf("k %Zd\n", k);
		if (mpz_sgn(k) != 1)
			error_msg("error in ell_mul k = 0");
		if (mpz_cmp_ui(k, 2) == 0)
			m_ell_duph2(n, e_C2, r, p, temp);
		else if (mpz_cmp_ui(k, 1) == 0)
			m_ellp_set(r, p);
	}
	mpz_temp_free(p_temp, 2);
}

static void m_ell_mul_ui(unsigned long k, const mpz_t n, const mpz_t e_C2, m_ellp * r, const m_ellp * p, m_ellp_temp * p_temp, mpz_temp * temp)
{
	m_ellp *u, *t;
	m_ellp_temp_get(u, p_temp);
	m_ellp_temp_get(t, p_temp);

	m_ellp_set(u, p);
	m_ell_duph2(n, e_C2, t, p, temp);
	unsigned long high;
	lazy_high_bit(high, k);
	DEBUG_STAMP("higher bit is in pos %lu\n", high);

	if (k > 2) {
		for (unsigned long i = high - 1; i > 0; i--)
			if (k & (1 << i))	//k[i] bit
			{
				m_ell_addh2(n, u, t, u, p, temp);
				m_ell_duph2(n, e_C2, t, t, temp);
			} else {
				m_ell_addh2(n, t, u, t, p, temp);
				m_ell_duph2(n, e_C2, u, u, temp);
			}
		if (k & 1)	//k[0]
			m_ell_addh2(n, r, u, t, p, temp);
		else
			m_ell_duph2(n, e_C2, r, u, temp);
	} else {
		gmp_printf("k %Zd\n", k);
		if (k <= 0)
			error_msg("error in ell_mul k = 0");
		if (k == 2)
			m_ell_duph2(n, e_C2, r, p, temp);
		else if (k == 1)
			m_ellp_set(r, p);
	}
	mpz_temp_free(p_temp, 2);
}

void m_ell_fase1_ui(const unsigned long k[], int size, const m_ellc * e, m_ellp * r, const m_ellp * p, m_ellp_temp * p_temp, mpz_temp * temp)
{
	m_ellp_set(r, p);
	for (int i = 0; i < size; i++) {
		m_ell_mul_ui(k[i], e->n, e->C2, r, r, p_temp, temp);
	}
}

void check_mul_ui(unsigned long k, const m_ellc * e, m_ellp * r, const m_ellp * p, m_ellp_temp * p_temp, mpz_temp * temp)
{
	struct timespec start, end, diff;
	mpz_t *t;
	mpz_temp_get(t, temp);
	mpz_set_ui(*t, k);
	get_current_time(&start);
	m_ell_mul(*t, e->n, e->C2, r, p, p_temp, temp);
	get_current_time(&end);
	timespec_diff(&start, &end, &diff);
	printf_timespec("mul", &diff);
	gmp_printf("X %Zd Z %Zd\n", r->X, r->Z);

	get_current_time(&start);
	m_ell_mul_ui(k, e->n, e->C2, r, p, p_temp, temp);
	get_current_time(&end);
	timespec_diff(&start, &end, &diff);
	printf_timespec("mul_ui", &diff);
	gmp_printf("X %Zd Z %Zd\n", r->X, r->Z);
}

void m_ell_diff(m_ellp_rep * rep, mpz_rep * beta, const mpz_t n, const mpz_t e_C2, const m_ellp * p, mpz_temp * temp)
{
	mpz_t *t1;
	mpz_temp_get(t1, temp);

	m_ell_duph2(n, e_C2, &(rep->p[0]), p, temp);
	mpz_mul(*t1, rep->p[0].X, rep->p[0].Z);
	mpz_mod(beta->v[0], *t1, n);
	m_ell_duph2(n, e_C2, &(rep->p[1]), &(rep->p[0]), temp);
	mpz_mul(*t1, rep->p[1].X, rep->p[1].Z);
	mpz_mod(beta->v[1], *t1, n);

	for (unsigned long i = 2; i < rep->lenght; i++) {
		m_ell_addh2(n, &(rep->p[i]), &(rep->p[i - 1]), &(rep->p[0]), &(rep->p[i - 2]), temp);
		mpz_mul(*t1, rep->p[i].X, rep->p[i].Z);	//beta[i] = rep[i].X * rep[i].Z
		mpz_mod(beta->v[i], *t1, n);
	}
	mpz_temp_free(temp, 1);
}

void m_ell_fase2(mpz_t g, unsigned long b1, unsigned long b2, const mpz_t n,
		 const mpz_t e_C2, const m_ellp * p, const m_ellp_rep rep,
		 const mpz_rep beta, const unsigned char vdiff[], m_ellp_temp * p_temp, mpz_temp * temp)
{
	if (mpz_temp_space(temp) < n_temp_fase2)
		error_msg("error in m_ell_fase2 temp full\n");
	if (m_ellp_temp_space(p_temp) < n_p_temp_fase2)
		error_msg("p_temp full in m_ell_fase2");

	const unsigned long max_diff = 2 * rep.lenght;	//max_diff = 2 * size
	unsigned long diff = 0, b = b1 - 1;
	int i = 1;
	mpz_t *alfa, *t1, *t2;
	m_ellp *T, *R;

	mpz_temp_get(alfa, temp);
	mpz_temp_get(t1, temp);
	mpz_temp_get(t2, temp);

	m_ellp_temp_get(T, p_temp);
	m_ellp_temp_get(R, p_temp);

	m_ell_mul_ui(b, n, e_C2, R, p, p_temp, temp);	// R = p * b
	m_ell_mul_ui(b - max_diff, n, e_C2, T, p, p_temp, temp);	//T = p * (b - max_diff)

	diff = vdiff[0];
	for (; b < b2; b += max_diff) {
		mpz_mul(*alfa, R->X, R->Z);
		mpz_mod(*alfa, *alfa, n);	//alfa = XR * ZR modn

		for (; diff <= max_diff; diff += vdiff[i], i++) {
			mpz_sub(*t1, R->X, rep.p[(diff / 2) - 1].X);	// t1 = (XR − XS(δ)
			mpz_add(*t2, R->Z, rep.p[(diff / 2) - 1].Z);	// t2 = (ZR + ZS(δ)
			mpz_mul(*t1, *t1, *t2);	// t1 =(XR − XS(δ))*(ZR + ZS(δ) 
			mpz_sub(*t1, *t1, *alfa);	// t1 = (XR − XS(δ))(ZR + ZS(δ)) − α
			mpz_add(*t1, *t1, beta.v[(diff / 2) - 1]);	// t1 = (XR − XS(δ))(ZR + ZS(δ)) − α + β(δ)
			mpz_mul(g, g, *t1);	//g *((XR − XS(δ))(ZR + ZS(δ)) − α + β(δ)) mod n
			mpz_mod(g, g, n);

			if (vdiff[i] == 0)
				goto end;
		}
		diff = diff - max_diff;
		m_ell_addh2(n, T, R, &(rep.p[rep.lenght - 1]), T, temp);	// T = R + s[D]
		mpz_swap(T->X, R->X);	// T = R
		mpz_swap(T->Z, R->Z);	// R = R + s[d]
	}
end:
	mpz_temp_free(temp, 3);
	m_ellp_temp_free(p_temp, 2);
}

mpz_t *m_ell_fact(gmp_randstate_t state, const mpz_t n, unsigned long b1, unsigned long b2, unsigned long max_iter, unsigned long *iter,  int *fase_found)
{
	m_ellp *p, *r;
	mpz_t *gcd;
	mpz_temp temp;
	m_ellp_temp p_temp;
	mpz_rep beta;
	m_ellp_rep rep;
	
	m_fact_param param;
	mpz_t *fact = NULL;
	unsigned long n_size = mpz_size(n) * mp_bits_per_limb;
	int vdiff_size = get_vdiff_size(b2);

	if (vdiff_size == -1) {
		perror("b2 to big\n");
		return NULL;
	}	

	param.vdiff = malloc(vdiff_size);
	if (param.vdiff == NULL) {
		perror("error in m_ell_fact malloc()");
		return NULL;
	}

	mpz_temp_init2(&temp, n_temp_fact, mpz_size(n) * mp_bits_per_limb);	
	m_ellp_temp_init2(&p_temp, n_p_temp_fact, n_size);
	m_ellp_rep_init2(&rep, FACT_REP_SIZE, n_size);
	mpz_rep_init2(&beta, FACT_REP_SIZE, n_size);
	mpz_init2(param.k, bigk_size_bits(b1));
	mpz_init(param.n);
	mpz_init2(param.e_C2, n_size);

	create_bigk(param.k, b1, &temp);
	get_prime_diff(b1, 1, b2, param.vdiff, &temp);
	mpz_realloc2(param.k, mpz_size(param.k) * mp_bits_per_limb);
	param.b1 = b1;
	param.b2 = b2;
	param.max_iter = max_iter;
	mpz_set(param.n, n);

	mpz_temp_get(gcd, &temp);
	m_ellp_temp_get(p, &p_temp);
	m_ellp_temp_get(r, &p_temp);
	*fase_found = -1;

	for (*iter = 0; *iter < param.max_iter; (*iter)++) {
		if (m_ell_setrand2(param.n, param.e_C2, p, state, &temp)) {	//TODO inversion can be avoided
			if (find_div_by_gcd(*gcd, p->X, param.n)) {
				*fase_found = 0;
				break;
			}
		} else {
			m_ell_mul(param.k, param.n, param.e_C2, r, p, &p_temp, &temp);	//FASE1
			if (find_div_by_gcd(*gcd, r->Z, param.n)) {
				*fase_found = 1;
				break;
			}
			mpz_set_ui(*gcd, 1);
			m_ell_diff(&rep, &beta, param.n, param.e_C2, r, &temp);
			m_ell_fase2(*gcd, param.b1, param.b2, param.n, param.e_C2, r, rep, beta, param.vdiff, &p_temp, &temp);
			if (find_div_by_gcd(*gcd, *gcd, param.n)) {
				*fase_found = 2;
				break;
			}
		}
	}

	if (*fase_found != -1) {
		fact = malloc(sizeof(mpz_t));
		if (fact == NULL) {
			perror("error in m_ell_fact malloc()");
		} else {
			mpz_init(*fact);
			mpz_set(*fact, *gcd);
		}
	}

	mpz_temp_free(&temp, 1);
	m_ellp_temp_free(&p_temp, 2);
	mpz_clears(param.k, NULL);
	free(param.vdiff);
	mpz_clear(param.n);
	return fact;
}