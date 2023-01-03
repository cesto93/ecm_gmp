#include "mm_ellcurv_lib.h"

#include <pthread.h>

typedef struct mm_fact_param {
	unsigned long b1;
	unsigned long b2;
	unsigned long max_iter;
	mform_data mdata;
	mpz_t k;
	unsigned char *vdiff;
	mpz_t e_C2;
} mm_fact_param;

static inline void mm_ell_addh_l(mp_limb_t * rx_l, mp_limb_t * rz_l,
				 const mp_limb_t * px_l, size_t px_s,
				 const mp_limb_t * pz_l, size_t pz_s,
				 const mp_limb_t * qx_l, size_t qx_s,
				 const mp_limb_t * qz_l, size_t qz_s,
				 const mp_limb_t * dx_l, size_t dx_s,
				 const mp_limb_t * dz_l, size_t dz_s,
				 const mform_data * mdata, mp_limb_t * t1_l, mp_limb_t * t2_l, mp_limb_t * t3_l, mp_limb_t * tmul_l);

static inline void mm_ell_duph_l(const mp_limb_t * e_C2_l, size_t e_C2_s,
				 mp_limb_t * rx_l, mp_limb_t * rz_l,
				 const mp_limb_t * px_l, size_t px_s,
				 const mp_limb_t * pz_l, size_t pz_s,
				 const mform_data * mdata, mp_limb_t * t1_l, mp_limb_t * t2_l, mp_limb_t * t3_l, mp_limb_t * tmul_l);

static inline void mm_ell_addh_l_n(mp_limb_t * rx_l, mp_limb_t * rz_l,
				   const mp_limb_t * px_l,
				   const mp_limb_t * pz_l,
				   const mp_limb_t * qx_l,
				   const mp_limb_t * qz_l,
				   const mp_limb_t * dx_l,
				   const mp_limb_t * dz_l, const mform_data * mdata, mp_limb_t * t1_l, mp_limb_t * t2_l, mp_limb_t * t3_l, mp_limb_t * tmul_l);

static inline void mm_ell_duph_l_n(const mp_limb_t * e_C2_l, mp_limb_t * rx_l,
				   mp_limb_t * rz_l, const mp_limb_t * px_l,
				   const mp_limb_t * pz_l, const mform_data * mdata, mp_limb_t * t1_l, mp_limb_t * t2_l, mp_limb_t * t3_l, mp_limb_t * tmul_l);

static void mm_ell_mul_ui(unsigned long k, const mpz_t e_C2, m_ellp * r, const m_ellp * p, const mform_data * mdata, m_ellp_temp * p_temp, mpz_temp * temp);

//LOW_LEVEL_MACRO

#define mm_ell_addh_t(r, p, q, d, mdata, t1_l, t2_l, t3_l, tmul_l)										\
do {																									\
	mp_limb_t *rz_l, *rx_l;																				\
	const mp_limb_t *px_l, *pz_l, *qx_l, *qz_l, *dx_l, *dz_l;											\
	size_t px_s, pz_s, qx_s, qz_s, dx_s, dz_s;															\
																										\
	rz_l = mpz_limbs_write((r)->Z, (mdata->n_s) * mp_bits_per_limb);									\
	rx_l = mpz_limbs_write((r)->X, (mdata->n_s) * mp_bits_per_limb);									\
																										\
	MPZ_READ_WITH_SIZE(px_l, px_s, (p)->X);																\
	MPZ_READ_WITH_SIZE(pz_l, pz_s, (p)->Z);																\
	MPZ_READ_WITH_SIZE(qx_l, qx_s, (q)->X);																\
	MPZ_READ_WITH_SIZE(qz_l, qz_s, (q)->Z);																\
	MPZ_READ_WITH_SIZE(dx_l, dx_s, (d)->X);																\
	MPZ_READ_WITH_SIZE(dz_l, dz_s, (d)->Z);																\
																										\
	mm_ell_addh_l(rx_l, rz_l, px_l, px_s, pz_l, pz_s, qx_l, qx_s, qz_l, qz_s, dx_l, dx_s, dz_l, dz_s, 	\
	mdata, t1_l, t2_l, t3_l, tmul_l);																	\
																										\
	mpz_limbs_finish((r)->X, mdata->n_s);																\
	mpz_limbs_finish((r)->Z, mdata->n_s);																\
} while (0)

#define mm_ell_duph_t(e_C2, r, p, mdata, t1, t2, t3, tmul)												\
do {																									\
	mp_limb_t *rx_l, *rz_l;																				\
	const mp_limb_t *px_l, *pz_l, *e_C2_l;																\
	size_t px_s, pz_s, e_C2_s;																			\
																										\
	rz_l = mpz_limbs_write((r)->Z, (mdata->n_s) * mp_bits_per_limb);									\
	rx_l = mpz_limbs_write((r)->X, (mdata->n_s) * mp_bits_per_limb);									\
																										\
	MPZ_READ_WITH_SIZE(px_l, px_s, (p)->X);																\
	MPZ_READ_WITH_SIZE(pz_l, pz_s, (p)->Z);																\
	MPZ_READ_WITH_SIZE(e_C2_l, e_C2_s, e_C2);															\
																										\
	mm_ell_duph_l(e_C2_l, e_C2_s, rx_l, rz_l, px_l, px_s, pz_l, pz_s, mdata, t1, t2, t3, tmul);			\
																										\
	mpz_limbs_finish((r)->X, mdata->n_s);																\
	mpz_limbs_finish((r)->Z, mdata->n_s);																\
} while (0)

//HIGH LEVEL FUNCTION

//USED ONCE FOR ITER
void mm_ell_mul(const mpz_t k, mpz_t e_C2, m_ellp * r, m_ellp * p, const mform_data * mdata, m_ellp_temp * p_temp, mpz_temp * temp)
{
	m_ellp *u, *t;
	mpz_t *t1, *t2, *t3, *tmul;

	mp_limb_t *t1_l, *t2_l, *t3_l, *tmul_l;
	m_ellp_temp_get(u, p_temp);
	m_ellp_temp_get(t, p_temp);

	mpz_temp_get(t1, temp);
	mpz_temp_get(t2, temp);
	mpz_temp_get(t3, temp);
	mpz_temp_get(tmul, temp);

	t1_l = mpz_limbs_write(*t1, (mdata->n_s) * mp_bits_per_limb);
	t2_l = mpz_limbs_write(*t2, (mdata->n_s) * mp_bits_per_limb);
	t3_l = mpz_limbs_write(*t3, (mdata->n_s) * mp_bits_per_limb);
	tmul_l = mpz_limbs_write(*tmul, 2 * (mdata->n_s) * mp_bits_per_limb);

	mp_limb_t *rx_l, *rz_l, *ux_l, *uz_l, *tx_l, *tz_l;
	mp_limb_t *px_l, *pz_l, *e_C2_l;	//CONST

	rx_l = mpz_limbs_write(r->X, (mdata->n_s) * mp_bits_per_limb);
	rz_l = mpz_limbs_write(r->Z, (mdata->n_s) * mp_bits_per_limb);
	ux_l = mpz_limbs_write(u->X, (mdata->n_s) * mp_bits_per_limb);
	uz_l = mpz_limbs_write(u->Z, (mdata->n_s) * mp_bits_per_limb);
	tx_l = mpz_limbs_write(t->X, (mdata->n_s) * mp_bits_per_limb);
	tz_l = mpz_limbs_write(t->Z, (mdata->n_s) * mp_bits_per_limb);

	MPZ_MODIFY_NSIZE(px_l, p->X, mdata->n_s);
	MPZ_MODIFY_NSIZE(pz_l, p->Z, mdata->n_s);
	MPZ_MODIFY_NSIZE(e_C2_l, e_C2, mdata->n_s);

	mpn_copyi(ux_l, px_l, mdata->n_s);	//SET_ELLP(u,p)
	mpn_copyi(uz_l, pz_l, mdata->n_s);

	mm_ell_duph_l_n(e_C2_l, tx_l, tz_l, px_l, pz_l, mdata, t1_l, t2_l, t3_l, tmul_l);

	unsigned long high;
	high = mpz_sizeinbase(k, 2) - 1;

	if (mpz_cmp_ui(k, 2) > 0) {
		for (unsigned long i = high - 1; i > 0; i--)
			if (mpz_tstbit(k, i)) {
				mm_ell_addh_l_n(ux_l, uz_l, tx_l, tz_l, ux_l, uz_l, px_l, pz_l, mdata, t1_l, t2_l, t3_l, tmul_l);
				mm_ell_duph_l_n(e_C2_l, tx_l, tz_l, tx_l, tz_l, mdata, t1_l, t2_l, t3_l, tmul_l);
			} else {
				mm_ell_addh_l_n(tx_l, tz_l, ux_l, uz_l, tx_l, tz_l, px_l, pz_l, mdata, t1_l, t2_l, t3_l, tmul_l);
				mm_ell_duph_l_n(e_C2_l, ux_l, uz_l, ux_l, uz_l, mdata, t1_l, t2_l, t3_l, tmul_l);
			}
		if (mpz_tstbit(k, 0))
			mm_ell_addh_l_n(rx_l, rz_l, ux_l, uz_l, tx_l, tz_l, px_l, pz_l, mdata, t1_l, t2_l, t3_l, tmul_l);
		else
			mm_ell_duph_l_n(e_C2_l, rx_l, rz_l, ux_l, uz_l, mdata, t1_l, t2_l, t3_l, tmul_l);

		mpz_limbs_finish(r->X, mdata->n_s);
		mpz_limbs_finish(r->Z, mdata->n_s);
	} else {
		gmp_printf("k %Zd\n", k);
		if (mpz_sgn(k) != 1)
			error_msg("error in ell_mul k = 0");
		if (mpz_cmp_ui(k, 2) == 0)
			mm_ell_duph_t(e_C2, r, p, mdata, t1_l, t2_l, t3_l, tmul_l);
		else if (mpz_cmp_ui(k, 1) == 0)
			m_ellp_set(r, p);
	}

	mpz_temp_free(temp, 4);
	m_ellp_temp_free(p_temp, 2);
}

//USED TWICE FOR ITER
static void mm_ell_mul_ui(unsigned long k, const mpz_t e_C2, m_ellp * r, const m_ellp * p, const mform_data * mdata, m_ellp_temp * p_temp, mpz_temp * temp)
{
	m_ellp *u, *t;
	mpz_t *t1, *t2, *t3, *tmul;

	mp_limb_t *t1_l, *t2_l, *t3_l, *tmul_l;

	mpz_temp_get(t1, temp);
	mpz_temp_get(t2, temp);
	mpz_temp_get(t3, temp);
	mpz_temp_get(tmul, temp);

	t1_l = mpz_limbs_write(*t1, (mdata->n_s) * mp_bits_per_limb);
	t2_l = mpz_limbs_write(*t2, (mdata->n_s) * mp_bits_per_limb);
	t3_l = mpz_limbs_write(*t3, (mdata->n_s) * mp_bits_per_limb);
	tmul_l = mpz_limbs_write(*tmul, 2 * (mdata->n_s) * mp_bits_per_limb);

	m_ellp_temp_get(u, p_temp);
	m_ellp_temp_get(t, p_temp);

	m_ellp_set(u, p);

	mm_ell_duph_t(e_C2, t, p, mdata, t1_l, t2_l, t3_l, tmul_l);
	unsigned long high;
	lazy_high_bit(high, k);

	if (k > 2) {
		for (unsigned long i = high - 1; i > 0; i--)
			if (k & (1 << i))	//k[i] bit
			{
				mm_ell_addh_t(u, t, u, p, mdata, t1_l, t2_l, t3_l, tmul_l);
				mm_ell_duph_t(e_C2, t, t, mdata, t1_l, t2_l, t3_l, tmul_l);
			} else {
				mm_ell_addh_t(t, u, t, p, mdata, t1_l, t2_l, t3_l, tmul_l);
				mm_ell_duph_t(e_C2, u, u, mdata, t1_l, t2_l, t3_l, tmul_l);
			}
		if (k & 1)	//k[0]
			mm_ell_addh_t(r, u, t, p, mdata, t1_l, t2_l, t3_l, tmul_l);
		else
			mm_ell_duph_t(e_C2, r, u, mdata, t1_l, t2_l, t3_l, tmul_l);
	} else {
		printf("k %lu\n", k);
		if (k <= 0)
			error_msg("error in ell_mul k = 0");
		if (k == 2)
			mm_ell_duph_t(e_C2, r, p, mdata, t1_l, t2_l, t3_l, tmul_l);
		else if (k == 1)
			m_ellp_set(r, p);
	}

	mpz_temp_free(temp, 4);
	mpz_temp_free(p_temp, 2);
}

//USED ONCE FOR ITER / VERY LIGHT
void mm_ell_diff(m_ellp * rep, mpz_t * beta, const unsigned long d, const mpz_t e_C2, const m_ellp * p, const mform_data * mdata, mpz_temp * temp)
{
	assert(mpz_temp_space(temp) >= n_temp_mdiff);

	mpz_t *t1, *t2, *t3, *tmul;

	mp_limb_t *t1_l, *t2_l, *t3_l, *tmul_l;

	mpz_temp_get(t1, temp);
	mpz_temp_get(t2, temp);
	mpz_temp_get(t3, temp);
	mpz_temp_get(tmul, temp);

	t1_l = mpz_limbs_write(*t1, (mdata->n_s) * mp_bits_per_limb);
	t2_l = mpz_limbs_write(*t2, (mdata->n_s) * mp_bits_per_limb);
	t3_l = mpz_limbs_write(*t3, (mdata->n_s) * mp_bits_per_limb);
	tmul_l = mpz_limbs_write(*tmul, 2 * (mdata->n_s) * mp_bits_per_limb);

	mm_ell_duph_t(e_C2, &(rep[0]), p, mdata, t1_l, t2_l, t3_l, tmul_l);
	mmul(beta[0], rep[0].X, rep[0].Z, tmul_l, mdata->n_l, mdata->n_inv, mdata->n_s);	// beta[0] = X0 * Z0
	mm_ell_duph_t(e_C2, &(rep[1]), &(rep[0]), mdata, t1_l, t2_l, t3_l, tmul_l);
	mmul(beta[1], rep[1].X, rep[1].Z, tmul_l, mdata->n_l, mdata->n_inv, mdata->n_s);	// beta[1] = X1 * Z1

	for (unsigned long i = 2; i < d; i++) {
		mm_ell_addh_t(&(rep[i]), &(rep[i - 1]), &(rep[0]), &(rep[i - 2]), mdata, t1_l, t2_l, t3_l, tmul_l);
		mmul(beta[i], rep[i].X, rep[i].Z, tmul_l, mdata->n_l, mdata->n_inv, mdata->n_s);	//beta[i] = rep[i].X * rep[i].Z
	}

	mpz_temp_free(temp, 4);
}

void mm_ell_fase2(mpz_t g, const unsigned long b1, const unsigned long b2,
		  const mpz_t e_C2, const m_ellp * p, m_ellp * S, mpz_t * beta,
		  const unsigned long d, const unsigned char vdiff[], const mform_data * mdata, m_ellp_temp * p_temp, mpz_temp * temp)
{
	assert(mpz_temp_space(temp) >= n_temp_mfase2);
	assert(m_ellp_temp_space(p_temp) >= n_p_temp_mfase2);

	const unsigned long max_diff = 2 * d;	//max_diff = 2 * d
	unsigned long diff = 0, b = b1 - 1;
	m_ellp *T, *R;
	mpz_t *alfa, *t1, *t2, *tmul_t, *t3_t;

	m_ellp_temp_get(T, p_temp);
	m_ellp_temp_get(R, p_temp);
	mpz_temp_get(alfa, temp);
	mpz_temp_get(t1, temp);
	mpz_temp_get(t2, temp);

	mm_ell_mul_ui(b, e_C2, R, p, mdata, p_temp, temp);	// R = p * b
	mm_ell_mul_ui(b - max_diff, e_C2, T, p, mdata, p_temp, temp);	//T = p * (b - max_diff)

	mpz_temp_get(t3_t, temp);
	mpz_temp_get(tmul_t, temp);

	mp_limb_t *t1_l, *t2_l, *t3_l, *tmul_l, *alfa_l, *g_l, *Rx_l, *Rz_l, *Tx_l, *Tz_l;
	mp_limb_t *Sx_l[d], *Sz_l[d], *beta_l[d];

	t1_l = mpz_limbs_write(*t1, (mdata->n_s) * mp_bits_per_limb);
	t2_l = mpz_limbs_write(*t2, (mdata->n_s) * mp_bits_per_limb);
	t3_l = mpz_limbs_write(*t3_t, (mdata->n_s) * mp_bits_per_limb);
	tmul_l = mpz_limbs_write(*tmul_t, 2 * (mdata->n_s) * mp_bits_per_limb);

	alfa_l = mpz_limbs_write(*alfa, (mdata->n_s) * mp_bits_per_limb);
	g_l = mpz_limbs_write(g, (mdata->n_s) * mp_bits_per_limb);

	MPZ_MODIFY_NSIZE(Rx_l, R->X, mdata->n_s);
	MPZ_MODIFY_NSIZE(Rz_l, R->Z, mdata->n_s);
	MPZ_MODIFY_NSIZE(Tx_l, T->X, mdata->n_s);
	MPZ_MODIFY_NSIZE(Tz_l, T->Z, mdata->n_s);

	for (unsigned int j = 0; j < d; j++) {
		MPZ_MODIFY_NSIZE(Sx_l[j], S[j].X, mdata->n_s);
		MPZ_MODIFY_NSIZE(Sz_l[j], S[j].Z, mdata->n_s);
		MPZ_MODIFY_NSIZE(beta_l[j], beta[j], mdata->n_s);
	}

	diff = vdiff[0];
	int i = 1;
	for (; b < b2; b += max_diff) {
		mmul_l_n(alfa_l, Rx_l, Rz_l, tmul_l, mdata->n_l, mdata->n_inv, mdata->n_s);	// α = XR * XZ

		for (; diff <= max_diff; diff += vdiff[i], i++) {
			sub_modR_l_n(t1_l, Rx_l, Sx_l[(diff / 2) - 1], mdata->n_l, mdata->n_s);	// T1 = (XR − XS(δ)
			add_modR_l_n(t2_l, Rz_l, Sz_l[(diff / 2) - 1], mdata->n_l, mdata->n_s);	// T2 = (ZR + ZS(δ)
			mmul_l_n(t1_l, t1_l, t2_l, tmul_l, mdata->n_l, mdata->n_inv, mdata->n_s);	// T1 =(XR − XS(δ))*(ZR + ZS(δ)

			sub_modR_l_n(t1_l, t1_l, alfa_l, mdata->n_l, mdata->n_s);	// T1 = (XR − XS(δ))(ZR + ZS(δ)) − α             
			add_modR_l_n(t1_l, t1_l, beta_l[(diff / 2) - 1], mdata->n_l, mdata->n_s);	// T1 = (XR − XS(δ))(ZR + ZS(δ)) − α + β(δ)
			mmul_l_n(g_l, g_l, t1_l, tmul_l, mdata->n_l, mdata->n_inv, mdata->n_s);	// G = G *((XR − XS(δ))(ZR + ZS(δ)) − α + β(δ)) mod n

			if (vdiff[i] == 0)
				goto end;
		}
		diff = diff - max_diff;

		mm_ell_addh_l_n(Tx_l, Tz_l, Rx_l, Rz_l, Sx_l[d - 1], Sz_l[d - 1], Tx_l, Tz_l, mdata, t1_l, t2_l, t3_l, tmul_l);	// T = R + s[D]
		SWAP_LIMB_PTR_W(Tx_l, Rx_l);	// T = R
		SWAP_LIMB_PTR_W(Tz_l, Rz_l);	// R = R + s[d]
	}
end:
	mpz_limbs_finish(g, mdata->n_s);

	mpz_temp_free(temp, 5);
	m_ellp_temp_free(p_temp, 2);
}

//LOW LEVEL FUNCTION / CALLED MANYYYYYYYYY TIMES

static inline void mm_ell_addh_l(mp_limb_t * rx_l, mp_limb_t * rz_l,
				 const mp_limb_t * px_l, size_t px_s,
				 const mp_limb_t * pz_l, size_t pz_s,
				 const mp_limb_t * qx_l, size_t qx_s,
				 const mp_limb_t * qz_l, size_t qz_s,
				 const mp_limb_t * dx_l, size_t dx_s,
				 const mp_limb_t * dz_l, size_t dz_s,
				 const mform_data * mdata, mp_limb_t * t1_l, mp_limb_t * t2_l, mp_limb_t * t3_l, mp_limb_t * tmul_l)
{
	sub_modR_l(t1_l, px_l, px_s, pz_l, pz_s, mdata->n_l, mdata->n_s);	//T3 = (x1-z1)(x2+z2)                                                                                                                                                                   
	add_modR_l(t2_l, qx_l, qx_s, qz_l, qz_s, mdata->n_l, mdata->n_s);	//              T1              T2
	mmul_l_n(t3_l, t1_l, t2_l, tmul_l, mdata->n_l, mdata->n_inv, mdata->n_s);

	add_modR_l(t1_l, px_l, px_s, pz_l, pz_s, mdata->n_l, mdata->n_s);	//T1 = (x1+z1)(x2-z2)
	sub_modR_l(t2_l, qx_l, qx_s, qz_l, qz_s, mdata->n_l, mdata->n_s);
	mmul_l_n(t1_l, t1_l, t2_l, tmul_l, mdata->n_l, mdata->n_inv, mdata->n_s);

	add_modR_l_n(t2_l, t3_l, t1_l, mdata->n_l, mdata->n_s);	//T2 = (x1-z1)(x2+z2) + (x1+z1)(x2-z2)
	sub_modR_l_n(t3_l, t3_l, t1_l, mdata->n_l, mdata->n_s);	//T3 = (x1-z1)(x2+z2) - (x1+z1)(x2-z2)
	msqr_l_n(t1_l, t2_l, tmul_l, mdata->n_l, mdata->n_inv, mdata->n_s);	//T1 = [(x1-z1)(x2+z2) + (x1+z1)(x2-z2)]^2
	msqr_l_n(t2_l, t3_l, tmul_l, mdata->n_l, mdata->n_inv, mdata->n_s);	//T2 = [(x1-z1)(x2+z2) - (x1+z1)(x2-z2)]^2
	mmul_l(t1_l, t1_l, mdata->n_s, dz_l, dz_s, tmul_l, mdata->n_l, mdata->n_inv, mdata->n_s);	//T1 = [((x1-z1)(x2+z2) + (x1+z1)(x2-z2))^2] * Z- = RX

	mmul_l(rz_l, t2_l, mdata->n_s, dx_l, dx_s, tmul_l, mdata->n_l, mdata->n_inv, mdata->n_s);	//RZ = [((x1-z1)(x2+z2) - (x1+z1)(x2-z2))^2] * X-

	mpn_copyi(rx_l, t1_l, mdata->n_s);	//RX = t1 
}

static inline void mm_ell_duph_l(const mp_limb_t * e_C2_l, size_t e_C2_s,
				 mp_limb_t * rx_l, mp_limb_t * rz_l,
				 const mp_limb_t * px_l, size_t px_s,
				 const mp_limb_t * pz_l, size_t pz_s,
				 const mform_data * mdata, mp_limb_t * t1_l, mp_limb_t * t2_l, mp_limb_t * t3_l, mp_limb_t * tmul_l)
{
	add_modR_l(t1_l, px_l, px_s, pz_l, pz_s, mdata->n_l, mdata->n_s);	// T2 = (x1+z1)^2                                                                                                                                                               
	msqr_l_n(t2_l, t1_l, tmul_l, mdata->n_l, mdata->n_inv, mdata->n_s);	//                      T1

	sub_modR_l(t1_l, px_l, px_s, pz_l, pz_s, mdata->n_l, mdata->n_s);	// T3 = (x1-z1)^2                                                                                                                                                       
	msqr_l_n(t3_l, t1_l, tmul_l, mdata->n_l, mdata->n_inv, mdata->n_s);	//                      T1

	mmul_l_n(rx_l, t2_l, t3_l, tmul_l, mdata->n_l, mdata->n_inv, mdata->n_s);	// RX = [(x1+z1)^2] * [(x1-z1)^2]

	sub_modR_l_n(t1_l, t2_l, t3_l, mdata->n_l, mdata->n_s);	//T1 = 4x1z1 = [(x1+z1)^2] - [(x1-z1)^2]
	mmul_l(t2_l, e_C2_l, e_C2_s, t1_l, mdata->n_s, tmul_l, mdata->n_l, mdata->n_inv, mdata->n_s);	//T2 = [c+2/4] * [4x1z1]
	add_modR_l_n(t3_l, t3_l, t2_l, mdata->n_l, mdata->n_s);	//T3 = [(x1-z1)^2] + [((c+2)/4) * 4x1z1]
	mmul_l_n(rz_l, t3_l, t1_l, tmul_l, mdata->n_l, mdata->n_inv, mdata->n_s);	//RZ = 4x1z1 * [(x1-z1)^2 + ((c+2)/4) * 4x1z1]
}

static inline void mm_ell_addh_l_n(mp_limb_t * rx_l, mp_limb_t * rz_l,
				   const mp_limb_t * px_l, const mp_limb_t * pz_l,
				   const mp_limb_t * qx_l, const mp_limb_t * qz_l,
				   const mp_limb_t * dx_l, const mp_limb_t * dz_l, 
				   const mform_data * mdata, mp_limb_t * t1_l, mp_limb_t * t2_l, mp_limb_t * t3_l, mp_limb_t * tmul_l)
{
	sub_modR_l_n(t1_l, px_l, pz_l, mdata->n_l, mdata->n_s);	//T3 = (x1-z1)(x2+z2)                                                                                                                                                                   
	add_modR_l_n(t2_l, qx_l, qz_l, mdata->n_l, mdata->n_s);	//              T1              T2
	mmul_l_n(t3_l, t1_l, t2_l, tmul_l, mdata->n_l, mdata->n_inv, mdata->n_s);

	add_modR_l_n(t1_l, px_l, pz_l, mdata->n_l, mdata->n_s);	//T1 = (x1+z1)(x2-z2)                                                                                                                                                           
	sub_modR_l_n(t2_l, qx_l, qz_l, mdata->n_l, mdata->n_s);	//              T1              T2
	mmul_l_n(t1_l, t1_l, t2_l, tmul_l, mdata->n_l, mdata->n_inv, mdata->n_s);

	add_modR_l_n(t2_l, t3_l, t1_l, mdata->n_l, mdata->n_s);	//T2 = (x1-z1)(x2+z2) + (x1+z1)(x2-z2)
	sub_modR_l_n(t3_l, t3_l, t1_l, mdata->n_l, mdata->n_s);	//T3 = (x1-z1)(x2+z2) - (x1+z1)(x2-z2)
	msqr_l_n(t1_l, t2_l, tmul_l, mdata->n_l, mdata->n_inv, mdata->n_s);	//T1 = [(x1-z1)(x2+z2) + (x1+z1)(x2-z2)]^2
	msqr_l_n(t2_l, t3_l, tmul_l, mdata->n_l, mdata->n_inv, mdata->n_s);	//T2 = [(x1-z1)(x2+z2) - (x1+z1)(x2-z2)]^2
	mmul_l_n(t1_l, t1_l, dz_l, tmul_l, mdata->n_l, mdata->n_inv, mdata->n_s);	//T1 = [((x1-z1)(x2+z2) + (x1+z1)(x2-z2))^2] * Z- = RX

	mmul_l_n(rz_l, t2_l, dx_l, tmul_l, mdata->n_l, mdata->n_inv, mdata->n_s);	//RZ = [((x1-z1)(x2+z2) - (x1+z1)(x2-z2))^2] * X-

	mpn_copyi(rx_l, t1_l, mdata->n_s);	//RX = t1 
}

static inline void mm_ell_duph_l_n(const mp_limb_t * e_C2_l, mp_limb_t * rx_l, mp_limb_t * rz_l, const mp_limb_t * px_l,
				   const mp_limb_t * pz_l, const mform_data * mdata, mp_limb_t * t1_l, mp_limb_t * t2_l, mp_limb_t * t3_l, mp_limb_t * tmul_l)
{
	add_modR_l_n(t1_l, px_l, pz_l, mdata->n_l, mdata->n_s);	// T2 = (x1+z1)^2                                                                                                                                                               
	msqr_l_n(t2_l, t1_l, tmul_l, mdata->n_l, mdata->n_inv, mdata->n_s);	//                      T1

	sub_modR_l_n(t1_l, px_l, pz_l, mdata->n_l, mdata->n_s);	// T3 = (x1-z1)^2                                                                                                                                                       
	msqr_l_n(t3_l, t1_l, tmul_l, mdata->n_l, mdata->n_inv, mdata->n_s);	//                      T1

	mmul_l_n(rx_l, t2_l, t3_l, tmul_l, mdata->n_l, mdata->n_inv, mdata->n_s);	// RX = [(x1+z1)^2] * [(x1-z1)^2]

	sub_modR_l_n(t1_l, t2_l, t3_l, mdata->n_l, mdata->n_s);	//T1 = 4x1z1 = [(x1+z1)^2] - [(x1-z1)^2]
	mmul_l_n(t2_l, e_C2_l, t1_l, tmul_l, mdata->n_l, mdata->n_inv, mdata->n_s);	//T2 = [c+2/4] * [4x1z1]
	add_modR_l_n(t3_l, t3_l, t2_l, mdata->n_l, mdata->n_s);	//T3 = [(x1-z1)^2] + [((c+2)/4) * 4x1z1]
	mmul_l_n(rz_l, t3_l, t1_l, tmul_l, mdata->n_l, mdata->n_inv, mdata->n_s);	//RZ = 4x1z1 * [(x1-z1)^2 + ((c+2)/4) * 4x1z1]
}

// WRAPPER FUNCTION FOR TESTING / DON'T USE TEMP
void mm_ell_mul_t(const mpz_t k, const mpz_t n, mpz_t e_C2, m_ellp * r, m_ellp * p)
{
	mform_data mdata;
	mpz_temp temp;
	m_ellp_temp p_temp;

	mpz_temp_init(&temp, n_temp_mmul);	//max(n_temp_addh, n_temp_duph)
	m_ellp_temp_init(&p_temp, n_p_temp_mmul);
	mform_data_init(&mdata, n, &temp);

	to_mform(e_C2, e_C2, &mdata, &temp);
	m_ellp_to_mform(p, p, &mdata, &temp);
	mm_ell_mul(k, e_C2, r, p, &mdata, &p_temp, &temp);
	m_ellp_from_mform(r, r, &mdata, &temp);
	m_ellp_from_mform(p, p, &mdata, &temp);
	from_mform(e_C2, e_C2, &mdata, &temp);

	mform_data_clear(&mdata);
	mpz_temp_clear(&temp);
	m_ellp_temp_clear(&p_temp);
}

mpz_t *mm_ell_fact(gmp_randstate_t state, const mpz_t n, unsigned long b1, unsigned long b2, unsigned long max_iter, unsigned long *iter, int *fase_found)
{
	m_ellp *p, *r;
	mpz_t *g, *g_r;
	mpz_temp temp;
	m_ellp_temp p_temp;
	m_ellp_rep rep;
	mpz_rep beta;

	mm_fact_param param;
	mpz_t *fact = NULL;
	const int FACT_REP_SIZE = 950;
	const unsigned long n_size = mpz_size(param.mdata.n) * mp_bits_per_limb;
	int vdiff_size = get_vdiff_size(b2);

	mpz_temp_init2(&temp, n_temp_mfact, n_size);
	m_ellp_temp_init2(&p_temp, n_p_temp_mfact, n_size);
	m_ellp_rep_init2(&rep, FACT_REP_SIZE, n_size);
	mpz_rep_init2(&beta, FACT_REP_SIZE, n_size);
	mpz_init2(param.k, bigk_size_bits(b1));
	mpz_init2(param.e_C2, n_size);

	
	mpz_temp_get(g, &temp);
	mpz_temp_get(g_r, &temp);
	m_ellp_temp_get(p, &p_temp);
	m_ellp_temp_get(r, &p_temp);
	*fase_found = -1;
	

	if (vdiff_size == -1) {
		perror("b2 to big\n");
		return NULL;
	}

	param.vdiff = malloc(vdiff_size);
	if (param.vdiff == NULL) {
		perror("error in mm_ell_fact malloc()");
		return NULL;
	}

	fact = malloc(sizeof(mpz_t));
	if (fact == NULL) {
		perror("error in m_ell_fact malloc()");
		free(param.vdiff);
		return NULL;
	}

	create_bigk(param.k, b1, &temp);
	get_prime_diff(b1, 1, b2, param.vdiff, &temp);
	mpz_realloc2(param.k, mpz_size(param.k) * mp_bits_per_limb);
	param.b1 = b1;
	param.b2 = b2;
	param.max_iter = max_iter;

	if (mform_data_init(&(param.mdata), n, &temp)){	// FOUND IN INVERTION OF R
		mpz_set(*fact, param.mdata.R2);
		*fase_found = 0;	//FASE 0
		goto found;
	}

	mpz_set_ui(*g, 1);
	to_mform(*g_r, *g, &(param.mdata), &temp);

	for (*iter = 0; *iter < param.max_iter; (*iter)++) {
		if (m_ell_setrand2(param.mdata.n, param.e_C2, p, state, &temp))	//TODO invertion can be avoited
		{
			if (find_div_by_gcd(*g, p->X, param.mdata.n)) {
				mpz_set(*fact, *g);
				*fase_found = 0;
				break;
			}
		} else {
			m_ellp_to_mform(p, p, &(param.mdata), &temp);
			to_mform(param.e_C2, param.e_C2, &(param.mdata), &temp);

			mm_ell_mul(param.k, param.e_C2, r, p, &(param.mdata), &p_temp, &temp);	//FASE1
			pthread_testcancel();

			m_ellp_from_mform(p, r, &(param.mdata), &temp);	// p = r from_mfrom
			if (find_div_by_gcd(*g, p->Z, param.mdata.n))	//check on p
			{
				mpz_set(*fact, *g);
				*fase_found = 1;
				break;
			}
			mpz_set(*g, *g_r);
			mm_ell_diff(rep.p, beta.v, rep.lenght, param.e_C2, r, &(param.mdata), &temp);
			mm_ell_fase2(*g, param.b1, param.b2, param.e_C2, r, rep.p, beta.v, rep.lenght, param.vdiff, &(param.mdata), &p_temp, &temp);
			pthread_testcancel();

			from_mform(*g, *g, &(param.mdata), &temp);
			if (find_div_by_gcd(*g, *g, param.mdata.n)) {
				mpz_set(*fact, *g);
				*fase_found = 2;
				break;
			}
		}
	}

	if (*fase_found == -1) {
		mpz_clear(*fact);
		fact = NULL;
	}

found:
	mpz_temp_free(&temp, 2);
	m_ellp_temp_free(&p_temp, 2);
	mpz_clears(param.k, NULL);
	free(param.vdiff);
	mform_data_clear(&(param.mdata));

	return fact;
}

