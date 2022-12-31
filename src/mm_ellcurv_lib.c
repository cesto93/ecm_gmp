#include "mm_ellcurv_lib.h"

static inline void mm_ell_addh_l(mp_limb_t * rx_l, mp_limb_t * rz_l,
				 const mp_limb_t * px_l, size_t px_s,
				 const mp_limb_t * pz_l, size_t pz_s,
				 const mp_limb_t * qx_l, size_t qx_s,
				 const mp_limb_t * qz_l, size_t qz_s,
				 const mp_limb_t * dx_l, size_t dx_s,
				 const mp_limb_t * dz_l, size_t dz_s,
				 const mform_data * mdata, mp_limb_t * t1_l,
				 mp_limb_t * t2_l, mp_limb_t * t3_l,
				 mp_limb_t * tmul_l);

static inline void mm_ell_duph_l(const mp_limb_t * e_C2_l, size_t e_C2_s,
				 mp_limb_t * rx_l, mp_limb_t * rz_l,
				 const mp_limb_t * px_l, size_t px_s,
				 const mp_limb_t * pz_l, size_t pz_s,
				 const mform_data * mdata, mp_limb_t * t1_l,
				 mp_limb_t * t2_l, mp_limb_t * t3_l,
				 mp_limb_t * tmul_l);

static inline void mm_ell_addh_l_n(mp_limb_t * rx_l, mp_limb_t * rz_l,
				   const mp_limb_t * px_l,
				   const mp_limb_t * pz_l,
				   const mp_limb_t * qx_l,
				   const mp_limb_t * qz_l,
				   const mp_limb_t * dx_l,
				   const mp_limb_t * dz_l,
				   const mform_data * mdata, mp_limb_t * t1_l,
				   mp_limb_t * t2_l, mp_limb_t * t3_l,
				   mp_limb_t * tmul_l);

static inline void mm_ell_duph_l_n(const mp_limb_t * e_C2_l, mp_limb_t * rx_l,
				   mp_limb_t * rz_l, const mp_limb_t * px_l,
				   const mp_limb_t * pz_l,
				   const mform_data * mdata, mp_limb_t * t1_l,
				   mp_limb_t * t2_l, mp_limb_t * t3_l,
				   mp_limb_t * tmul_l);

static void mm_ell_mul_ui(unsigned long k, const mpz_t e_C2, m_ellp * r,
			  const m_ellp * p, const mform_data * mdata,
			  m_ellp_temp * p_temp, mpz_temp * temp);

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
void mm_ell_mul(const mpz_t k, mpz_t e_C2, m_ellp * r, m_ellp * p,
		const mform_data * mdata, m_ellp_temp * p_temp, mpz_temp * temp)
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

	////////////////////////////////////////////////////////////////////////// LOW LEVEL
	mp_limb_t *rx_l, *rz_l, *ux_l, *uz_l, *tx_l, *tz_l;
	mp_limb_t *px_l, *pz_l, *e_C2_l;	//CONST

	rx_l = mpz_limbs_write((r)->X, (mdata->n_s) * mp_bits_per_limb);
	rz_l = mpz_limbs_write((r)->Z, (mdata->n_s) * mp_bits_per_limb);
	ux_l = mpz_limbs_write((u)->X, (mdata->n_s) * mp_bits_per_limb);
	uz_l = mpz_limbs_write((u)->Z, (mdata->n_s) * mp_bits_per_limb);
	tx_l = mpz_limbs_write((t)->X, (mdata->n_s) * mp_bits_per_limb);
	tz_l = mpz_limbs_write((t)->Z, (mdata->n_s) * mp_bits_per_limb);

	MPZ_MODIFY_NSIZE(px_l, p->X, mdata->n_s);
	MPZ_MODIFY_NSIZE(pz_l, p->Z, mdata->n_s);
	MPZ_MODIFY_NSIZE(e_C2_l, e_C2, mdata->n_s);

	mpn_copyi(ux_l, px_l, mdata->n_s);	//SET_ELLP(u,p)
	mpn_copyi(uz_l, pz_l, mdata->n_s);

	mm_ell_duph_l_n(e_C2_l, tx_l, tz_l, px_l, pz_l, mdata, t1_l, t2_l, t3_l,
			tmul_l);
	///////////////////////////////////////////////////////////////////////////////////

	unsigned long int high;
	high = mpz_sizeinbase(k, 2) - 1;

	if (mpz_cmp_ui(k, 2) > 0) {
		for (unsigned long i = high - 1; i > 0; i--)
			if (mpz_tstbit(k, i)) {
				mm_ell_addh_l_n(ux_l, uz_l, tx_l, tz_l, ux_l,
						uz_l, px_l, pz_l, mdata, t1_l,
						t2_l, t3_l, tmul_l);
				mm_ell_duph_l_n(e_C2_l, tx_l, tz_l, tx_l, tz_l,
						mdata, t1_l, t2_l, t3_l,
						tmul_l);
			} else {
				mm_ell_addh_l_n(tx_l, tz_l, ux_l, uz_l, tx_l,
						tz_l, px_l, pz_l, mdata, t1_l,
						t2_l, t3_l, tmul_l);
				mm_ell_duph_l_n(e_C2_l, ux_l, uz_l, ux_l, uz_l,
						mdata, t1_l, t2_l, t3_l,
						tmul_l);
			}
		if (mpz_tstbit(k, 0))
			mm_ell_addh_l_n(rx_l, rz_l, ux_l, uz_l, tx_l, tz_l,
					px_l, pz_l, mdata, t1_l, t2_l, t3_l,
					tmul_l);
		else
			mm_ell_duph_l_n(e_C2_l, rx_l, rz_l, ux_l, uz_l, mdata,
					t1_l, t2_l, t3_l, tmul_l);

		mpz_limbs_finish(r->X, mdata->n_s);
		mpz_limbs_finish(r->Z, mdata->n_s);
	} else {
		gmp_printf("k %Zd\n", k);
		if (mpz_sgn(k) != 1)
			error_msg("error in ell_mul k = 0");
		if (mpz_cmp_ui(k, 2) == 0)
			mm_ell_duph_t(e_C2, r, p, mdata, t1_l, t2_l, t3_l,
				      tmul_l);
		else if (mpz_cmp_ui(k, 1) == 0)
			m_ellp_set(r, p);
	}

	mpz_temp_free(temp, 4);
	m_ellp_temp_free(p_temp, 2);
}

//USED TWICE FOR ITER
static void mm_ell_mul_ui(unsigned long k, const mpz_t e_C2, m_ellp * r,
			  const m_ellp * p, const mform_data * mdata,
			  m_ellp_temp * p_temp, mpz_temp * temp)
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
				mm_ell_addh_t(u, t, u, p, mdata, t1_l, t2_l,
					      t3_l, tmul_l);
				mm_ell_duph_t(e_C2, t, t, mdata, t1_l, t2_l,
					      t3_l, tmul_l);
			} else {
				mm_ell_addh_t(t, u, t, p, mdata, t1_l, t2_l,
					      t3_l, tmul_l);
				mm_ell_duph_t(e_C2, u, u, mdata, t1_l, t2_l,
					      t3_l, tmul_l);
			}
		if (k & 1)	//k[0]
			mm_ell_addh_t(r, u, t, p, mdata, t1_l, t2_l, t3_l,
				      tmul_l);
		else
			mm_ell_duph_t(e_C2, r, u, mdata, t1_l, t2_l, t3_l,
				      tmul_l);
	} else {
		printf("k %lu\n", k);
		if (k <= 0)
			error_msg("error in ell_mul k = 0");
		if (k == 2)
			mm_ell_duph_t(e_C2, r, p, mdata, t1_l, t2_l, t3_l,
				      tmul_l);
		else if (k == 1)
			m_ellp_set(r, p);
	}

	mpz_temp_free(temp, 4);
	mpz_temp_free(p_temp, 2);
}

//USED ONCE FOR ITER / VERY LIGHT
void mm_ell_diff(m_ellp * rep, mpz_t * beta, const unsigned long d, mpz_t e_C2,
		 const m_ellp * p, const mform_data * mdata, mpz_temp * temp)
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
	mm_ell_duph_t(e_C2, &(rep[1]), &(rep[0]), mdata, t1_l, t2_l, t3_l,
		      tmul_l);
	mmul(beta[1], rep[1].X, rep[1].Z, tmul_l, mdata->n_l, mdata->n_inv, mdata->n_s);	// beta[1] = X1 * Z1

	for (unsigned long i = 2; i < d; i++) {
		mm_ell_addh_t(&(rep[i]), &(rep[i - 1]), &(rep[0]),
			      &(rep[i - 2]), mdata, t1_l, t2_l, t3_l, tmul_l);
		mmul(beta[i], rep[i].X, rep[i].Z, tmul_l, mdata->n_l, mdata->n_inv, mdata->n_s);	//beta[i] = rep[i].X * rep[i].Z
	}

	mpz_temp_free(temp, 4);
}

void mm_ell_fase2(mpz_t g, const unsigned long b1, const unsigned long b2,
		  const mpz_t e_C2, const m_ellp * p, m_ellp * S, mpz_t * beta,
		  const unsigned long d, const unsigned char vdiff[],
		  const mform_data * mdata, m_ellp_temp * p_temp,
		  mpz_temp * temp)
{
	assert(mpz_temp_space(temp) >= n_temp_mfase2);
	assert(m_ellp_temp_space(p_temp) >= n_p_temp_mfase2);

	const unsigned long int max_diff = 2 * d;	//max_diff = 2 * d
	unsigned long int diff = 0, b = b1 - 1;
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

	/////////////////////////////////////////////// LOW LEVEL
	mp_limb_t *t1_l, *t2_l, *t3_l, *tmul_l, *alfa_l, *g_l, *Rx_l, *Rz_l,
	    *Tx_l, *Tz_l;
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
	//////////////////////////////////////////////////////////////////////////////////

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
				 const mform_data * mdata, mp_limb_t * t1_l,
				 mp_limb_t * t2_l, mp_limb_t * t3_l,
				 mp_limb_t * tmul_l)
{
	sub_modR_l(t1_l, px_l, px_s, pz_l, pz_s, mdata->n_l, mdata->n_s);	//T3 = (x1-z1)(x2+z2)                                                                                                                                                                   
	add_modR_l(t2_l, qx_l, qx_s, qz_l, qz_s, mdata->n_l, mdata->n_s);	//              T1              T2
	mmul_l_n(t3_l, t1_l, t2_l, tmul_l, mdata->n_l, mdata->n_inv,
		 mdata->n_s);

	add_modR_l(t1_l, px_l, px_s, pz_l, pz_s, mdata->n_l, mdata->n_s);	//T1 = (x1+z1)(x2-z2)
	//              T1              T2
	sub_modR_l(t2_l, qx_l, qx_s, qz_l, qz_s, mdata->n_l, mdata->n_s);
	mmul_l_n(t1_l, t1_l, t2_l, tmul_l, mdata->n_l, mdata->n_inv,
		 mdata->n_s);

	add_modR_l_n(t2_l, t3_l, t1_l, mdata->n_l, mdata->n_s);	//T2 = (x1-z1)(x2+z2) + (x1+z1)(x2-z2)
	//                      T3                              T1
	sub_modR_l_n(t3_l, t3_l, t1_l, mdata->n_l, mdata->n_s);	//T3 = (x1-z1)(x2+z2) - (x1+z1)(x2-z2)
	//                      T3                              T1
	msqr_l_n(t1_l, t2_l, tmul_l, mdata->n_l, mdata->n_inv, mdata->n_s);	//T1 = [(x1-z1)(x2+z2) + (x1+z1)(x2-z2)]^2
	//                                              T2
	msqr_l_n(t2_l, t3_l, tmul_l, mdata->n_l, mdata->n_inv, mdata->n_s);	//T2 = [(x1-z1)(x2+z2) - (x1+z1)(x2-z2)]^2
	//                                      T3
	mmul_l(t1_l, t1_l, mdata->n_s, dz_l, dz_s, tmul_l, mdata->n_l, mdata->n_inv, mdata->n_s);	//T1 = [((x1-z1)(x2+z2) + (x1+z1)(x2-z2))^2] * Z- = RX

	mmul_l(rz_l, t2_l, mdata->n_s, dx_l, dx_s, tmul_l, mdata->n_l, mdata->n_inv, mdata->n_s);	//RZ = [((x1-z1)(x2+z2) - (x1+z1)(x2-z2))^2] * X-

	mpn_copyi(rx_l, t1_l, mdata->n_s);	//RX = t1 
}

static inline void mm_ell_duph_l(const mp_limb_t * e_C2_l, size_t e_C2_s,
				 mp_limb_t * rx_l, mp_limb_t * rz_l,
				 const mp_limb_t * px_l, size_t px_s,
				 const mp_limb_t * pz_l, size_t pz_s,
				 const mform_data * mdata, mp_limb_t * t1_l,
				 mp_limb_t * t2_l, mp_limb_t * t3_l,
				 mp_limb_t * tmul_l)
{
	add_modR_l(t1_l, px_l, px_s, pz_l, pz_s, mdata->n_l, mdata->n_s);	// T2 = (x1+z1)^2                                                                                                                                                               
	msqr_l_n(t2_l, t1_l, tmul_l, mdata->n_l, mdata->n_inv, mdata->n_s);	//                      T1

	sub_modR_l(t1_l, px_l, px_s, pz_l, pz_s, mdata->n_l, mdata->n_s);	// T3 = (x1-z1)^2                                                                                                                                                       
	msqr_l_n(t3_l, t1_l, tmul_l, mdata->n_l, mdata->n_inv, mdata->n_s);	//                      T1

	mmul_l_n(rx_l, t2_l, t3_l, tmul_l, mdata->n_l, mdata->n_inv, mdata->n_s);	// RX = [(x1+z1)^2] * [(x1-z1)^2]

	sub_modR_l_n(t1_l, t2_l, t3_l, mdata->n_l, mdata->n_s);	//T1 = 4x1z1 = [(x1+z1)^2] - [(x1-z1)^2]
	//                              T2                              T3
	mmul_l(t2_l, e_C2_l, e_C2_s, t1_l, mdata->n_s, tmul_l, mdata->n_l, mdata->n_inv, mdata->n_s);	//T2 = [c+2/4] * [4x1z1]
	//              EC_2            T1
	add_modR_l_n(t3_l, t3_l, t2_l, mdata->n_l, mdata->n_s);	//T3 = [(x1-z1)^2] + [((c+2)/4) * 4x1z1]
	//                      T3                              T2
	mmul_l_n(rz_l, t3_l, t1_l, tmul_l, mdata->n_l, mdata->n_inv, mdata->n_s);	//RZ = 4x1z1 * [(x1-z1)^2 + ((c+2)/4) * 4x1z1]
}

static inline void mm_ell_addh_l_n(mp_limb_t * rx_l, mp_limb_t * rz_l,
				   const mp_limb_t * px_l,
				   const mp_limb_t * pz_l,
				   const mp_limb_t * qx_l,
				   const mp_limb_t * qz_l,
				   const mp_limb_t * dx_l,
				   const mp_limb_t * dz_l,
				   const mform_data * mdata, mp_limb_t * t1_l,
				   mp_limb_t * t2_l, mp_limb_t * t3_l,
				   mp_limb_t * tmul_l)
{
	sub_modR_l_n(t1_l, px_l, pz_l, mdata->n_l, mdata->n_s);	//T3 = (x1-z1)(x2+z2)                                                                                                                                                                   
	add_modR_l_n(t2_l, qx_l, qz_l, mdata->n_l, mdata->n_s);	//              T1              T2
	mmul_l_n(t3_l, t1_l, t2_l, tmul_l, mdata->n_l, mdata->n_inv,
		 mdata->n_s);

	add_modR_l_n(t1_l, px_l, pz_l, mdata->n_l, mdata->n_s);	//T1 = (x1+z1)(x2-z2)                                                                                                                                                           
	sub_modR_l_n(t2_l, qx_l, qz_l, mdata->n_l, mdata->n_s);	//              T1              T2
	mmul_l_n(t1_l, t1_l, t2_l, tmul_l, mdata->n_l, mdata->n_inv,
		 mdata->n_s);

	add_modR_l_n(t2_l, t3_l, t1_l, mdata->n_l, mdata->n_s);	//T2 = (x1-z1)(x2+z2) + (x1+z1)(x2-z2)
	//                      T3                              T1
	sub_modR_l_n(t3_l, t3_l, t1_l, mdata->n_l, mdata->n_s);	//T3 = (x1-z1)(x2+z2) - (x1+z1)(x2-z2)
	//                      T3                              T1
	msqr_l_n(t1_l, t2_l, tmul_l, mdata->n_l, mdata->n_inv, mdata->n_s);	//T1 = [(x1-z1)(x2+z2) + (x1+z1)(x2-z2)]^2
	//                                              T2
	msqr_l_n(t2_l, t3_l, tmul_l, mdata->n_l, mdata->n_inv, mdata->n_s);	//T2 = [(x1-z1)(x2+z2) - (x1+z1)(x2-z2)]^2
	//                                      T3
	mmul_l_n(t1_l, t1_l, dz_l, tmul_l, mdata->n_l, mdata->n_inv, mdata->n_s);	//T1 = [((x1-z1)(x2+z2) + (x1+z1)(x2-z2))^2] * Z- = RX

	mmul_l_n(rz_l, t2_l, dx_l, tmul_l, mdata->n_l, mdata->n_inv, mdata->n_s);	//RZ = [((x1-z1)(x2+z2) - (x1+z1)(x2-z2))^2] * X-

	mpn_copyi(rx_l, t1_l, mdata->n_s);	//RX = t1 
}

static inline void mm_ell_duph_l_n(const mp_limb_t * e_C2_l, mp_limb_t * rx_l,
				   mp_limb_t * rz_l, const mp_limb_t * px_l,
				   const mp_limb_t * pz_l,
				   const mform_data * mdata, mp_limb_t * t1_l,
				   mp_limb_t * t2_l, mp_limb_t * t3_l,
				   mp_limb_t * tmul_l)
{
	add_modR_l_n(t1_l, px_l, pz_l, mdata->n_l, mdata->n_s);	// T2 = (x1+z1)^2                                                                                                                                                               
	msqr_l_n(t2_l, t1_l, tmul_l, mdata->n_l, mdata->n_inv, mdata->n_s);	//                      T1

	sub_modR_l_n(t1_l, px_l, pz_l, mdata->n_l, mdata->n_s);	// T3 = (x1-z1)^2                                                                                                                                                       
	msqr_l_n(t3_l, t1_l, tmul_l, mdata->n_l, mdata->n_inv, mdata->n_s);	//                      T1

	mmul_l_n(rx_l, t2_l, t3_l, tmul_l, mdata->n_l, mdata->n_inv, mdata->n_s);	// RX = [(x1+z1)^2] * [(x1-z1)^2]

	sub_modR_l_n(t1_l, t2_l, t3_l, mdata->n_l, mdata->n_s);	//T1 = 4x1z1 = [(x1+z1)^2] - [(x1-z1)^2]
	//                              T2                              T3
	mmul_l_n(t2_l, e_C2_l, t1_l, tmul_l, mdata->n_l, mdata->n_inv, mdata->n_s);	//T2 = [c+2/4] * [4x1z1]
	//              EC_2            T1
	add_modR_l_n(t3_l, t3_l, t2_l, mdata->n_l, mdata->n_s);	//T3 = [(x1-z1)^2] + [((c+2)/4) * 4x1z1]
	//                      T3                              T2
	mmul_l_n(rz_l, t3_l, t1_l, tmul_l, mdata->n_l, mdata->n_inv, mdata->n_s);	//RZ = 4x1z1 * [(x1-z1)^2 + ((c+2)/4) * 4x1z1]
}

//WRAPPER FUNCTION FOR TESTING / DON'T USE TEMP

void mm_ell_mul_t(const mpz_t k, const mpz_t n, mpz_t e_C2, m_ellp * r,
		  m_ellp * p)
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

//when used first trasform to mform
/*static inline void mm_ell_addh(m_ellp *r, const m_ellp *p, const m_ellp *q, const m_ellp *diff, const mform_data *m_data, mpz_temp *temp) 
{
	mpz_t *t1, *t2, *t3, *tmul;
	mp_limb_t *tmul_l;
	
	mpz_temp_get(t1, temp);
	mpz_temp_get(t2, temp);
	mpz_temp_get(t3, temp);
	mpz_temp_get(tmul, temp);
	
	tmul_l = mpz_limbs_write(*tmul, 2 * (m_data->n_s) * mp_bits_per_limb);
	
	mpz_sub(*t1, p->X, p->Z);		//sub
	abs_modn(*t1, m_data->n);
	//sub_modR(*t1, p->X, p->Z, m_data->n_l, m_data->n_s);
	
	mpz_add(*t2, q->X, q->Z);
	//add_modR(*t2, q->X, q->Z, m_data->n_s);
	
	//t3 = (x1-z1)(x2+z2)
	mmul(*t3, *t1, *t2, tmul_l, m_data->n_l, m_data->n_first, m_data->n_s);
	
	mpz_add(*t1, p->X, p->Z);
	//add_modR(*t1, p->X, p->Z, m_data->n_s);
	
	mpz_sub(*t2, q->X, q->Z); //sub
	abs_modn(*t2, m_data->n);
	//sub_modR(*t2, q->X, q->Z, m_data->n_l, m_data->n_s);
	
	//t1 = (x1+z1)(x2-z2)
	mmul(*t1, *t1, *t2, tmul_l, m_data->n_l, m_data->n_first, m_data->n_s);
	
	mpz_add(*t2, *t3, *t1); //t2 = (x1-z1)(x2+z2) + (x1+z1)(x2-z2)
	//add_modR(*t2, *t3, *t1, m_data->n_s);
	
	mpz_sub(*t3, *t3, *t1); //t3 = (x1-z1)(x2+z2) - (x1+z1)(x2-z2) no abs_modn needed
	//sub_modR(*t3, *t3, *t1, m_data->n_l, m_data->n_s);

	//t1 = t2^2 = [(x1-z1)(x2+z2) + (x1+z1)(x2-z2)]^2
	msqr(*t1, *t2, tmul_l, m_data->n_l, m_data->n_first, m_data->n_s);
	
	//t2 = t3^2 = [(x1-z1)(x2+z2) - (x1+z1)(x2-z2)]^2
	msqr(*t2, *t3, tmul_l, m_data->n_l, m_data->n_first, m_data->n_s); 
	 
	//t1 = [((x1-z1)(x2+z2) + (x1+z1)(x2-z2))^2] * z- = RX
	mmul(*t1, *t1, diff->Z, tmul_l, m_data->n_l, m_data->n_first, m_data->n_s);
	 
	//RZ = [((x1-z1)(x2+z2) - (x1+z1)(x2-z2))^2] * x-
	mmul(r->Z, *t2, diff->X, tmul_l, m_data->n_l, m_data->n_first, m_data->n_s);
	
	mpz_set(r->X, *t1);  //rx = t1

	mpz_temp_free(temp, 4);
}

//when used first transform to m_form
static inline void mm_ell_duph(const mpz_t e_C2, m_ellp *r, const m_ellp *p, const mform_data *m_data, mpz_temp *temp)
{
	mpz_t *t1, *t2, *t3, *tmul;
	mp_limb_t *tmul_l;
	mpz_temp_get(t1, temp);
	mpz_temp_get(t2, temp);
	mpz_temp_get(t3, temp);
	mpz_temp_get(tmul, temp);
	
	tmul_l = mpz_limbs_write(*tmul, 2 * (m_data->n_s) * mp_bits_per_limb);
	
	// t2 = (x1+z1)^2
	mpz_add(*t1, p->X, p->Z);
	
	msqr(*t2, *t1, tmul_l, m_data->n_l, m_data->n_first, m_data->n_s);
	
	// t3 = (x1-z1)^2
	mpz_sub(*t1, p->X, p->Z); //t1 no abs_modn needed square
	
	msqr(*t3, *t1, tmul_l, m_data->n_l, m_data->n_first, m_data->n_s);
	
	// rx = [(x1+z1)^2] * [(x1-z1)^2]
	mmul(r->X, *t2, *t3, tmul_l, m_data->n_l, m_data->n_first, m_data->n_s);
	
	mpz_sub(*t1, *t2, *t3); //t1 = 4x1z1 = [(x1+z1)^2] - [(x1-z1)^2]
	abs_modn(*t1, m_data->n);
	
	//t2 = [c+2/4] * [4x1z1]
	mmul(*t2, e_C2, *t1, tmul_l, m_data->n_l, m_data->n_first, m_data->n_s);
	
	mpz_add(*t3, *t3, *t2); //t3 = [(x1-z1)^2] + [((c+2)/4) * 4x1z1]
	
	//rz = 4x1z1 * [(x1-z1)^2 + ((c+2)/4) * 4x1z1]
	mmul(r->Z, *t3, *t1, tmul_l, m_data->n_l, m_data->n_first, m_data->n_s);
	
	mpz_temp_free(temp, 4);
}*/