#ifndef MPN_L_H_		/* Include guard */
#define MPN_L_H_

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <errno.h>
#include <gmp.h>

#include "base_lib.h"
#include "matbase_lib.h"
#include "m_ellcurv_struct.h"

typedef struct mpl {
	mp_limb_t *l;
	size_t s;
} mpl;

typedef struct m_ellp_l {
	mp_limb_t *X;
	mp_limb_t *Z;
} m_ellp_l;

typedef struct m_ellp_lro {
	const mp_limb_t *X;
	size_t X_S;
	const mp_limb_t *Z;
	size_t Z_S;
} m_ellp_lro;

#define to_mpl_read(rop, op)			\
do {									\
	(rop)->l = mpz_limbs_read(*(op));	\
	(rop)->s = mpz_size(*(op));			\
} while (0)

#define to_mpl_write(rop, op, new_size)				\
do {												\
	(rop)->l = mpz_limbs_write(*(op), new_size);	\
	(rop)->s = mpz_size(*(op));						\
} while (0)

#define to_mpl_modify(rop, op, new_size)			\
do {												\
	(rop)->l = mpz_limbs_modify(*(op), new_size);	\
	(rop)->s = mpz_size(*(op));						\
} while (0)

typedef struct mform_data {
	mpz_t n;
	const mp_limb_t *n_l;
	mp_limb_t n_inv;	// -1/n mod limb
	mpz_t R2;		// R^2 mod n
	size_t n_s;
} mform_data;

int mform_data_init(mform_data * m_data, const mpz_t n, mpz_temp * temp);	//return 1 if can't invert R mod n
#define mform_data_clear(m_data) mpz_clears((m_data)->n, (m_data)->R2, NULL)

void mred(mpz_t res, const mpz_t op, const mp_limb_t * n_l,
	  const mp_limb_t n_inv, const size_t n_s, mpz_temp * temp);
void mred_destroy(mpz_t res, mpz_t op, const mp_limb_t * n_l, const mp_limb_t n_inv, const size_t n_s);	//op != res destroy op

void msqr_l(mp_limb_t * rop_l, const mp_limb_t * op1_l, const size_t op1_s,
	    mp_limb_t * t_l, const mp_limb_t * n_l, mp_limb_t n_first,
	    const size_t n_s);
		
void mmul_l(mp_limb_t * rop_l, const mp_limb_t * op1_l, size_t op1_s,
	    const mp_limb_t * op2_l, size_t op2_s, mp_limb_t * t_l,
	    const mp_limb_t * n_l, mp_limb_t n_inv, const size_t n_s);

void add_modR_l(mp_limb_t * rop_l, const mp_limb_t * op1_l, size_t op1_s,
		const mp_limb_t * op2_l, size_t op2_s, const mp_limb_t * n_l,
		const size_t n_s);

void sub_modR_l(mp_limb_t * rop_l, const mp_limb_t * op1_l, size_t op1_s,
		const mp_limb_t * op2_l, size_t op2_s, const mp_limb_t * n_l,
		const size_t n_s);

void msqr_l_n(mp_limb_t * rop_l, const mp_limb_t * op1_l, mp_limb_t * t_l,
	      const mp_limb_t * n_l, mp_limb_t n_first, const size_t n_s);
void mmul_l_n(mp_limb_t * rop_l, const mp_limb_t * op1_l,
	      const mp_limb_t * op2_l, mp_limb_t * t_l, const mp_limb_t * n_l,
	      mp_limb_t n_first, const size_t n_s);
void add_modR_l_n(mp_limb_t * rop_l, const mp_limb_t * op1_l,
		  const mp_limb_t * op2_l, const mp_limb_t * n_l,
		  const size_t n_s);
void sub_modR_l_n(mp_limb_t * rop_l, const mp_limb_t * op1_l,
		  const mp_limb_t * op2_l, const mp_limb_t * n_l,
		  const size_t n_s);

#define from_mform(rop, op, m_data, temp) mred(rop, op, (m_data)->n_l, (m_data)->n_inv, (m_data)->n_s, temp)
#define from_mform_d(rop, op, m_data) mred_destroy(rop, op, (m_data)->n_l, (m_data)->n_inv, (m_data)->n_s)
#define to_mform(rop, op, m_data, temp) mmul_t((rop), (op), (m_data)->R2, (m_data)->n_l, (m_data)->n_inv, (m_data)->n_s, (temp))

#define add_modR(rop, op1, op2, n_l, n_s)										\
do {																			\
	mp_limb_t *rop_l;															\
	const mp_limb_t *op1_l, *op2_l;												\
																				\
	rop_l = mpz_limbs_write(rop, n_s * mp_bits_per_limb);						\
	op1_l = mpz_limbs_read(op1);												\
	op2_l = mpz_limbs_read(op2);												\
	add_modR_l(rop_l, op1_l, mpz_size(op1), op2_l, mpz_size(op2), n_l, n_s);	\
	mpz_limbs_finish(rop, n_s);													\
																				\
} while (0)

#define sub_modR(rop, op1, op2, n_l , n_s)										\
do {																			\
	mp_limb_t *rop_l;															\
	const mp_limb_t *op1_l, *op2_l;												\
																				\
	rop_l = mpz_limbs_write(rop, n_s * mp_bits_per_limb);						\
	op1_l = mpz_limbs_read(op1);												\
	op2_l = mpz_limbs_read(op2);												\
	sub_modR_l(rop_l, op1_l, mpz_size(op1), op2_l, mpz_size(op2), n_l, n_s);	\
	mpz_limbs_finish(rop, n_s);													\
																				\
} while (0)

#define mmul(rop, op1, op2, t_l, n_l, n_inv, n_s)										\
do {																					\
	mp_limb_t *rop_l; 																	\
	const mp_limb_t *op1_l, *op2_l;														\
																						\
	rop_l = mpz_limbs_write(rop, n_s * mp_bits_per_limb);								\
	op1_l = mpz_limbs_read(op1);														\
	op2_l = mpz_limbs_read(op2);														\
																						\
	mmul_l(rop_l, op1_l, mpz_size(op1), op2_l, mpz_size(op2), t_l, n_l, n_inv, n_s);	\
	mpz_limbs_finish(rop, n_s);															\
																						\
} while (0)

#define mmul_t(rop, op1, op2, n_l, n_inv, n_s, temp)									\
do {																					\
	mpz_t *t;																			\
	mp_limb_t *t_l;																		\
	mpz_temp_get(t, temp);																\
	t_l = mpz_limbs_write(*t, 2 * n_s * mp_bits_per_limb);								\
	mmul(rop, op1, op2, t_l, n_l, n_inv, n_s);											\
	mpz_temp_free(temp, 1);																\
} while (0)

#define msqr(rop, op1, t_l, n_l, n_inv, n_s)											\
do {																					\
	mp_limb_t *rop_l; 																	\
	const mp_limb_t *op1_l;																\
																						\
	rop_l = mpz_limbs_write(rop, n_s * mp_bits_per_limb);								\
	op1_l = mpz_limbs_read(op1);														\
																						\
	msqr_l(rop_l, op1_l, mpz_size(op1), t_l, n_l, n_inv, n_s);							\
	mpz_limbs_finish(rop, n_s);															\
} while (0)

#define MPZ_MODIFY_NSIZE(rop, op, n_s)							\
do {															\
	rop = mpz_limbs_modify(op, n_s * mp_bits_per_limb);			\
	size_t size = mpz_size(op);									\
	if (size <= n_s)											\
		mpn_zero(rop + size, n_s - size);						\
} while (0)

#define MPZ_READ_WITH_SIZE(rop, rop_s, op)						\
do {															\
	rop = mpz_limbs_read((op));									\
	rop_s = mpz_size(op);										\
} while (0)

#define NORMALIZE(l_ptr, size)			\
do {									\
	while ((size) > 0)					\
	{									\
		if ((l_ptr) [(size) - 1] != 0)	\
			break;						\
		(size) --;						\
	}									\
} while (0)

#define SWAP_SIZE_T(x, y)		\
do { 							\
	x = x ^ y;					\
	y = x ^ y;					\
	x = x ^ y;					\
} while (0)

#define SWAP_LIMB_PTR(x, y)		\
do { 							\
	const mp_limb_t *t = (x);	\
	(x) = (y);					\
	(y) = t;					\
} while (0)

#define SWAP_LIMB_PTR_W(x, y)	\
do { 							\
	mp_limb_t *t = (x);			\
	(x) = (y);					\
	(y) = t;					\
} while (0)

#endif //MPN_L_H_
