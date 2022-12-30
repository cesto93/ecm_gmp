#include "mpn_l.h"

static inline void mred_l(mp_limb_t *res_l, mp_limb_t *t_l, const mp_limb_t *n_l, mp_limb_t n_first, const size_t R_s);

int mform_data_init(mform_data *m_data, const mpz_t n, mpz_temp *temp)	//return 1 if can't invert R mod n
{
	mpz_t *R;
	mpz_temp_get(R, temp);
	mpz_inits(m_data->n, m_data->R2, NULL);
	
	mpz_set(m_data->n, n);
	m_data->n_l = mpz_limbs_read(m_data->n);
	m_data->n_s = mpz_size(n); // k = limbs(n)
	
	mpz_set_ui(*R, 1);
	mpz_mul_2exp(*R, *R, (m_data->n_s) * mp_bits_per_limb);// R = 2^k
	
	if (mpz_invert(m_data->R2, n, *R) == 0) //tr2 = 1/n mod R
	{
		mpz_set(m_data->R2, *R);
		mpz_temp_free(temp, 1);
		return 1;
	}
	
	mpz_neg(m_data->R2, m_data->R2); //tr2 = -1/n mod R1
	abs_modn(m_data->R2, *R);
	m_data->n_inv = mpz_getlimbn(m_data->R2, 0);
	DEBUG_STAMP("n_inv %lu\n", m_data->n_inv);
	
	mpz_mul(m_data->R2, *R, *R);
	mpz_mod(m_data->R2, m_data->R2, n); // R2 = R^2 mod n
	
	mpz_temp_free(temp, 1);
	return 0;
}

void msqr_l(mp_limb_t *rop_l, const mp_limb_t *op1_l, const size_t op1_s, mp_limb_t *t_l, const mp_limb_t *n_l, mp_limb_t n_first, 
			const size_t n_s)
{
	assert(op1_s <= n_s);
	mpn_sqr(t_l, op1_l, op1_s);
	
	if (op1_s < n_s)
		mpn_zero(t_l + (2 * op1_s), (2 * n_s) - (2 * op1_s));	
		
	mred_l(rop_l, t_l, n_l, n_first, n_s);
}

void mmul_l(mp_limb_t *rop_l, const mp_limb_t *op1_l, size_t op1_s, const mp_limb_t *op2_l, size_t op2_s, 
					mp_limb_t *t_l, const mp_limb_t *n_l, mp_limb_t n_first, const size_t n_s)
{
	//FIX_ME crash if tested on fase2 with prime n becouse we do mul g * 0
	if (op1_s < op2_s)
	{
		SWAP_LIMB_PTR(op1_l, op2_l);
		SWAP_SIZE_T(op1_s, op2_s);
	}
	assert(op1_s <= n_s);
	mpn_mul(t_l, op1_l, op1_s, op2_l, op2_s);
		
	if (op2_s != n_s)	//op2 < n_s => op1 + op2 < 2 * n_s
		mpn_zero(t_l + op1_s + op2_s, (2 * n_s) - op1_s - op2_s);		
		
	mred_l(rop_l, t_l, n_l, n_first, n_s);
}

void add_modR_l(mp_limb_t *rop_l, const mp_limb_t *op1_l, size_t op1_s, const mp_limb_t *op2_l, size_t op2_s, const mp_limb_t *n_l, 
				const size_t n_s)
{
	mp_limb_t cy;
	if (op1_s < op2_s) 							//mpn_add want 1st addend size >= 2nd addend size
	{
		SWAP_LIMB_PTR(op1_l, op2_l);
		SWAP_SIZE_T(op1_s, op2_s);
	}
	assert(op1_s <= n_s);
	cy = mpn_add(rop_l, op1_l, op1_s, op2_l, op2_s);
	
	if (op1_s != n_s)
	{
		mpn_zero(rop_l + op1_s, n_s - op1_s);
		if (cy)
		{
			mpn_add_1(rop_l + op1_s, rop_l + op1_s, n_s - op1_s, cy);
			cy = 0;
			if ((op1_s + 1 == n_s) && (mpn_cmp(rop_l, n_l, n_s) >= 0))	//we reduce to less than n 		
				mpn_sub_n(rop_l, rop_l, n_l, n_s);							
		}
	}
	else
		if ((cy) || (mpn_cmp(rop_l, n_l, n_s) >= 0))	//case (carry with n_s) => a + b = cy * R + rop => we reduce to less than R number
			cy -= mpn_sub_n(rop_l, rop_l, n_l, n_s);	
								
	assert(cy == 0);
	assert(mpn_cmp(rop_l, n_l, n_s) < 0);		
}

void sub_modR_l(mp_limb_t *rop_l, const mp_limb_t *op1_l, size_t op1_s, const mp_limb_t *op2_l, size_t op2_s, const mp_limb_t *n_l, 
				const size_t n_s)
{
	mp_limb_t cy;
	//gmp_printf("op1 = %Nd\nop2 = %Nd\nn = %Nd\n", op1_l, op1_s, op2_l, op2_s, n_l, n_s);
	
	if ((op1_s < op2_s) || ((op1_s == op2_s) && (mpn_cmp(op1_l, op2_l, op1_s) < 0))) // op1_l < op2_l //MOD_OP
	{
		assert(rop_l != op2_l); // WARNING rop_l used as temp and can't be 2nd operator
		assert(op2_s <= n_s);

		cy = mpn_add(rop_l, n_l, n_s, op1_l, op1_s);		
		cy -= mpn_sub(rop_l, rop_l, n_s, op2_l, op2_s); // res = op1 + n - op2			
	}
	else  // op1_l >= op2_l => mpn_sub positive / no borrow
	{
		assert(op1_s <= n_s);
		cy = mpn_sub(rop_l, op1_l, op1_s, op2_l, op2_s);	// res = op1 - op2
				
		if (op1_s != n_s)
			mpn_zero(rop_l + op1_s, n_s - op1_s);
	}
	assert(cy == 0);		
}

void msqr_l_n(mp_limb_t *rop_l, const mp_limb_t *op1_l, mp_limb_t *t_l, const mp_limb_t *n_l, mp_limb_t n_first, const size_t n_s)
{
	mpn_sqr(t_l, op1_l, n_s);
	mred_l(rop_l, t_l, n_l, n_first, n_s);
}

void mmul_l_n(mp_limb_t *rop_l, const mp_limb_t *op1_l, const mp_limb_t *op2_l, mp_limb_t *t_l, 
				const mp_limb_t *n_l, mp_limb_t n_first, const size_t n_s)
{	
	mpn_mul_n(t_l, op1_l ,op2_l, n_s);
	mred_l(rop_l, t_l, n_l, n_first, n_s);
}

void add_modR_l_n(mp_limb_t *rop_l, const mp_limb_t *op1_l, const mp_limb_t *op2_l, const mp_limb_t *n_l, const size_t n_s)
{
	mp_limb_t cy;
	cy = mpn_add_n(rop_l, op1_l, op2_l, n_s);
		
	if ((cy) || (mpn_cmp(rop_l, n_l, n_s) >= 0))	//case (carry with n_s) => a + b = cy * R + rop => we reduce to less than n
		cy -= mpn_sub_n(rop_l, rop_l, n_l, n_s);
	
	assert(mpn_cmp(rop_l, n_l, n_s) < 0);
	assert(cy == 0);		
}

void sub_modR_l_n(mp_limb_t *rop_l, const mp_limb_t *op1_l, const mp_limb_t *op2_l, const mp_limb_t *n_l, const size_t n_s)
{
	mp_limb_t bw;
	bw = mpn_sub_n(rop_l, op1_l, op2_l, n_s); // res = op1 + n - op2
		
	if (bw) //mod => this is used if op2_s is bigger than n => should be because of mmul
		bw -= mpn_add_n(rop_l, rop_l, n_l, n_s);
	assert(bw == 0);				
}

static inline void mred_l(mp_limb_t *res_l, mp_limb_t *t_l, const mp_limb_t *n_l, mp_limb_t n_first, const size_t n_s)
{
	mp_limb_t d;
	for (size_t i = n_s; i != 0; i--)
	{
		d = mpn_addmul_1(t_l, n_l, n_s, (n_first * t_l[0])); //t = t + N * (n_f * op[i] mod b)
		assert(t_l[0] == 0);
		t_l[0] = d;
		t_l++;
	}
	d = mpn_add_n(res_l, t_l, t_l - n_s, n_s); //no overlap
	if ((d != 0) || (mpn_cmp(res_l, n_l, n_s) >= 0)) 
		mpn_sub_n(res_l, res_l, n_l, n_s);

	assert(mpn_cmp(res_l, n_l, n_s) < 0);	
}

void mred_destroy(mpz_t res, mpz_t op, const mp_limb_t *n_l, const mp_limb_t n_first, const size_t n_s) //op != res destroy op
{
	mp_limb_t *op_l, *res_l; 
	const size_t op_s = mpz_size(op);
	
	assert(mpz_sgn(op) != -1);
	assert(mpz_size(op) <= (n_s * 2));
	
	op_l = mpz_limbs_modify(op, 2 * n_s * mp_bits_per_limb);
	res_l = mpz_limbs_write(res, n_s * mp_bits_per_limb);
	
	if (op_s < (2 * n_s))
		mpn_zero(op_l + op_s, (2 * n_s) - op_s);
	
	mred_l(res_l, op_l, n_l, n_first, n_s);
	mpz_limbs_finish(res, n_s);
}

void mred(mpz_t res, const mpz_t op, const mp_limb_t *n_l, const mp_limb_t n_first, const size_t n_s, mpz_temp *temp)
{
	mpz_t *t;
	mp_limb_t *t_l, *res_l; 
	const size_t op_s = mpz_size(op);
	
	assert(mpz_sgn(op) != -1);
	assert(op_s <= (n_s * 2));
	
	mpz_temp_get(t, temp);
	mpz_set(*t, op);
	t_l = mpz_limbs_modify(*t, 2 * n_s * mp_bits_per_limb);
	res_l = mpz_limbs_write(res, n_s * mp_bits_per_limb);
	
	if (op_s < (2 * n_s))
		mpn_zero(t_l + op_s, (2 * n_s) - op_s);		
	
	mred_l(res_l, t_l, n_l, n_first, n_s);
	
	mpz_limbs_finish(res, n_s);
	mpz_temp_free(temp, 1);	
}

