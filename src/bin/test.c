#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>

#include "../base_lib.h"
#include "../m_ellcurv_struct.h"
#include "../m_ellcurv_fact.h"
#include "../m_ellcurv_lib.h"
#include "../mm_ellcurv_lib.h"

#define CHECK_ARGC(param_n) if (argc < (param_n + 1)) \
do { \
	perror("usage: ./test.o  mod \ntentative digits_start:inc:size b_start:inc:size\n"); \
	exit(EXIT_FAILURE); \
} while (0)

#define DEFAULT_COEFF 100	// b2 = b1 * coeff
#define MAX_ITER 2000ul

#define FACT_IS_CORRECT(p, q, f1, f2) ((mpz_equal(p, f1) && mpz_equal(q, f2)) || ((mpz_equal(q, f1) && mpz_equal(p, f2))))

static void test_fact_times(const char *tent_s, const char *digits_s, const char *b_s, const char *path_t, const char *path_f);
static void test_fact(FILE * log_t, FILE * log_f, unsigned long tent, unsigned long digits, unsigned long b1,
		      unsigned long b2, gmp_randstate_t state);
			  
void testpoint(const char *n);
void testdup(const char *n);
void testmul(const char *n, const char *k_s);
void test_diff(const char *n_s);
void test_prime_diff(const char *start_s, const char *end_s);

void test_mul_ui(const char *n_s, const char *k_s);
void test_mmul(const char *n_s, const char *k_s);
void test_mfase2(const char *n_s, const char *b1_s, const char *b2_s);

int main(int argc, char *argv[])
{
	CHECK_ARGC(1);
	if (strcmp(argv[1], "fact") == 0) {
		CHECK_ARGC(6);
		test_fact_times(argv[2], argv[3], argv[4], argv[5], argv[6]);
		exit(EXIT_SUCCESS);
	}
	if (strcmp(argv[1], "mmul") == 0) {
		CHECK_ARGC(3);
		test_mmul(argv[2], argv[3]);
		exit(EXIT_SUCCESS);
	}
	if (strcmp(argv[1], "mul") == 0) {
		CHECK_ARGC(3);
		testmul(argv[2], argv[3]);
		exit(EXIT_SUCCESS);
	}
	if (strcmp(argv[1], "mfase2") == 0) {
		CHECK_ARGC(4);
		test_mfase2(argv[2], argv[3], argv[4]);
		exit(EXIT_SUCCESS);
	}
	if (strcmp(argv[1], "mul_ui") == 0) {
		CHECK_ARGC(3);
		test_mul_ui(argv[2], argv[3]);
		exit(EXIT_SUCCESS);
	}
	if (strcmp(argv[1], "prime_diff") == 0) {
		CHECK_ARGC(3);
		test_prime_diff(argv[2], argv[3]);
		exit(EXIT_SUCCESS);
	}
}

static void test_fact_times(const char *tent_s, const char *digits_s, const char *b_s, const char *path_t, const char *path_f)
{
	unsigned long b1_inc, b2_inc, b_start, b1, b2, digits;
	unsigned int b_size;
	unsigned long tent;
	FILE *log_t, *log_f;
	gmp_randstate_t state;

	gmp_randinit_default(state);
	gmp_randseed_ui(state, time(NULL));
	tent = str2l(tent_s);
	digits = str2l(digits_s);

	get_arg(b_s, &b_start, &b1_inc, &b_size);
	b2_inc = b1_inc * DEFAULT_COEFF;

	log_t = fopen(path_t, "a");
	if (log_t == NULL) {
		perror("error in fopen");
		return;
	}

	log_f = fopen(path_f, "a");
	if (log_f == NULL) {
		fclose(log_t);
		perror("error in fopen");
		return;
	}

	b1 = b_start;
	b2 = b1 * DEFAULT_COEFF;
	for (unsigned int i = 0; i < b_size; i++) {
		test_fact(log_t, log_f, tent, digits, b1, b2, state);
		fputs("\n", log_f);
		fflush(log_t);
		fflush(log_f);

		b1 += b1_inc;
		b2 += b2_inc;
	}
	fclose(log_t);
	fclose(log_f);
}

void testpoint(const char *n)
{
	m_ellc e;
	m_ellp p;
	gmp_randstate_t state;

	gmp_randinit_default(state);
	m_ellp_init(&p);
	m_ellc_init(&e);

	gmp_randseed_ui(state, time(NULL));
	mpz_set_str(e.n, n, 10);
	if (m_ell_setrand2_t(e.n, e.C2, &p, state))
		error_msg("error in setrand at testpoint\n");
	gmp_printf("c_2 = %Zd\t X = %Zd\t Z = %Zd\n", e.C2, p.X, p.Z);
}

void testdup(const char *n)
{
	m_ellc e;
	m_ellp p, p2, p3, p4d, p4a;
	gmp_randstate_t state;

	gmp_randinit_default(state);
	m_ellp_init(&p);
	m_ellp_init(&p2);
	m_ellp_init(&p3);
	m_ellp_init(&p4a);
	m_ellp_init(&p4d);
	m_ellc_init(&e);

	gmp_randseed_ui(state, time(NULL));
	mpz_set_str(e.n, n, 10);
	if (m_ell_setrand2_t(e.n, e.C2, &p, state))
		error_msg("error in set_rand at testdup\n");
	gmp_printf("c_2 = %Zd\nX1 = %Zd\t Z1 = %Zd\n", e.C2, p.X, p.Z);

	m_ell_duph_t(&e, &p2, &p);
	gmp_printf("X2 = %Zd\t Z2 = %Zd\n", p2.X, p2.Z);
	m_ell_addh_t(&e, &p3, &p, &p2, &p);
	gmp_printf("X3 = %Zd\t Z3 = %Zd\n", p3.X, p3.Z);
	m_ell_addh_t(&e, &p4a, &p3, &p, &p2);
	gmp_printf("X4A = %Zd\t Z4A = %Zd\n", p4a.X, p4a.Z);
	m_ell_duph_t(&e, &p4d, &p2);
	gmp_printf("X4D = %Zd\t Z4D = %Zd\n", p4d.X, p4d.Z);
}

void testmul(const char *n, const char *k_s)
{
	mpz_t k;
	m_ellc e;
	m_ellp p, p2, p3, p4, pk;
	struct timespec start, end, diff;
	gmp_randstate_t state;

	gmp_randinit_default(state);
	mpz_init_set_str(k, k_s, 10);
	m_ellc_init(&e);
	m_ellp_init(&p);
	m_ellp_init(&p2);
	m_ellp_init(&p3);
	m_ellp_init(&p4);
	m_ellp_init(&pk);

	mpz_set_str(e.n, n, 10);
	gmp_randseed_ui(state, time(NULL));

	if (m_ell_setrand2_t(e.n, e.C2, &p, state))
		error_msg("error in setrand at testmul\n");
	gmp_printf("c_2 = %Zd\nX1 = %Zd\t Z1 = %Zd\n", e.C2, p.X, p.Z);

	m_ell_duph_t(&e, &p2, &p);
	gmp_printf("X2 = %Zd\t Z2 = %Zd\n", p2.X, p2.Z);
	m_ell_addh_t(&e, &p3, &p, &p2, &p);
	gmp_printf("X3 = %Zd\t Z3 = %Zd\n", p3.X, p3.Z);
	m_ell_addh_t(&e, &p4, &p3, &p, &p2);
	gmp_printf("X4 = %Zd\t Z4 = %Zd\n", p4.X, p4.Z);

	get_current_time(&start);
	m_ell_mul_t(k, &e, &pk, &p);
	get_current_time(&end);
	timespec_diff(&start, &end, &diff);
	gmp_printf("Xk = %Zd\t Zk = %Zd\n", pk.X, pk.Z);
	printf("mul time: %lu[s]\t%lu[ns]\n", diff.tv_sec, diff.tv_nsec);
}

static void test_fact(FILE * log_t, FILE * log_f, unsigned long tent,
		      unsigned long digits, unsigned long b1, unsigned long b2,
		      gmp_randstate_t state)
{
	struct timespec start, end, diff, mean = {.tv_sec = 0,.tv_nsec = 0 };
	mpz_t p, q, n, rnd_rng1, rnd_rng2, fact[2];
	unsigned long tot_iter = 0;
	long res;

	if (tent == 0)
		return;

	mpz_inits(p, q, n, rnd_rng1, rnd_rng2, fact[0], fact[1], NULL);
	mpz_ui_pow_ui(rnd_rng1, 10, digits);
	mpz_ui_pow_ui(rnd_rng2, 10, digits + 10);
	mpz_mul_ui(rnd_rng1, rnd_rng1, 4);	// rand_range = 4*10^c
	mpz_mul_ui(rnd_rng2, rnd_rng2, 4);	// rand_range = 4*10^c

	
	for (unsigned int i = 0; i < tent; i++) {
		get_randprime(p, rnd_rng1, rnd_rng1, state);
		get_randprime(q, rnd_rng2, rnd_rng2, state);
		mpz_mul(n, p, q);

		get_current_time(&start);
		res = factorize(fact, n, b1, b2, MAX_ITER);
		get_current_time(&end);
		timespec_diff(&start, &end, &diff);
		timespec_sum(&mean, &mean, &diff);

		if (res != ELL_FACT_NOT_FOUND) {
			gmp_fprintf(log_f,
				    "n = %Zd\tf1 = %Zd\tf2 = %Zd\tFASE = %d\t ITER = %d\n",
				    n, fact[0], fact[1], ell_fact_FASE(res,
								       DEFAULT_COEFF),
				    ell_fact_ITER(res, MAX_ITER));
			tot_iter += ell_fact_ITER(res, MAX_ITER);
			if (!FACT_IS_CORRECT(p, q, fact[0], fact[1])) {
				fputs("FACTORIZATION NOT CORRECT!!\n\n\n",
				      log_f);
				error_msg("FACTORIZATION NOT CORRECT!!\n");
			}
		} else {
			tot_iter += MAX_ITER;
			fputs("FACT NOT FOUND\n", log_f);
			puts("FACT NOT FOUND\n");
		}
	}
	timespec_div(&mean, &mean, tent);
	tot_iter = tot_iter / tent;
	gmp_fprintf(log_t, "%lu,%lu,%lu,%lu,%lu.%lu\n", digits, b1, b2, tot_iter,
		    diff.tv_sec, diff.tv_nsec / 1000000);
	mpz_clears(p, q, rnd_rng1, rnd_rng2, n, fact[0], fact[1], NULL);
}

void test_diff(const char *n_s)
{
#define rep_size 8
	m_ellc e;
	m_ellp p, p2, p4, p6;
	mpz_rep beta;
	m_ellp_rep rep;
	gmp_randstate_t state;

	gmp_randinit_default(state);
	m_ellc_init(&e);
	m_ellp_init(&p);
	m_ellp_init(&p2);
	m_ellp_init(&p4);
	m_ellp_init(&p6);
	m_ellp_rep_init(&rep, rep_size);
	mpz_rep_init(&beta, rep_size);

	mpz_set_str(e.n, n_s, 10);
	gmp_randseed_ui(state, time(NULL));
	if (m_ell_setrand2_t(e.n, e.C2, &p, state))
		error_msg("error in setrand at test_diff\n");

	m_ell_duph_t(&e, &p2, &p);
	m_ell_duph_t(&e, &p4, &p2);
	m_ell_addh_t(&e, &p6, &p4, &p2, &p2);
	gmp_printf
	    ("2p: X = %Zd and Z = %Zd\n4p: X = %Zd and Z = %Zd\n6p: X = %Zd and Z = %Zd\n",
	     p2.X, p2.Z, p4.X, p4.Z, p6.X, p6.Z);

	m_ell_diff_t(&rep, &beta, &e, &p);
	for (int i = 0; i < rep_size; i++)
		gmp_printf("%lu * p = (%Zd, %Zd)\t beta = %Zd\n", (i + 1) * 2,
			   rep.p[i].X, rep.p[i].Z, beta.v[i]);
}

void test_prime_diff(const char *start_s, const char *end_s)
{
#define bo 100000000
	unsigned char *v;
	unsigned long start, end;
	mpz_temp temp;

	start = str2l(start_s);
	end = str2l(end_s);
	mpz_temp_init(&temp, 3);
	v = allocate(sizeof(unsigned char) * bo);

	get_prime_diff(start, 0, end, v, &temp);
	int i;
	for (i = 0; i < bo; i++) {
		if (!v[i])
			break;
		start = start + v[i];
		printf("prime : %lu diff : %u\n", start, v[i]);
	}
	printf("n_prime : %u\n", i);
}

void test_mul_ui(const char *n_s, const char *k_s)
{
	unsigned long k;
	mpz_t n;
	m_ellc e;
	m_ellp p, r;
	mpz_temp temp;
	m_ellp_temp p_temp;
	gmp_randstate_t state;

	mpz_init(n);
	gmp_randinit_default(state);
	m_ellc_init(&e);
	m_ellp_init(&p);
	m_ellp_init(&r);
	mpz_temp_init(&temp, n_temp_mul + 1);
	m_ellp_temp_init(&p_temp, n_p_temp_mul);

	mpz_set_str(n, n_s, 10);
	k = str2l(k_s);

	gmp_randseed_ui(state, time(NULL));
	mpz_set(e.n, n);

	m_ell_setrand2_t(e.n, e.C2, &p, state);
	check_mul_ui(k, &e, &r, &p, &p_temp, &temp);

	mpz_clear(n);
}

void test_mmul(const char *n_s, const char *b_s)
{
	unsigned long b;
	mpz_t n, k;
	m_ellc e;
	m_ellp p, r, r_m;
	mpz_temp temp;
	gmp_randstate_t state;
	struct timespec time1, time2, time3;

	mpz_inits(n, k, NULL);
	gmp_randinit_default(state);
	m_ellc_init(&e);
	m_ellp_init(&p);
	m_ellp_init(&r);
	m_ellp_init(&r_m);
	mpz_temp_init(&temp, n_temp_create_bigk);

	mpz_set_str(n, n_s, 10);
	b = str2l(b_s);

	gmp_randseed_ui(state, time(NULL));
	mpz_set(e.n, n);

	m_ell_setrand2_t(e.n, e.C2, &p, state);
	create_bigk(k, b, &temp);

	get_current_time(&time1);
	m_ell_mul_t(k, &e, &r, &p);
	get_current_time(&time2);
	timespec_diff(&time1, &time2, &time3);
	printf("ell_mul time: %lu[s]\t%lu[ns]\n", time3.tv_sec, time3.tv_nsec);

	get_current_time(&time1);
	mm_ell_mul_t(k, e.n, e.C2, &r_m, &p);
	get_current_time(&time2);
	timespec_diff(&time1, &time2, &time3);
	printf("ell_mulm time: %lu[s]\t%lu[ns]\n", time3.tv_sec, time3.tv_nsec);
	gmp_printf("R\nX %Zd\nZ %Zd\nR2\n X %Zd\n Z %Zd\n", r.X, r.Z, r_m.X,
		   r_m.Z);
	assert(mpz_cmp(r.X, r_m.X) == 0);
	assert(mpz_cmp(r.Z, r_m.Z) == 0);

	mpz_clears(n, k, NULL);
}

void test_mfase2(const char *n_s, const char *b1_s, const char *b2_s)
{
	unsigned long b1, b2;
	mpz_t n, g, g_r;
	m_ellc e, e_r;
	m_ellp p, p_r;
	mpz_rep beta;
	m_ellp_rep rep;
	unsigned char *vdiff;
	mpz_temp temp;
	m_ellp_temp p_temp;
	mform_data mdata;
	gmp_randstate_t state;
	struct timespec time1, time2, time3;

	mpz_inits(n, g, g_r, NULL);
	gmp_randinit_default(state);
	m_ellc_init(&e);
	m_ellc_init(&e_r);
	m_ellp_init(&p);
	m_ellp_init(&p_r);
	m_ellp_rep_init(&rep, FACT_REP_SIZE);
	mpz_rep_init(&beta, FACT_REP_SIZE);
	mpz_temp_init(&temp, n_temp_mfase2);
	m_ellp_temp_init(&p_temp, n_p_temp_mfase2);

	mpz_set_str(n, n_s, 10);
	b1 = str2l(b1_s);
	b2 = str2l(b2_s);

	const int vdiff_size = get_vdiff_size(b2);
	if (vdiff_size == -1)
		error_msg("b2 to big\n");
	if ((vdiff = malloc(vdiff_size)) == NULL)
		error_msg("error in malloc at m_ell_fact\n");
	get_prime_diff(b1, 1, b2, vdiff, &temp);

	gmp_randseed_ui(state, time(NULL));
	mpz_set(e.n, n);
	mpz_set_ui(g, 1);

	m_ell_setrand2_t(e.n, e.C2, &p, state);

	mform_data_init(&mdata, n, &temp);
	to_mform(g_r, g, &mdata, &temp);
	m_ellc_to_mform(&e_r, &e, &mdata, &temp);
	m_ellp_to_mform(&p_r, &p, &mdata, &temp);

	get_current_time(&time1);
	m_ell_diff(&rep, &beta, e.n, e.C2, &p, &temp);
	m_ell_fase2(g, b1, b2, e.n, e.C2, &p, rep, beta, vdiff, &p_temp, &temp);
	get_current_time(&time2);
	timespec_diff(&time1, &time2, &time3);
	printf("m_ell_fase2 time: %lu[s]\t%lu[ns]\n", time3.tv_sec,
	       time3.tv_nsec);

	get_current_time(&time1);
	mm_ell_diff(rep.p, beta.v, rep.lenght, e_r.C2, &p_r, &mdata, &temp);
	mm_ell_fase2(g_r, b1, b2, e_r.C2, &p_r, rep.p, beta.v, rep.lenght,
		     vdiff, &mdata, &p_temp, &temp);
	get_current_time(&time2);
	timespec_diff(&time1, &time2, &time3);
	printf("mm_ell_fase2 time: %lu[s]\t%lu[ns]\n", time3.tv_sec,
	       time3.tv_nsec);
	from_mform(g_r, g_r, &mdata, &temp);

	gmp_printf("g %Zd\ng_r %Zd\n", g, g_r);
	assert(mpz_cmp(g, g_r) == 0);

	mpz_clears(n, g, g_r, NULL);
}
