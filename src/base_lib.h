#ifndef BASE_LIB_H_		/* Include guard */
#define BASE_LIB_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MINIMAL
#define TIME_PRECISION 1000000000

static inline void error_msg(const char *msg)
{
	perror(msg);
	exit(EXIT_FAILURE);
}

static inline void *allocate(size_t size)
{
	void *p = malloc(size);
	if (p == NULL) {
		perror("error in malloc at ell_point_fact");
		return NULL;
	}
	memset(p, 0, size);
	return p;
}

void get_arg(const char *coded, unsigned long *start, unsigned long *inc, unsigned int *size);

void get_current_time(struct timespec *curr_t);

static inline void timespec_div(struct timespec *rop, const struct timespec *op, unsigned long div)
{
	rop->tv_nsec = ((op->tv_sec % div) * TIME_PRECISION / div) + (op->tv_nsec / div);
	rop->tv_sec = op->tv_sec / div;
}

static inline void timespec_sum(struct timespec *rop, const struct timespec *op1, const struct timespec *op2)
{
	if (op1->tv_nsec + op2->tv_nsec >= TIME_PRECISION) {
		rop->tv_sec = op1->tv_sec + op2->tv_sec + 1;
		rop->tv_nsec = op1->tv_nsec + op2->tv_nsec - TIME_PRECISION;
	} else {
		rop->tv_sec = op1->tv_sec + op2->tv_sec;
		rop->tv_nsec = op1->tv_nsec + op2->tv_nsec;
	}
}

static inline void timespec_diff(const struct timespec *start, const struct timespec *stop, struct timespec *result)
{
	if ((stop->tv_nsec - start->tv_nsec) < 0) {
		result->tv_sec = stop->tv_sec - start->tv_sec - 1;
		result->tv_nsec = stop->tv_nsec - start->tv_nsec + TIME_PRECISION;
	} else {
		result->tv_sec = stop->tv_sec - start->tv_sec;
		result->tv_nsec = stop->tv_nsec - start->tv_nsec;
	}
}

#define printf_timespec(s, ts) printf("%s\t %lu [s]\t %lu[ns]\n", s, (ts)->tv_sec, (ts)->tv_nsec)

#ifdef MINIMAL
#define STAMP(...)
#else
#define STAMP(...) printf(__VA_ARGS__)
#endif

#ifdef DEBUG
#define GMP_DEBUG_STAMP(...) gmp_printf(__VA_ARGS__)
#else
#define GMP_DEBUG_STAMP(...)
#endif

#ifdef DEBUG
#define DEBUG_STAMP(...) printf(__VA_ARGS__)
#else
#define DEBUG_STAMP(...)
#endif

#endif //BASE_LIB_H_
