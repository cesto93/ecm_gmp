#ifndef BASE_LIB_H_   /* Include guard */
#define BASE_LIB_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

//#define DEBUG
#define MINIMAL
#define TIME_PRECISION 1000000000

void error_msg(char *msg);
void *allocate(int size);

void get_arg(char *coded, unsigned long int *start, unsigned long int *inc, unsigned int *size);

void get_current_time(struct timespec *curr_t);
void timespec_diff(struct timespec *start, struct timespec *stop, struct timespec *result);
void timespec_sum(struct timespec *rop, struct timespec *op1, struct timespec *op2);
void timespec_div(struct timespec *rop, struct timespec *op, unsigned int div);

#define CHECK_ARGC(param_n) if (argc < (param_n + 1)) \
do { \
	perror("usage: ./test.o  mod \ntentative cifer_start:inc:size b_start:inc:size\n"); \
	exit(EXIT_FAILURE); \
} while (0)

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
