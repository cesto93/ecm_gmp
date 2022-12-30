#include "base_lib.h"

void error_msg(char *msg)
{
	perror(msg);
	exit(EXIT_FAILURE);
}

void *allocate(int size)
{
	void *p;
	if ((p = malloc(size)) == NULL)
		error_msg("error in malloc at ell_point_fact");
	memset(p, 0, size);
	return p;
}

void get_arg(char *coded, unsigned long int *start, unsigned long int *inc,
	     unsigned int *size)
{
	if (sscanf(coded, "%lu:%lu:%u", start, inc, size) != 3)
		error_msg("error in getarg");
}

void get_current_time(struct timespec *curr_t)
{
	if (clock_gettime(CLOCK_REALTIME, curr_t) == -1)
		error_msg("error in gettime");
}

void timespec_diff(struct timespec *start, struct timespec *stop,
		   struct timespec *result)
{
	if ((stop->tv_nsec - start->tv_nsec) < 0) {
		result->tv_sec = stop->tv_sec - start->tv_sec - 1;
		result->tv_nsec =
		    stop->tv_nsec - start->tv_nsec + TIME_PRECISION;
	} else {
		result->tv_sec = stop->tv_sec - start->tv_sec;
		result->tv_nsec = stop->tv_nsec - start->tv_nsec;
	}
}

void timespec_div(struct timespec *rop, struct timespec *op, unsigned int div)
{
	rop->tv_nsec =
	    ((op->tv_sec % div) * TIME_PRECISION / div) + (op->tv_nsec / div);
	rop->tv_sec = op->tv_sec / div;
}

void timespec_sum(struct timespec *rop, struct timespec *op1,
		  struct timespec *op2)
{
	if (op1->tv_nsec + op2->tv_nsec >= TIME_PRECISION) {
		rop->tv_sec = op1->tv_sec + op2->tv_sec + 1;
		rop->tv_nsec = op1->tv_nsec + op2->tv_nsec - TIME_PRECISION;
	} else {
		rop->tv_sec = op1->tv_sec + op2->tv_sec;
		rop->tv_nsec = op1->tv_nsec + op2->tv_nsec;
	}
}
