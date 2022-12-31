#include "base_lib.h"

void error_msg(const char *msg)
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

void get_arg(const char *coded, unsigned long *start, unsigned long *inc,
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

