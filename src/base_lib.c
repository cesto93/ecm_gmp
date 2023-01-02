#include "base_lib.h"

void get_arg(const char *coded, unsigned long *start, unsigned long *inc, unsigned int *size)
{
	if (sscanf(coded, "%lu:%lu:%u", start, inc, size) != 3)
		error_msg("error in getarg");
}

void get_current_time(struct timespec *curr_t)
{
	if (clock_gettime(CLOCK_REALTIME, curr_t) == -1)
		error_msg("error in gettime");
}
