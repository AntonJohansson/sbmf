#pragma once

#include <time.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>

#define PROFILE_ENABLE 1
#define PROFILE_MAX_DATA 100
#define PROFILE_MAX_NAME_LEN 100
#define PROFILE_DATA_SAMPLE_LEN 10

typedef struct profile_data_t {
	char name[PROFILE_MAX_NAME_LEN];
	size_t delta_num_samples;

	struct timespec start;
	struct timespec end;
	struct timespec deltas[PROFILE_DATA_SAMPLE_LEN];
} profile_data_t;

extern size_t profile_used_data;
extern profile_data_t profile_data[PROFILE_MAX_DATA];

extern void profile_begin(char const name[]);
extern void profile_end(char const name[]);

extern void profile_print_results_impl();

#if PROFILE_ENABLE
	#define PROFILE_BEGIN(name)\
		profile_begin(name)
	#define PROFILE_END(name)\
		profile_end(name)

	// Returns the result of func_call via the
	// comma operator. Expressions such as:
	//
	// 	int x = PROFILE_FUNC(compute_x())
	// 
	// should therefore be possible.
	#define PROFILE_FUNC(func_call)\
		(PROFILE_BEGIN( #func_call ), func_call); PROFILE_END( #func_call )

	#define profile_print_results() \
		profile_print_results_impl()
#else
	#define PROFILE_BEGIN(name)
	#define PROFILE_END(name)
	#define PROFILE_FUNC(func_call)\
		func_call
	#define profile_print_results()
#endif
