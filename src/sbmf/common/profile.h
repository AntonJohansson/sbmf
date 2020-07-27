#pragma once

#include <time.h>
#include <sbmf/common/common.h>

#define PROFILE_ENABLE 1
#define PROFILE_MAX_ENTRIES 100
#define PROFILE_MAX_NAME_LEN 100
#define PROFILE_MAX_SAMPLE_COUNT 10

struct profile_entry {
	char name[PROFILE_MAX_NAME_LEN];

	struct timespec start;
	struct timespec end;

	struct timespec elapsed[PROFILE_MAX_SAMPLE_COUNT];
	u32 sample_count;
};

void profile_begin(char const name[]);
void profile_end(char const name[]);
void profile_print_results_impl();

#if PROFILE_ENABLE
	extern u32 profile_entry_count;
	extern struct profile_entry profile_entries[PROFILE_MAX_ENTRIES];

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
