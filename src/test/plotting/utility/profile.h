#pragma once

#include <time.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include <debug/log.h>

#define PROFILE_MAX_DATA 100
#define PROFILE_MAX_NAME_LEN 100
#define PROFILE_DATA_SAMPLE_LEN 10

static inline void profile_begin(char const name[]);
static inline void profile_end(char const name[]);

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

typedef struct profile_data_t {
	char name[PROFILE_MAX_NAME_LEN];
	size_t delta_num_samples;

	struct timespec start;
	struct timespec end;
	struct timespec deltas[PROFILE_DATA_SAMPLE_LEN];
} profile_data_t;

extern size_t profile_used_data;
extern profile_data_t profile_data[PROFILE_MAX_DATA];

static inline profile_data_t* profile_find_data_by_name(char const name[]) {
	for(int i = 0; i < PROFILE_MAX_DATA; i++) {
		if(strcmp(profile_data[i].name, name) == 0) {
			return &profile_data[i];
		}
	}

	return 0;
}

static inline void profile_begin(char const name[]) {
	profile_data_t* data = profile_find_data_by_name(name);
	if(!data) {
		data = &profile_data[profile_used_data];
		size_t len = strlen(name)+1;
		size_t size = (len > PROFILE_MAX_NAME_LEN) ? PROFILE_MAX_NAME_LEN : len;
		memcpy(data->name, name, size);
		profile_used_data++;
	}

	if(clock_gettime(CLOCK_REALTIME, &data->start) != 0) {
		log_error("clock_gettime(...): failed in profile_begin(...) [errno: %s]\n", strerror(errno));
		data->start = (struct timespec){0};
	}
}

static inline void profile_end(char const name[]) {
	profile_data_t* data = profile_find_data_by_name(name);
	if(!data) {
		log_error("profile_end(...) failed, profile data not found.\n");
		return;
	}

	if(clock_gettime(CLOCK_REALTIME, &data->end) != 0) {
		log_error("clock_gettime(...) failed in profile_end(...) [errno: %s]\n",
						strerror(errno));
		data->end = (struct timespec){0};
		return;
	}

	struct timespec* delta = &data->deltas[data->delta_num_samples];
	delta->tv_sec  = data->end.tv_sec - data->start.tv_sec;
	delta->tv_nsec = data->end.tv_nsec - data->start.tv_nsec;
	data->delta_num_samples++;
	data->delta_num_samples %= PROFILE_DATA_SAMPLE_LEN;
}
