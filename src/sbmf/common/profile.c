#include "profile.h"

size_t profile_used_data = 0;
profile_data_t profile_data[PROFILE_MAX_DATA] = {0};

static inline profile_data_t* profile_find_data_by_name(char const name[]) {
	for (int i = 0; i < profile_used_data; ++i) {
		if(strcmp(profile_data[i].name, name) == 0) {
			return &profile_data[i];
		}
	}

	return 0;
}

void profile_begin(char const name[]) {
	profile_data_t* data = profile_find_data_by_name(name);
	if(!data) {
		data = &profile_data[profile_used_data];
		size_t len = strlen(name)+1;
		size_t size = (len > PROFILE_MAX_NAME_LEN) ? PROFILE_MAX_NAME_LEN : len;
		memcpy(data->name, name, size);
		profile_used_data++;
	}

	if(clock_gettime(CLOCK_REALTIME, &data->start) != 0) {
		fprintf(stderr, "clock_gettime(...): failed in profile_begin(...) [errno: %s]\n", strerror(errno));
		data->start = (struct timespec){0};
	}
}

void profile_end(char const name[]) {
	profile_data_t* data = profile_find_data_by_name(name);
	if(!data) {
		fprintf(stderr, "profile_end(...) failed, profile data not found.\n");
		return;
	}

	if(clock_gettime(CLOCK_REALTIME, &data->end) != 0) {
		fprintf(stderr, "clock_gettime(...) failed in profile_end(...) [errno: %s]\n",
						strerror(errno));
		data->end = (struct timespec){0};
		return;
	}

	struct timespec* delta = &data->deltas[data->delta_num_samples];
	delta->tv_sec  = data->end.tv_sec - data->start.tv_sec;
	delta->tv_nsec = data->end.tv_nsec - data->start.tv_nsec;
	data->delta_num_samples = (data->delta_num_samples+1) % PROFILE_DATA_SAMPLE_LEN;
}

void profile_print_results_impl() {
	printf("Timing results:\n");
	for (size_t i = 0; i < profile_used_data; ++i) {
		profile_data_t data = profile_data[i];

		long avg = 0;
		for (size_t j = 0; j < PROFILE_DATA_SAMPLE_LEN; ++j) {
			avg += 1000000000*data.deltas[j].tv_sec + data.deltas[j].tv_nsec;
		}
		avg /= PROFILE_DATA_SAMPLE_LEN;

		printf("%10s -- %ld (ms)\n", data.name, avg/1000000);
	}
}
