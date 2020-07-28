#include "profile.h"
#include <sbmf/debug/log.h>
#include <errno.h>
#include <string.h>

#if PROFILE_ENABLE
	u32 profile_entry_count = 0;
	struct profile_entry profile_entries[PROFILE_MAX_ENTRIES] = {0};
#endif

static inline struct profile_entry* profile_find_entry_by_name(char const name[]) {
	for (int i = 0; i < profile_entry_count; ++i) {
		if(strcmp(profile_entries[i].name, name) == 0) {
			return &profile_entries[i];
		}
	}

	return 0;
}

void profile_begin(char const name[]) {
	struct profile_entry* entry = profile_find_entry_by_name(name);
	if(!entry) {
		if (profile_entry_count >= PROFILE_MAX_ENTRIES)
			return;

		entry = &profile_entries[profile_entry_count];
		u32 len = strlen(name)+1;
		u32 size = (len > PROFILE_MAX_NAME_LEN) ? PROFILE_MAX_NAME_LEN : len;
		memcpy(entry->name, name, size);
		profile_entry_count++;
	}

	if(clock_gettime(CLOCK_REALTIME, &entry->start) != 0) {
		log_error("clock_gettime(): failed in profile_begin() [errno: %s]", strerror(errno));
		entry->start = (struct timespec){0};
	}
}

void profile_end(char const name[]) {
	struct profile_entry* entry = profile_find_entry_by_name(name);
	if(!entry) {
		log_error("profile_end() failed, entry %s not found", name);
		return;
	}

	if(clock_gettime(CLOCK_REALTIME, &entry->end) != 0) {
		log_error("clock_gettime(): failed in profile_end() [errno: %s]", strerror(errno));
		entry->end = (struct timespec){0};
		return;
	}

	struct timespec* elapsed = &entry->elapsed[entry->sample_count % PROFILE_MAX_SAMPLE_COUNT];
	elapsed->tv_sec  = entry->end.tv_sec - entry->start.tv_sec;
	elapsed->tv_nsec = entry->end.tv_nsec - entry->start.tv_nsec;
	entry->total_ms += (1000000000*elapsed->tv_sec + elapsed->tv_nsec)/(f64)1000000;
	entry->sample_count++;
}

void profile_print_results_impl() {
	if (profile_entry_count == 0)
		return;

	log_info("Timing results:");
	log_info("%20s | %10s | %15s ", "entry", "avg. (us)", "tot. (ms)");
	log_info("---------------------+------------+----------------");
	for (size_t i = 0; i < profile_entry_count; ++i) {
		struct profile_entry entry = profile_entries[i];

		long avg = 0;
		u32 size = fmin(entry.sample_count, PROFILE_MAX_SAMPLE_COUNT);
		for (size_t j = 0; j < size; ++j) {
			avg += 1000000000*entry.elapsed[j].tv_sec + entry.elapsed[j].tv_nsec;
		}
		avg /= size;

		log_info("%20s | %10ld | %.5lf", entry.name, avg/1000, entry.total_ms);
	}
}
