#include "profile.h"
#include "log.h"
#include <sbmf/types.h>

#include <time.h>
#include <errno.h>
#include <string.h>

#if PROFILE_ENABLE
	#define PROFILE_MAX_ENTRIES 100
	#define PROFILE_MAX_NAME_LEN 100
	#define PROFILE_MAX_SAMPLE_COUNT 100

	struct profile_entry {
		char name[PROFILE_MAX_NAME_LEN];

		struct timespec start;
		struct timespec end;

		struct timespec elapsed[PROFILE_MAX_SAMPLE_COUNT];
		u32 sample_count;

		f64 total_ms;
	};

	static u32 profile_entry_count = 0;
	static struct profile_entry profile_entries[PROFILE_MAX_ENTRIES] = {0};
#endif

static inline struct profile_entry* profile_find_entry_by_name(char const name[]) {
	for (int i = 0; i < profile_entry_count; ++i) {
		if(strcmp(profile_entries[i].name, name) == 0) {
			return &profile_entries[i];
		}
	}

	return 0;
}

void profile_begin_impl(char const name[]) {
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

void profile_end_impl(char const name[]) {
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
	//log_info("%s: %lu", name, elapsed->tv_nsec);
	entry->total_ms += (1000000000*elapsed->tv_sec + elapsed->tv_nsec)/(f64)1000000;
	entry->sample_count++;
}

void profile_print_results_impl() {
	if (profile_entry_count == 0)
		return;

	log_info("Timing results:");
	log_info("%40s | %10s | %10s | %10s | %15s ", "entry", "avg. (ns)","min (ns)", "max (ns)", "tot (ms)");
	log_info("-----------------------------------------+------------+------------+------------+----------------");
	for (u32 i = 0; i < profile_entry_count; ++i) {
		struct profile_entry entry = profile_entries[i];

		u64 avg = 0;
		u32 size = fmin(entry.sample_count, PROFILE_MAX_SAMPLE_COUNT);
		u64 ns_min = -1;
		u64 ns_max = 0;
		for (u32 j = 0; j < size; ++j) {
			u64 ns = 1000000000*entry.elapsed[j].tv_sec + entry.elapsed[j].tv_nsec;
			if (ns < ns_min)
				ns_min = ns;
			else if (ns > ns_max)
				ns_max = ns;
			avg += 	ns;
		}
		avg /= size;

		log_info("%40s | %10ld | %10ld | %10ld | %10.5lf ", entry.name, avg, ns_min, ns_max, entry.total_ms);
	}
}
