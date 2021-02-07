static f64 elapsed_time(struct timespec t0, struct timespec t1) {
	u64 elapsed_ns = (t1.tv_nsec - t0.tv_nsec) + (t1.tv_sec - t0.tv_sec)*(u64)1e9;
	return elapsed_ns/((f64)1e9);
}

static inline struct timespec current_time() {
	struct timespec t;
	if (clock_gettime(CLOCK_REALTIME, &t) != 0) {
		sbmf_log_error("clock_gettime failed!");
	}
	return t;
}
