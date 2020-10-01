#include "sbmf.h"
#include <sbmf/memory/stack_allocator.h>
#include <sbmf/debug/log.h>
#include <sbmf/debug/profile.h>

#include <omp.h>

#include <signal.h>
#include <stdio.h>
#include <time.h>

#define THREAD_COUNT 4

static struct {
	//struct stack_allocator* main_stack;
	struct stack_allocator* stacks[THREAD_COUNT];
	bool initialized;
	FILE* memory_log_fd;

	struct timespec start_time;
} _state;

static void interrupt_handler(int dummy) {
	SBMF_UNUSED(dummy);
	sbmf_shutdown();
	exit(1);
}

static f64 elapsed_time() {
	struct timespec current_time;
	u64 elapsed_ns = 0;
	if (clock_gettime(CLOCK_REALTIME, &current_time) == 0) {
		elapsed_ns = (current_time.tv_nsec - _state.start_time.tv_nsec) + (current_time.tv_sec - _state.start_time.tv_sec)*1e9;
	} else {
		log_error("sbmf_stack_push_impl(): clock_gettime failed!");
	}

	return elapsed_ns/((f64)1e9);
}

void sbmf_init() {
	signal(SIGINT, interrupt_handler);
	signal(SIGABRT, interrupt_handler);
	/* log_set_level(LOG_LEVEL_WARNING); */
	log_open_and_clear_file("sbmf.log");

	_state.memory_log_fd = fopen("memory.log", "w");
	if (!_state.memory_log_fd) {
		log_error("Unable to open memory log file!");
	} else {
		fprintf(_state.memory_log_fd, "#time\ttop\tsize\n");
	}

	for (u32 i = 0; i < THREAD_COUNT; ++i) {
		_state.stacks[i] = sa_make(512*1024);
	}

	if (clock_gettime(CLOCK_REALTIME, &_state.start_time) != 0) {
		log_error("sbmf_init(): clock_gettime failed!");
	}

	_state.initialized = true;
}

void sbmf_shutdown() {
	log_close_file();
	if (_state.memory_log_fd) {
		fclose(_state.memory_log_fd);
	}

	for (u32 i = 0; i < THREAD_COUNT; ++i) {
		sa_destroy(_state.stacks[i]);
	}

	_state.initialized = false;

	/* profile_print_results(); */
}

u8* sbmf_stack_push_impl(u32 size_in_bytes, const u32 linenumber, const char file[], const char func[]) {
	i32 threadid = omp_get_thread_num();
	assert(threadid < THREAD_COUNT);

	/* Dump memory usage after allocating */
	if (_state.memory_log_fd) {
		f64 elapsed = elapsed_time();
		fprintf(_state.memory_log_fd, "%lf\t%u\t%u\t[%s:%d in %s()]\n", elapsed, _state.stacks[threadid]->top, _state.stacks[threadid]->size, file, linenumber, func);
	}

	u8* ptr = sa_push(_state.stacks[threadid], size_in_bytes);
	return ptr;
}

u32 sbmf_stack_marker() {
	i32 threadid = omp_get_thread_num();
	return _state.stacks[threadid]->top;
}

void sbmf_stack_free_to_marker(u32 marker) {
	i32 threadid = omp_get_thread_num();
	_state.stacks[threadid]->top = marker;
	/* Dump memory usage after freeing */
	if (_state.memory_log_fd) {
		f64 elapsed = elapsed_time();

		fprintf(_state.memory_log_fd, "%lf\t%u\t%u\n", elapsed, _state.stacks[threadid]->top, _state.stacks[threadid]->size);
	}
}
