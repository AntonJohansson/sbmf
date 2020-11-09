#include <sbmf/sbmf.h>
#include <sbmf/memory/stack_allocator.h>
#include <sbmf/debug/log.h>
#include <sbmf/debug/profile.h>

#include <omp.h>

#include <signal.h>
#include <stdio.h>
#include <time.h>

struct thread_local_storage {
	struct stack_allocator* stack;
};

static struct {
	struct thread_local_storage* thread_storage;
	u32 thread_count;

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
	if (_state.initialized) {
		log_error("sbmf_init(): Already intialized!");
		return;
	}

	log_info("Initializing");

	signal(SIGINT,  interrupt_handler);
	signal(SIGABRT, interrupt_handler);

	/* Setup logging */
	{
		/* log_set_level(LOG_LEVEL_WARNING); */
		log_open_and_clear_file("output/sbmf.log");

		_state.memory_log_fd = fopen("output/memory.log", "w");
		if (!_state.memory_log_fd) {
			log_error("Unable to open memory log file!");
		} else {
			fprintf(_state.memory_log_fd, "#time\ttop\tsize\n");
		}
	}

	/* Setup thread local storage */
	{
		_state.thread_count = (u32)omp_get_max_threads();
		log_info("\t- OpenMP using %u threads", _state.thread_count);

		_state.thread_storage = xmalloc(_state.thread_count*sizeof(struct thread_local_storage));
		for (u32 i = 0; i < _state.thread_count; ++i) {
			_state.thread_storage[i].stack = sa_make(32*1024*1024);
		}
	}

	if (clock_gettime(CLOCK_REALTIME, &_state.start_time) != 0) {
		log_error("sbmf_init(): clock_gettime failed!");
	}

	_state.initialized = true;
}

void sbmf_shutdown() {
	if (!_state.initialized) {
		log_error("sbmf_shutdown(): Not initialized!");
		return;
	}

	log_info("Shutting down");

	log_close_file();
	if (_state.memory_log_fd) {
		fclose(_state.memory_log_fd);
	}

	for (u32 i = 0; i < _state.thread_count; ++i) {
		sa_destroy(_state.thread_storage[i].stack);
	}

	_state.initialized = false;

	profile_print_results();
}

u8* sbmf_stack_push_impl(u32 size_in_bytes, const u32 linenumber, const char file[], const char func[]) {
	u32 threadid = (u32)omp_get_thread_num();
	assert(threadid < _state.thread_count);

	/* Dump memory usage after allocating */
	//if (_state.memory_log_fd) {
	//	f64 elapsed = elapsed_time();
	//	fprintf(_state.memory_log_fd, "%lf\t%u\t%u\t[%s:%d in %s()]\n",
	//			elapsed,
	//			_state.thread_storage[threadid].stack->top,
	//			_state.thread_storage[threadid].stack->size,
	//			file, linenumber, func);
	//}

	u8* ptr = sa_push(_state.thread_storage[threadid].stack, size_in_bytes);
	return ptr;
}

u32 sbmf_stack_marker() {
	u32 threadid = (u32)omp_get_thread_num();
	return _state.thread_storage[threadid].stack->top;
}

void sbmf_stack_free_to_marker(u32 marker) {
	u32 threadid = (u32)omp_get_thread_num();
	_state.thread_storage[threadid].stack->top = marker;
	/* Dump memory usage after freeing */
	//if (_state.memory_log_fd) {
	//	f64 elapsed = elapsed_time();

	//	fprintf(_state.memory_log_fd, "%lf\t%u\t%u\n",
	//			elapsed,
	//			_state.thread_storage[threadid].stack->top,
	//			_state.thread_storage[threadid].stack->size);
	//}
}
