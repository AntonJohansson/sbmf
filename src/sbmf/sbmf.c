#include <sbmf/sbmf.h>
#include <sbmf/memory/stack_allocator.h>
#include <sbmf/debug/profile.h>

#include <omp.h>

#include <signal.h>
#include <stdio.h>
#include <time.h>

struct thread_local_storage {
	struct stack_allocator* stack;
};

static struct {
	bool initialized;

	struct thread_local_storage* thread_storage;
	u32 thread_count;


	struct timespec start_time;

	sbmf_log_callback_func* log_callback;
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
		sbmf_log_error("sbmf_stack_push_impl(): clock_gettime failed!");
	}

	return elapsed_ns/((f64)1e9);
}

void sbmf_init() {
	if (_state.initialized) {
		sbmf_log_error("sbmf_init(): Already intialized!");
		return;
	}

	sbmf_log_info("Initializing");

	signal(SIGINT,  interrupt_handler);
	signal(SIGABRT, interrupt_handler);

	/* Setup thread local storage */
	{
		_state.thread_count = (u32)omp_get_max_threads();
		sbmf_log_info("\t- OpenMP using %u threads", _state.thread_count);

		_state.thread_storage = xmalloc(_state.thread_count*sizeof(struct thread_local_storage));
		for (u32 i = 0; i < _state.thread_count; ++i) {
			_state.thread_storage[i].stack = sa_make(32*1024*1024);
		}
	}

	if (clock_gettime(CLOCK_REALTIME, &_state.start_time) != 0) {
		sbmf_log_error("sbmf_init(): clock_gettime failed!");
	}

	_state.initialized = true;
}

void sbmf_shutdown() {
	if (!_state.initialized) {
		sbmf_log_error("sbmf_shutdown(): Not initialized!");
		return;
	}

	sbmf_log_info("Shutting down");

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
	//if (_state.memory_sbmf_log_fd) {
	//	f64 elapsed = elapsed_time();
	//	fprintf(_state.memory_sbmf_log_fd, "%lf\t%u\t%u\t[%s:%d in %s()]\n",
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
	//if (_state.memory_sbmf_log_fd) {
	//	f64 elapsed = elapsed_time();

	//	fprintf(_state.memory_sbmf_log_fd, "%lf\t%u\t%u\n",
	//			elapsed,
	//			_state.thread_storage[threadid].stack->top,
	//			_state.thread_storage[threadid].stack->size);
	//}
}







/*
 * Logging definitions
 */

#include <stdarg.h>

#define MAX_LOG_MSG_LEN 128

void sbmf_set_log_callback(sbmf_log_callback_func* func) {
	_state.log_callback = func;
}

static void sbmf_log(enum sbmf_log_level log_level, const char* fmt, va_list args) {
	if (!_state.log_callback)
		return;

	static char msg_buffer[MAX_LOG_MSG_LEN];

	vsnprintf(msg_buffer, MAX_LOG_MSG_LEN, fmt, args);
	_state.log_callback(log_level, msg_buffer);
}

void sbmf_log_info(const char* fmt, ...) {
	va_list args;
	va_start(args, fmt);
	sbmf_log(SBMF_LOG_LEVEL_INFO, fmt, args);
	va_end(args);
}

void sbmf_log_warning(const char* fmt, ...) {
	va_list args;
	va_start(args, fmt);
	sbmf_log(SBMF_LOG_LEVEL_WARNING, fmt, args);
	va_end(args);
}

void sbmf_log_error(const char* fmt, ...) {
	va_list args;
	va_start(args, fmt);
	sbmf_log(SBMF_LOG_LEVEL_ERROR, fmt, args);
	va_end(args);
}

void sbmf_log_panic(const char* fmt, ...) {
	va_list args;
	va_start(args, fmt);
	sbmf_log(SBMF_LOG_LEVEL_PANIC, fmt, args);
	va_end(args);
}
