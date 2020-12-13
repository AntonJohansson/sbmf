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
}

#define sbmf_stack_push(size_in_bytes) \
	sbmf_stack_push_impl(size_in_bytes, __LINE__, __FILE__, __func__)

void* sbmf_stack_push_impl(u32 size_in_bytes, const u32 linenumber, const char file[], const char func[]) {
	u32 threadid = (u32)omp_get_thread_num();
	assert(threadid < _state.thread_count);
	SBMF_UNUSED(linenumber);
	SBMF_UNUSED(file);
	SBMF_UNUSED(func);

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
	return (void*)ptr;
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
