#include "sbmf.h"
#include <sbmf/memory/stack_allocator.h>
#include <sbmf/debug/log.h>
#include <sbmf/debug/profile.h>
#include <signal.h>
#include <stdio.h>

static struct {
	struct stack_allocator* main_stack;
	bool initialized;
	FILE* memory_log_fd;
} _state;

static void inthandler(int dummy) {
	sbmf_shutdown();
	exit(1);
}

void sbmf_init() {
	signal(SIGINT, inthandler);
	/* log_set_level(LOG_LEVEL_WARNING); */
	log_open_and_clear_file("sbmf.log");

	_state.memory_log_fd = fopen("memory.log", "w");
	if (!_state.memory_log_fd) {
		log_error("Unable to open memory log file!");
	} else {
		fprintf(_state.memory_log_fd, "#top\tsize\n");
	}

	_state.main_stack = sa_make(32*1024*1024);
	_state.initialized = true;
}

void sbmf_shutdown() {
	log_close_file();
	if (_state.memory_log_fd) {
		fclose(_state.memory_log_fd);
	}

	sa_destroy(_state.main_stack);
	_state.initialized = false;

	/* profile_print_results(); */
}

u8* sbmf_stack_push(u32 size_in_bytes) {
	u8* ptr = sa_push(_state.main_stack, size_in_bytes);
	if (_state.memory_log_fd) {
		fprintf(_state.memory_log_fd, "%u\t%u\n", _state.main_stack->top, _state.main_stack->size);
	}
	return ptr;
}

u32 sbmf_stack_marker() {
	return _state.main_stack->top;
}

void sbmf_stack_free_to_marker(u32 marker) {
	_state.main_stack->top = marker;
	if (_state.memory_log_fd) {
		fprintf(_state.memory_log_fd, "%u\t%u\n", _state.main_stack->top, _state.main_stack->size);
	}
}
