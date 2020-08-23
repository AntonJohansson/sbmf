#include "sbmf.h"
#include <sbmf/memory/stack_allocator.h>
#include <sbmf/debug/log.h>
#include <sbmf/debug/profile.h>
#include <signal.h>

static struct {
	struct stack_allocator* main_stack;
	bool initialized;
} _state;

static void inthandler(int dummy) {
	sbmf_shutdown();
	exit(1);
}

void sbmf_init() {
	signal(SIGINT, inthandler);
	log_set_level(LOG_LEVEL_WARNING);
	log_open_and_clear_file("sbmf.log");
	_state.main_stack = sa_make(32*1024*1024);
	_state.initialized = true;
}

void sbmf_shutdown() {
	log_close_file();
	sa_destroy(_state.main_stack);
	//profile_print_results();
	_state.initialized = false;
}

u8* sbmf_stack_push(u32 size_in_bytes) {
	return sa_push(_state.main_stack, size_in_bytes);
}

u32 sbmf_stack_marker() {
	return _state.main_stack->top;
}

void sbmf_stack_free_to_marker(u32 marker) {
	_state.main_stack->top = marker;
}
