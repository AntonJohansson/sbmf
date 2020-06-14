#include "common.h"
#include <sbmf/memory/stack_allocator.h>

sbmf_state _sbmf;

void sbmf_init() {
	_sbmf.main_stack = sa_make(32*1024*1024);
}

void sbmf_shutdown() {
	sa_destroy(_sbmf.main_stack);
}
