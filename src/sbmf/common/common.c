#include "common.h"
#include <sbmf/memory/stack_allocator.h>
#include <sbmf/debug/log.h>

sbmf_state _sbmf;

void sbmf_init() {
	//log_set_level(LOG_LEVEL_WARNING);
	log_open_and_clear_file("sbmf.log");
	_sbmf.main_stack = sa_make(32*1024*1024);
}

void sbmf_shutdown() {
	log_close_file();
	sa_destroy(_sbmf.main_stack);
}
