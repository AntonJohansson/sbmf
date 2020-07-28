#include "common.h"
#include <sbmf/memory/stack_allocator.h>
#include <sbmf/debug/log.h>
#include <sbmf/common/profile.h>
#include <signal.h>

sbmf_state _sbmf;

static void inthandler(int dummy) {
	sbmf_shutdown();
	exit(1);
}

void sbmf_init() {
	signal(SIGINT, inthandler);
	//log_set_level(LOG_LEVEL_WARNING);
	log_open_and_clear_file("sbmf.log");
	_sbmf.main_stack = sa_make(32*1024*1024);
	_sbmf.initialized = true;
}

void sbmf_shutdown() {
	log_close_file();
	sa_destroy(_sbmf.main_stack);
	profile_print_results();
	_sbmf.initialized = false;
}
