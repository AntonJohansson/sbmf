#pragma once

#include "types.h"

#define SBMF_UNUSED(x) \
	(void)x

void sbmf_init();
void sbmf_shutdown();

/* Memory */
u8* sbmf_stack_push_impl(u32 size_in_bytes, const u32 linenumber, const char file[], const char func[]);

#define sbmf_stack_push(size_in_bytes) \
	sbmf_stack_push_impl(size_in_bytes, __LINE__, __FILE__, __func__)

u32  sbmf_stack_marker();
void sbmf_stack_free_to_marker(u32 marker);

/* Logging */

enum sbmf_log_level {
	SBMF_LOG_LEVEL_INFO    = 0,
	SBMF_LOG_LEVEL_WARNING = 1,
	SBMF_LOG_LEVEL_ERROR   = 2,
	SBMF_LOG_LEVEL_PANIC   = 3,
};

typedef void sbmf_log_callback_func(enum sbmf_log_level, const char*);

void sbmf_set_log_callback(sbmf_log_callback_func* func);

/* Common macros to simplify typing */
void sbmf_log_info(const char* fmt, ...);
void sbmf_log_warning(const char* fmt, ...);
void sbmf_log_error(const char* fmt, ...);
void sbmf_log_panic(const char* fmt, ...);
