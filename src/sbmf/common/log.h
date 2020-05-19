#pragma once
#include "common.h"

typedef enum {
	LOG_MODE_INFO = 0,
	LOG_MODE_WARNING,
	LOG_MODE_ERROR,
} log_mode;

extern void log_info(const char* fmt, ...);
extern void log_warning(const char* fmt, ...);
extern void log_error(const char* fmt, ...);
