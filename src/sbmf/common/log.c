#include "log.h"
#include <stdio.h>
#include <stdarg.h>

static void log_impl(log_mode mode, const char* fmt, va_list args) {
	FILE* fd = stdout;
	switch (mode) {
		case LOG_MODE_INFO:
			fputs("\033[1;34m[info] \033[0m ", fd);
			break;
		case LOG_MODE_WARNING:
			fputs("\033[1;33m[warning]\033[0m ", fd);
			break;
		case LOG_MODE_ERROR:
			fputs("\033[1;31m[error]\033[0m ", fd);
			break;
	};

	vfprintf(fd, fmt, args);
}

void log_info(const char* fmt, ...) {
	va_list args;
	va_start(args, fmt);
	log_impl(LOG_MODE_INFO, fmt, args);
	va_end(args);
}

void log_warning(const char* fmt, ...) {
	va_list args;
	va_start(args, fmt);
	log_impl(LOG_MODE_WARNING, fmt, args);
	va_end(args);
}

void log_error(const char* fmt, ...) {
	va_list args;
	va_start(args, fmt);
	log_impl(LOG_MODE_ERROR, fmt, args);
	va_end(args);
}
