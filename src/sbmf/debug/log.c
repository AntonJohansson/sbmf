#include "log.h"
#include <sbmf/common/common.h>
#include <stdio.h>
#include <stdarg.h>
#include <errno.h>
#include <string.h>

#define MAX_LOG_LEN 256
#define LOG_BUFFER_LEN 16

static struct {
	FILE* fd;
	enum log_level current_log_level;
	char buffer[LOG_BUFFER_LEN][MAX_LOG_LEN];
	u32 current_log_entry;
} _log_state = {
	.fd = 0,
	.current_log_level = LOG_LEVEL_INFO,
	.buffer = {{0}},
	.current_log_entry = 0,
};

static void dump_buffer() {
	for (u32 i = 0; i < _log_state.current_log_entry; ++i) {
		fprintf(_log_state.fd, "%s\n", &_log_state.buffer[i][0]);
	}
	_log_state.current_log_entry = 0;
}

void log_set_level(enum log_level level) {
	_log_state.current_log_level = level;
}

void log_open_and_clear_file(const char* file) {
	if (_log_state.fd != 0)
		log_close_file();

	_log_state.fd = fopen(file, "w");
	if (!_log_state.fd) {
		log_error("Unable to open log file: %s", file);
		log_error("errno: (%d) %s", errno, strerror(errno));
	}
}

void log_close_file() {
	if (_log_state.fd == 0) {
		log_error("Trying to close unopened log file!");
		return;
	}

	if (_log_state.current_log_entry > 0)
		dump_buffer();

	fclose(_log_state.fd);
	_log_state.fd = 0;
}

static void log_impl(enum log_level level, const char* fmt, va_list args) {
	if (_log_state.fd != 0) {
		if (_log_state.current_log_entry >= LOG_BUFFER_LEN) {
			dump_buffer();
		}
	} else {
		_log_state.current_log_entry = 0;
	}

	u32 i = _log_state.current_log_entry;
	static const char* info_str 	= "\033[1;34m[info]\033[0m ";
	static const char* warning_str 	= "\033[1;33m[warning]\033[0m ";
	static const char* error_str 	= "\033[1;31m[error]\033[0m ";
	static const char* unknown_str 	= "\033[1;31m[unknown]\033[0m ";

	const char* level_str = 0;
	switch (level) {
		case LOG_LEVEL_INFO:
			level_str = info_str;
			break;
		case LOG_LEVEL_WARNING:
			level_str = warning_str;
			break;
		case LOG_LEVEL_ERROR:
			level_str = error_str;
			break;
		default:
			level_str = unknown_str;
	};

	u32 level_str_len = strlen(level_str);
	strncpy(&_log_state.buffer[i][0], level_str, level_str_len);

	vsnprintf(&_log_state.buffer[i][level_str_len], MAX_LOG_LEN-level_str_len, fmt, args);
	//_log_state.buffer[i][bytes] = 0;

	if (_log_state.current_log_level <= level) {
		printf("%s\n", &_log_state.buffer[i][0]);
	}

	_log_state.current_log_entry++;
}

void log_info(const char* fmt, ...) {
	va_list args;
	va_start(args, fmt);
	log_impl(LOG_LEVEL_INFO, fmt, args);
	va_end(args);
}
void log_warning(const char* fmt, ...) {
	va_list args;
	va_start(args, fmt);
	log_impl(LOG_LEVEL_WARNING, fmt, args);
	va_end(args);
}
void log_error(const char* fmt, ...) {
	va_list args;
	va_start(args, fmt);
	log_impl(LOG_LEVEL_ERROR, fmt, args);
	va_end(args);
}
