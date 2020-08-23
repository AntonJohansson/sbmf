#pragma once

enum log_level {
	LOG_LEVEL_INFO = 0,
	LOG_LEVEL_WARNING = 1,
	LOG_LEVEL_ERROR = 2,
};

void log_set_level(enum log_level level);

/* Open/close the passed in log file, only one log file is
 * kept track of.
 */
void log_open_and_clear_file(const char* file);
void log_close_file();

/* Same interface as printf. log_<level> logs the passed in message
 * at log level <level>. A message will be mirrored to stdout
 * provided that the set log level LOG_LEVEL_<level> is less strict
 * than the called function, e.g. if LOG_LEVEL_WARNING is used
 * log_warning(...) and log_error(...) will be mirrored to stdout,
 * whilst log_info(...) will only show up in a log file if such a
 * file is opened.
 */
void log_info(const char* fmt, ...);
void log_warning(const char* fmt, ...);
void log_error(const char* fmt, ...);
void log_panic(const char* fmt, ...);
