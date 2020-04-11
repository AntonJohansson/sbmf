#include "log.h"
#include <stdio.h>
#include <stdarg.h>

log_item items[MAX_ITEMS] = {0};
log_size_t item_top = 0;
log_size_t item_bottom = 0;

static void log_impl(log_mode mode, const char* fmt, va_list args) {
	log_item* item = &items[item_top];
	item->mode = mode;

	int len = vsnprintf(item->buf, MAX_ITEM_BUF_LEN, fmt, args);
	item->buf[len] = 0;

	item_top = (item_top + 1) % MAX_ITEMS;
	//if (item_top >= item_bottom) {
	//	item_bottom = (item_top + 1) % MAX_ITEMS;
	//}

		FILE* fd = stdout;
	switch (item->mode) {
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
		fputs(item->buf, fd);
		fputs("\n", fd);
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

void log_mirror_to_fd(FILE* fd) {
	for (log_size_t i = item_bottom; i < item_top; i = (i+1) % MAX_ITEMS) {
		log_item* item = &items[i];
		switch (item->mode) {
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
		fputs(item->buf, fd);
		fputs("\n", fd);
	}
}
