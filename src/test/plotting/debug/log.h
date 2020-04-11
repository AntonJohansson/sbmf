#pragma once
#include <stdint.h>
#include <stdio.h>

#define MAX_ITEM_BUF_LEN 1024
#define MAX_ITEMS 128
typedef uint8_t log_size_t;

typedef enum {
	LOG_MODE_INFO = 0,
	LOG_MODE_WARNING,
	LOG_MODE_ERROR,
} log_mode;

typedef struct {
	log_mode mode;
	char buf[MAX_ITEM_BUF_LEN];
} log_item;

extern log_item items[MAX_ITEMS];
extern log_size_t item_top;
extern log_size_t item_bottom;

extern void log_info(const char* fmt, ...);
extern void log_warning(const char* fmt, ...);
extern void log_error(const char* fmt, ...);

extern void log_mirror_to_fd(FILE* fd);
