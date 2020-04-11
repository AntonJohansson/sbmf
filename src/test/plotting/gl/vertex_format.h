#pragma once

#include <glad/glad.h>

typedef struct vertex_attribute{
	const char* name;
	GLenum type;
	size_t size;
} vertex_attribute;


typedef struct vertex_format{
} vertex_format;

extern vertex_format make_vertex_format();
extern void free_vertex_format(vertex_format* format);

extern void add_attribute(vertex_format* format, const char* name, GLenum type, size_t size)
