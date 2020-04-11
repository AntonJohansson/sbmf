#pragma once

#include <glad/glad.h>
#include <commonstd.h>

extern GLuint program_begin();
extern void program_end();
extern void program_link();

extern void program_attach_shader_from_file(const char* file, GLenum type);
extern void program_attach_shader_from_buffer(const char* buffer, GLenum type);

typedef struct {
	GLuint location;
	GLint size;
	GLenum type;
	GLboolean normalized;
	GLsizei stride;
	GLsizei offset;
} vertex_attribute;

static vertex_attribute format_xy_uv[] = {
	{0,2,GL_FLOAT,GL_FALSE,4*sizeof(GL_FLOAT),0}, // position
	{1,2,GL_FLOAT,GL_FALSE,4*sizeof(GL_FLOAT),2*sizeof(GL_FLOAT)}, // texcoord
};

static vertex_attribute format_xy[] = {
	{0,2,GL_FLOAT,GL_FALSE,2*sizeof(GL_FLOAT),0}, // position
};

extern void bind_vertex_attribute(vertex_attribute attrib);
extern void bind_vertex_format(vertex_attribute format[], uint8_t size);
