#pragma once

#include <glad/glad.h>

#define MAX_SHADERCOUNT 8

typedef struct {
	GLuint handle;

	GLuint current_shader_index;
	GLuint shaders[MAX_SHADERCOUNT];

} ShaderProgram;

typedef struct {
	GLuint location;
	GLint size;
	GLenum type;
	GLboolean normalized;
	GLsizei stride;
	GLsizei offset;
} VertexAttribute;

//static VertexAttribute format_xy_uv[] = {
//	{0,2,GL_FLOAT,GL_FALSE,4*sizeof(GL_FLOAT),0}, // position
//	{1,2,GL_FLOAT,GL_FALSE,4*sizeof(GL_FLOAT),2*sizeof(GL_FLOAT)}, // texcoord
//};
//
//static VertexAttribute format_xy[] = {
//	{0,2,GL_FLOAT,GL_FALSE,2*sizeof(GL_FLOAT),0}, // position
//};

extern ShaderProgram program_make();
extern void program_attach_shader_from_buffer(ShaderProgram* program, const char* buffer, GLenum type);
extern void program_link(ShaderProgram* program);
extern void program_bind_fragdata_location(ShaderProgram* program, const char* loc);
extern void program_bind_vertex_attribute(VertexAttribute attrib);
extern void program_bind_vertex_format(VertexAttribute format[], unsigned int size);
