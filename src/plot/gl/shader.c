#include "shader.h"
#include <assert.h>

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <assert.h>

static GLuint compile_shader(const char* buffer, GLenum shader_type){
	GLuint handle = glCreateShader(shader_type);
	// NULL is passed as last parameter as buffer is zero-terminated
	// from read_file_to_buffer(...)
	glShaderSource(handle, 1, (const GLchar**)&buffer, NULL);
	glCompileShader(handle);

	GLint succesful = 0;
	glGetShaderiv(handle, GL_COMPILE_STATUS, &succesful);
	if(!succesful){
		GLint log_size;
		glGetShaderiv(handle, GL_INFO_LOG_LENGTH, &log_size);

		char* log_buffer = (char*)malloc(log_size);
		glGetShaderInfoLog(handle, log_size, &log_size, log_buffer);

		fprintf(stderr, "Shader compile error: %s\n for buffer:\n%s\n", log_buffer, buffer);

		free(log_buffer);

		glDeleteShader(handle);
	}

	return handle;
}



ShaderProgram program_make() {
	return (ShaderProgram) {
		.handle = glCreateProgram()
	};
}

void program_attach_shader_from_buffer(ShaderProgram* program, const char* buffer, GLenum type) {
	assert(program->current_shader_index + 1 < MAX_SHADERCOUNT);
	program->shaders[program->current_shader_index] = compile_shader(buffer, type);
	glAttachShader(program->handle, program->shaders[program->current_shader_index]);
	program->current_shader_index++;
}

extern void program_link(ShaderProgram* program) {
	glLinkProgram(program->handle);

	GLint link_succesful = 0;
	glGetProgramiv(program->handle, GL_LINK_STATUS, &link_succesful);
	if (!link_succesful) {
		GLint log_size;
		glGetProgramiv(program->handle, GL_INFO_LOG_LENGTH, &log_size);

		char log_buffer[log_size+1];
		glGetProgramInfoLog(program->handle, log_size, &log_size, log_buffer);
		log_buffer[log_size] = 0;

		fprintf(stderr, "Program link error (size: %i):%s\n", log_size+1, log_buffer);

		glDeleteProgram(program->handle);
		for (int i = 0; i < program->current_shader_index; i++)
			glDeleteShader(program->shaders[i]);
	}

	for (size_t i = 0; i < program->current_shader_index; i++) {
		glDetachShader(program->handle, program->shaders[i]);
		glDeleteShader(program->shaders[i]);
	}
}

void program_bind_fragdata_location(ShaderProgram* program, const char* loc){
	glBindFragDataLocation(program->handle, 0, loc);
}

void program_bind_vertex_attribute(VertexAttribute attrib) {
	const GLvoid* offset_ptr = (const char*)(0) + attrib.offset;
	glEnableVertexAttribArray(attrib.location);
	glVertexAttribPointer(attrib.location, attrib.size, attrib.type,
												attrib.normalized, attrib.stride, offset_ptr);
}

void program_bind_vertex_format(VertexAttribute format[], unsigned int size) {
	for (uint8_t i = 0; i < size; ++i) {
		program_bind_vertex_attribute(format[i]);
	}
}
