#include "shader.h"
#include <utility/file_io.h>
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

//shaderprogram make_program(const char* vertex_file, const char* fragment_file){
//	char* vertex_buffer = read_file_to_buffer(vertex_file);
//	char* fragment_buffer = read_file_to_buffer(fragment_buffer);
//
//	shaderprogram program;
//	program.vertex_shader 	= compile_shader(vertex_buffer, GL_VERTEX_SHADER);
//	program.fragment_shader = compile_shader(fragment_shader, GL_FRAGMENT_SHADER);
//
//	progam.handle = glCreateProgram();
//	glAttachShader(progam.handle, program.vertex_shader);
//	glAttachShader(progam.handle, program.fragment_shader);
//
//	return progam;
//}

#define MAX_SHADERCOUNT 8
static GLuint _program;
static GLuint _shaders[MAX_SHADERCOUNT];
static size_t _shader_index = 0;

GLuint program_begin(){
	_program = glCreateProgram();
	return _program;
}

void program_end(){
	glBindFragDataLocation(_program, 0, "out_color");
	_shader_index = 0;
}

void program_link(){
	glLinkProgram(_program);

	GLint link_succesful = 0;
	glGetProgramiv(_program, GL_LINK_STATUS, &link_succesful);
	if(!link_succesful){
		GLint log_size;
		glGetProgramiv(_program, GL_INFO_LOG_LENGTH, &log_size);

		char* log_buffer = (char*)malloc(log_size+1);
		glGetProgramInfoLog(_program, log_size, &log_size, log_buffer);
		log_buffer[log_size] = 0;

		fprintf(stderr, "Program link error (size: %i):%s\n", log_size+1, log_buffer);

		free(log_buffer);

		glDeleteProgram(_program);
		for(int i = 0; i < MAX_SHADERCOUNT; i++)
			glDeleteShader(_shaders[i]);
	}

	for(size_t i = 0; i < _shader_index; i++){
		glDetachShader(_program, _shaders[i]);
		glDeleteShader(_shaders[i]);
	}
}

void program_attach_shader_from_file(const char* file, GLenum type){
	char* buffer = read_file_to_buffer(file);
	program_attach_shader_from_buffer(buffer, type);
	free(buffer);
}

void program_attach_shader_from_buffer(const char* buffer, GLenum type){
	assert(_shader_index+1 < MAX_SHADERCOUNT);
	_shaders[_shader_index] = compile_shader(buffer, type);
	glAttachShader(_program, _shaders[_shader_index]);
	++_shader_index;
}

void bind_vertex_attribute(vertex_attribute attrib) {
	const GLvoid* offset_ptr = (const char*)(0) + attrib.offset;
	glEnableVertexAttribArray(attrib.location);
	glVertexAttribPointer(attrib.location, attrib.size, attrib.type,
												attrib.normalized, attrib.stride, offset_ptr);
}

void bind_vertex_format(vertex_attribute format[], uint8_t size) {
	for (uint8_t i = 0; i < size; ++i) {
		bind_vertex_attribute(format[i]);
	}
}
