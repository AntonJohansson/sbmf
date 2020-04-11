#pragma once

#include <debug/log.h>
#include <glad/glad.h>

#define GL_DUMP_ERRORS() 														\
	do {  																						\
		GLenum e; 																			\
		while((e = glGetError()) != GL_NO_ERROR) {			\
			log_error("GL: %s", glGetString(e)); 					\
		}																								\
	}while(0)

extern void setup_gl_error_callback();
