#pragma once

#include <glad/glad.h>
#include <stdio.h>

#define GL_DUMP_ERRORS() 														\
	do {  																						\
		GLenum e; 																			\
		while((e = glGetError()) != GL_NO_ERROR) {			\
			fprintf(stderr, "GL: %s", glGetString(e)); 		\
		}																								\
	}while(0)

extern void setup_gl_error_callback();
