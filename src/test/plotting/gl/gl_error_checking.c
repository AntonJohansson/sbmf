#include "gl_error_checking.h"
#include <glad/glad.h>
#include <string.h>
#include <stdio.h>

static const char* APIENTRY gl_error_source_to_string(GLenum source){
  switch(source){
  case GL_DEBUG_SOURCE_API:
    return "api"; break;
  case GL_DEBUG_SOURCE_WINDOW_SYSTEM:
    return "window system"; break;
  case GL_DEBUG_SOURCE_SHADER_COMPILER:
    return "shader compiler"; break;
  case GL_DEBUG_SOURCE_APPLICATION:
    return "shader compiler"; break;
  case GL_DEBUG_SOURCE_THIRD_PARTY:
    return "third party"; break;
  case GL_DEBUG_SOURCE_OTHER:
    return "other"; break;
  default:
    return "unknown source";
  }
}

static const char* APIENTRY gl_error_type_to_string(GLenum type){
  switch(type){
  case GL_DEBUG_TYPE_ERROR:
    return "ERROR"; break;
  case GL_DEBUG_TYPE_UNDEFINED_BEHAVIOR:
    return "undefined behaviour"; break;
  case GL_DEBUG_TYPE_DEPRECATED_BEHAVIOR:
    return "deprecated behaviour"; break;
  case GL_DEBUG_TYPE_PORTABILITY:
    return "portability"; break;
  case GL_DEBUG_TYPE_PERFORMANCE:
    return "performance"; break;
  case GL_DEBUG_TYPE_MARKER:
    return "marker"; break;
  case GL_DEBUG_TYPE_PUSH_GROUP:
    return "push group"; break;
  case GL_DEBUG_TYPE_POP_GROUP:
    return "pop group"; break;
  case GL_DEBUG_TYPE_OTHER:
    return "other"; break;
  default:
    return "unknown message type";
  };
}

static const char* APIENTRY gl_error_severity_to_string(GLenum severity){
  switch(severity){
  case GL_DEBUG_SEVERITY_NOTIFICATION: 	return "\033[1;32mnotification\033[0m"; break;
  case GL_DEBUG_SEVERITY_LOW: 					return "\033[1;34mlow\033[0m"; break;
  case GL_DEBUG_SEVERITY_MEDIUM: 				return "\033[1;33mmedium\033[0m"; break;
  case GL_DEBUG_SEVERITY_HIGH: 					return "\033[1;31mhigh\033[0m"; break;
  default: 															return "unknown severity";
  };
}

static void APIENTRY gl_error_callback( GLenum source,
                        GLenum type,
                        GLuint id,
                        GLenum severity,
                        GLsizei length,
                        const GLchar* message,
                        const void* userParam ){
  (void)id;
  (void)length;
  (void)userParam;

  if(severity == GL_DEBUG_SEVERITY_NOTIFICATION)return;
	fprintf(stderr, "OpenGL debug message\nsource: %s\ntype: %s\nseverity: %s\nmessage: %s\n",
			gl_error_source_to_string(source),
			gl_error_type_to_string(type),
			gl_error_severity_to_string(severity),
			message);
}

void setup_gl_error_callback(){
	// TODO: check if this functionality is available
	glEnable(GL_DEBUG_OUTPUT);
  glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS);
  glDebugMessageControl(GL_DONT_CARE, GL_DONT_CARE, GL_DONT_CARE, 0, NULL, GL_TRUE);
	glDebugMessageCallback((GLDEBUGPROC) gl_error_callback, 0);
}
