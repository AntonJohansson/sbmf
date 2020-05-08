#pragma once

#include <glad/glad.h>

static inline void bind_texture(GLuint id, unsigned int index){
	glActiveTexture(GL_TEXTURE0+index);
	glBindTexture(GL_TEXTURE_2D, id);
}

static inline GLuint create_texture(GLuint index, unsigned int w, unsigned int h, GLuint mode, GLuint type){
	GLuint texture;
	glGenTextures(1, &texture);
	bind_texture(texture, index);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	(void)type;
	(void)mode;
	(void)h;
	(void)w;
	glTexImage2D(GL_TEXTURE_2D, 0, mode, w, h, 0, mode, type, NULL);

	return texture;
}

//glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, w, h, 0, GL_RGBA, GL_FLOAT, NULL);
