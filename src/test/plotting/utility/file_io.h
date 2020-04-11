#pragma once

#include <stdio.h>
#include <commonstd.h>
#include <debug/log.h>
#include <stdlib.h>

static char* read_file_to_buffer(const char filename[]){
	FILE* fd = fopen(filename, "rb");
	if(!fd){
		log_error("Loading file %s not found", filename);
		return 0;
	}

	fseek(fd, 0, SEEK_END);
	size_t size = ftell(fd);
	char* buffer = (char*)malloc(size+1);

	fseek(fd, 0, SEEK_SET);
	fread(buffer, size, 1, fd);
	fclose(fd);

	buffer[size] = 0;

	return buffer;
}
