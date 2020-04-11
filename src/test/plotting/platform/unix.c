#include "unix.h"
#include <debug/log.h>
#include <dlfcn.h>
#include <stdio.h>
#include <commonstd.h>

void* plt_dlopen(const char* name) {
	void* handle = dlopen(name, RTLD_NOW);
	if (!handle) {
		log_error("%s\n", dlerror());
		return 0;
	}

	return handle;
}

void plt_dlsym(void** var, void* handle, const char* sym) {
	*var = dlsym(handle, sym);
}

void plt_dlclose(void* handle) {
	dlclose(handle);
}

void plt_shellcmd(const char* command) {
	FILE* pd = popen(command, "w");
	if (!pd) {
		log_error("Could not execute command %s", command);
		return;
	}
	
	// Read entire file into buffer
	// and print contents if necessary
	{
		fseek(pd, 0, SEEK_END);

		size_t size = ftell(pd);

		if (size != (size_t)-1 && size > 0) {
			char buffer[size+1];

			fseek(pd, 0, SEEK_SET);
			fread(buffer, size, 1, pd);

			buffer[size] = 0;

			log_error("Output of command \"%s\":\n%s:", command, buffer);
		}
	}

	pclose(pd);
}
