#pragma once

// Dynamic library loading
extern void* plt_dlopen(const char* name);
extern void plt_dlsym(void** var, void* handle, const char* sym);
extern void plt_dlclose(void* handle);
extern void plt_shellcmd(const char* command);
