#pragma once

#ifdef SBMF_MALLOC
	#include <stdlib.h> // for malloc
	
	#define SBMF_MALLOC(size_in_bytes) \
		malloc(size_in_bytes)
#endif
