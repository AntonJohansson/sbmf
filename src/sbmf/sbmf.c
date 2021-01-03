#include <sbmf/sbmf.h>

#define SBMF_UNUSED(x) \
	(void)x

#include <omp.h>
#include <cblas.h>
#include <lapacke.h>
#include <arpack/arpack.h>
/* #include <arpack/debug_c.h> */

/* c stdlib */
#include <stdio.h>  /* sprintf, fopen, fprintf, fread */
#include <stdlib.h> /* malloc */
#include <string.h> /* memcpy,memset */
#include <stdarg.h>
#include <assert.h>
#include <float.h>  /* isinf, isnan */

/* unix headers */
#include <signal.h>

/*
 * Getting into actual code
 */

void sbmf_log_info(const char* fmt, ...);
void sbmf_log_warning(const char* fmt, ...);
void sbmf_log_error(const char* fmt, ...);
void sbmf_log_panic(const char* fmt, ...);

#include "memory/stack_allocator.c"
#include "global_state.c"
#include "memory/bucketarray.c"
#include "memory/prioqueue.c"
#include "memory/linprob_hashmap.c"
#include "math/functions.c"
#include "math/matrix.c"
#include "math/find_eigenpairs.c"
#include "math/basis.c"
#include "methods/quadgk.c"
#include "methods/diis.c"
#include "methods/nlse_solver.c"
#include "methods/grosspitaevskii.c"
#include "methods/best_meanfield.c"
#include "methods/perturbation_theory.c"

/*
 * Logging definitions
 */

#define MAX_LOG_MSG_LEN 128

void sbmf_set_log_callback(sbmf_log_callback_func* func) {
	_state.log_callback = func;
}

static void sbmf_log(enum sbmf_log_level log_level, const char* fmt, va_list args) {
	if (!_state.log_callback)
		return;

	static char msg_buffer[MAX_LOG_MSG_LEN];

	vsnprintf(msg_buffer, MAX_LOG_MSG_LEN, fmt, args);
	_state.log_callback(log_level, msg_buffer);
}

void sbmf_log_info(const char* fmt, ...) {
	va_list args;
	va_start(args, fmt);
	sbmf_log(SBMF_LOG_LEVEL_INFO, fmt, args);
	va_end(args);
}

void sbmf_log_warning(const char* fmt, ...) {
	va_list args;
	va_start(args, fmt);
	sbmf_log(SBMF_LOG_LEVEL_WARNING, fmt, args);
	va_end(args);
}

void sbmf_log_error(const char* fmt, ...) {
	va_list args;
	va_start(args, fmt);
	sbmf_log(SBMF_LOG_LEVEL_ERROR, fmt, args);
	va_end(args);
}

void sbmf_log_panic(const char* fmt, ...) {
	va_list args;
	va_start(args, fmt);
	sbmf_log(SBMF_LOG_LEVEL_PANIC, fmt, args);
	va_end(args);
}
