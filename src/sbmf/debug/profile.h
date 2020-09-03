#pragma once

#define PROFILE_ENABLE 1

/* Implementations of the functions
 * wrapped in macros below.
 * Never use these directly.
 */
void profile_begin_impl(char const name[]);
void profile_end_impl(char const name[]);
void profile_print_results_impl();
void profile_clear_impl();

#if PROFILE_ENABLE
	#define PROFILE_BEGIN(name)\
		profile_begin_impl(name)
	#define PROFILE_END(name)\
		profile_end_impl(name)

	/* Returns the result of func_call via the
	 * comma operator. Expressions such as:
	 *
	 * 	int x = PROFILE_FUNC(compute_x())
	 *
	 * should therefore be possible.
	 */
	#define PROFILE_FUNC(func_call)\
		(PROFILE_BEGIN( #func_call ), func_call); PROFILE_END( #func_call )

	#define profile_print_results() \
		profile_print_results_impl()
	#define profile_clear() \
		profile_clear_impl()
#else
	#define PROFILE_BEGIN(name)
	#define PROFILE_END(name)
	#define PROFILE_FUNC(func_call)\
		func_call
	#define profile_print_results()
	#define profile_clear()
#endif
