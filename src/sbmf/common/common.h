#pragma once

#include <complex.h>
#include <float.h>
#include <stdint.h>
#include <math.h>
#include <stdbool.h>

#define ARRLEN(arr) \
	(sizeof(arr)/sizeof(arr[0]))

// Definitions of common numerical types used.
typedef int8_t 		i8;
typedef uint8_t 	u8;
typedef int16_t 	i16;
typedef uint16_t 	u16;
typedef int32_t 	i32;
typedef uint32_t 	u32;
typedef int64_t 	i64;
typedef uint64_t 	u64;
typedef __int128_t  i128;
typedef __uint128_t u128;
typedef float 		f32;
typedef double 		f64;
typedef long double f128;
typedef float complex 		c32;
typedef double complex 		c64;
typedef long double complex c128;

static inline bool f64_compare(f64 a, f64 b, f64 epsilon) {
	return (fabs(a-b) <= epsilon);
}

static inline bool f64_is_valid(f64 f) {
	return !isinf(f) && !isnan(f);
}


struct stack_allocator;

typedef struct sbmf_state {
	struct stack_allocator* main_stack;
	bool initialized;
} sbmf_state;

extern sbmf_state _sbmf;

void sbmf_init();
void sbmf_shutdown();
