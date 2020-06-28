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
typedef float 		f32;
typedef double 		f64;
typedef float complex 	c32;
typedef double complex 	c64;

#include <stdio.h>
static inline bool float_compare(f64 a, f64 b, f64 epsilon) {
	return (fabs(a-b) <= epsilon);
}

struct stack_allocator;

typedef struct sbmf_state {
	struct stack_allocator* main_stack;
} sbmf_state;

extern sbmf_state _sbmf;

void sbmf_init();
void sbmf_shutdown();
