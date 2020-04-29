#pragma once

// <complex.h> needs C99 compliance but makes fftw use the native complex type.
#include <complex.h> 
#include <float.h>
#include <stdint.h>
#include <math.h>

#define ARRLEN(arr) \
	(sizeof(arr)/sizeof(arr[0]))

// Definitions of common numerical types used.
typedef int32_t int_t;
typedef double real_t;
typedef double complex complex_t;

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

