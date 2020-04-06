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
