#pragma once

#include <sbmf/types.h>
#include <math.h>

static inline f64 gaussian(f64 x, f64 mu, f64 sigma) {
	return 1.0/(sigma*sqrt(M_2_PI)) * exp(-(x-mu)*(x-mu)/(2*sigma*sigma));
}

static inline void f64_normalize(f64* data, u32 size) {
	f64 sum = 0.0;
	for (u32 i = 0; i < size; ++i) {
		f64 absval = fabs(data[i]);
		sum += absval*absval;
	}
	for (u32 i = 0; i < size; ++i) {
		data[i] *= (1.0/sqrt(sum));
	}
}

static inline void c64_normalize(c64* data, u32 size) {
	f64 sum = 0.0;
	for (u32 i = 0; i < size; ++i) {
		f64 absval = cabs(data[i]);
		sum += absval*absval;
	}
	for (u32 i = 0; i < size; ++i) {
		data[i] *= (1.0/sqrt(sum));
	}
}
