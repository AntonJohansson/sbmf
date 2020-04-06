#pragma once

#include <sbmf/common/common.h>
#include <sbmf/common/grid.h>

typedef struct {
	grid g;

	//int_t resolution;

	int_t max_iterations;
	real_t error_tol;

	//real_t length_x;
	//real_t length_y;

	int_t measure_every;

	// Specific for aitem
	real_t c;
	real_t A;
	real_t dt;
} gss_settings;

typedef struct {
	int_t count;
	int_t current;

	real_t* error;
	real_t* iteration_time;
	complex_t* wavefunction;
} gss_debug;

typedef struct {
	gss_settings settings;
	gss_debug debug;

	//int_t rows;
	//int_t cols;

	//real_t* X;
	//real_t* Y;
	complex_t* wavefunction;

	real_t error;
	int_t iterations;
} gss_result;

typedef real_t gss_potential_func(real_t* v, int_t n, complex_t u);
typedef complex_t  gss_guess_func(real_t* v, int_t n);

extern void gss_free_result(gss_result res);

extern gss_result aitem_execute(gss_settings settings,
																gss_potential_func* potential, gss_guess_func* guess);

extern gss_result niter_execute(gss_settings settings,
																gss_potential_func* potential, gss_guess_func* guess);
extern gss_result item_execute(gss_settings settings,
																gss_potential_func* potential, gss_guess_func* guess);
