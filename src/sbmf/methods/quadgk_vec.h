#pragma once

#include <sbmf/types.h>

#define MAX_GAUSS_POINTS 20
struct gk_data {
	f64 kronod_nodes[MAX_GAUSS_POINTS+1];
	f64 kronod_weights[MAX_GAUSS_POINTS+1];
	f64 gauss_weights[(MAX_GAUSS_POINTS+1)/2];
	u32 kronod_size;
	u32 gauss_size;
};

extern struct gk_data gk7;
extern struct gk_data gk10;
extern struct gk_data gk15;
extern struct gk_data gk20;

typedef f64 integrand(f64,void*);
typedef struct integration_settings {
	struct gk_data gk;

	// Aboslute and relative error tolarences.
	f64 abs_error_tol;
	f64 rel_error_tol;

	// The maximum allowed function evaluations of
	// the supplied integrand.
	i32 max_evals;

	void* userdata;
} integration_settings;

typedef struct integration_result {
	f64 integral;
	f64 error;
	i32 performed_evals;
	bool converged;
} integration_result;

typedef void integrand_vec(f64*,f64*,u32,void*);

integration_result quadgk_vec(integrand_vec* f, f64 start, f64 end, integration_settings settings);
