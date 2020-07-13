#include <sbmf/common/common.h>

#define MAX_GAUSS_POINTS 7
typedef struct {
	f64 kronod_nodes[MAX_GAUSS_POINTS+1];
	f64 kronod_weights[MAX_GAUSS_POINTS+1];
	f64 gauss_weights[(MAX_GAUSS_POINTS+1)/2];
	u32 kronod_size;
	u32 gauss_size;
} gk_data;

extern gk_data gk7;

typedef f64 integrand(f64, void*);
typedef struct integration_settings {
	gk_data gk;

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

integration_result quadgk(integrand* f, f64 start, f64 end, integration_settings settings);
