#include <sbmf/common/common.h>

typedef f64 integrand(f64, void*);

typedef struct integration_settings {
	// Order of the Gauss-Kronod method used.
	// This setting is currently ignored.
	i32 order;

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
