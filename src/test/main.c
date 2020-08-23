#define SNOW_ENABLED
#include "snow.h"

#include <sbmf/methods/find_groundstate.h>
#include <sbmf/math/grid.h>
#include <sbmf/math/harmonic_oscillator.h>

#include <sbmf/math/find_eigenpairs.h>

#include <plot/plot.h>

static inline f64 gaussian(f64 x, f64 mu, f64 sigma) {
	return 1.0/(sigma*sqrt(M_2_PI)) * exp(-x*x/(2*sigma*sigma));
}
void f64_normalize(f64* data, u32 size) {
	f64 sum = 0.0;
	for (u32 i = 0; i < size; ++i) {
		f64 absval = fabs(data[i]);
		sum += absval*absval;
	}
	for (u32 i = 0; i < size; ++i) {
		data[i] *= (1.0/sqrt(sum));
	}
}
void c64_normalize(c64* data, u32 size) {
	f64 sum = 0.0;
	for (u32 i = 0; i < size; ++i) {
		f64 absval = cabs(data[i]);
		sum += absval*absval;
	}
	for (u32 i = 0; i < size; ++i) {
		data[i] *= (1.0/sqrt(sum));
	}
}

static inline f64 ho_integrand(f64 x, void* data) {
	u32* params = data;
	return
		ho_eigenfunction((i32[]){params[0]}, &x, 1) *
		ho_potential(&x,1,0)*
		ho_eigenfunction((i32[]){params[1]}, &x, 1);
}

static inline f64 ho_perturbed_potential(f64* x, i32 n, void* data) {
	return ho_potential(x,1,0) + gaussian(*x,0,0.2);
}

static inline f64 ho_perturbed_integrand(f64 x, void* data) {
	u32* params = data;
	return
		ho_eigenfunction((i32[]){params[0]}, &x, 1) *
		ho_perturbed_potential(&x,1,0) *
		ho_eigenfunction((i32[]){params[1]}, &x, 1);
}

//#define PLOT_SPARSE_1D_HO_FDM 		0
//#define PLOT_SPARSE_2D_HO_FDM 		0
//#define PLOT_SPARSE_1D_PB_FDM 		0
//#define PLOT_SPARSE_2D_PB_FDM 		0
//#define PLOT_SPARSE_1D_PB_HOBASIS 	0

#include "integration.c"
//#include "particle_in_a_box.c"
//#include "fdm_ho_comp.c"
//#include "fdm_solving.c"
//#include "ho_basis_func.c"
#include "find_groundstate.c"

snow_main();
