#include <sbmf/groundstate_solver/groundstate_solver.h>
#include <sbmf/common/profile.h>
#include <sbmf/common/fdm.h>
#include <sbmf/common/grid.h>

#include <lapacke.h>

#include <stdio.h>
#include <assert.h>
#include <unistd.h>

#include "plotting/plt.h"

static PlotState* state;

real_t potential(real_t* v, int_t n, complex_t u) {
	//real_t temp = 0.0f;
	//for (int_t i = 0; i < n; ++i)
	//	temp += cos(v[i])*cos(v[i]);
	//return 3*temp + cabs(u)*cabs(u);

	real_t temp = 0.0f;
	for (int_t i = 0; i < n; ++i) {
		temp += v[i]*v[i];
	}
	return 0.5*temp;
}

complex_t initial_guess(real_t* v, int_t n) {
	//return (1.0/(10.0*m_pi*10.0*m_pi + 10.0*m_pi*10.0*m_pi))*exp(-(x*x + y*y));
	return 1.0/(10.0*M_PI*10.0*M_PI);
}

void dbgcallback(grid g, complex_t* wf) {
	float x[g.total_pointcount];
	float y[g.total_pointcount];
	float v[g.total_pointcount];
	for (int_t i = 0; i < g.total_pointcount; ++i) {
		x[i] = (float)g.points[i];
		y[i] = (float)cabs(wf[i]);
		v[i] = (float)potential(&g.points[i], 1, y[i]);
	}

	plt_1d(state, 0, x, v, g.total_pointcount);
	plt_1d(state, 1, x, y, g.total_pointcount);

	plt_update(state);
}

int main(int argc, char** argv) {
	const int_t N = 128;

	grid g = generate_grid(1, 
			(real_t[]){-5*M_PI, -5*M_PI},
			(real_t[]){ 5*M_PI,  5*M_PI},
			 (int_t[]){ N,  	 	 N     }
			);

	gss_settings settings = {
		.g = g,
		.max_iterations = 2000,
		.error_tol = 1e-10,
		.dt = 0.01,

		.measure_every = 1,
		.dbgcallback = dbgcallback,
	};

	printf("Finding groundstate\n");

	state = plt_init();
	PROFILE_BEGIN("item");

	gss_result res = item_execute(settings,
																potential,
																initial_guess);
	PROFILE_END("item");
	plt_update_until_closed(state);
	plt_shutdown(state);

	//printf("Generating FD matrix\n");
	//PROFILE_BEGIN("gen fd");
	//bandmat fdm = generate_fd_matrix(settings.resolution, 4, 2, (real_t[]){
	//		settings.length_x/settings.resolution,
	//		settings.length_y/settings.resolution});

	//for (int_t i = 0; i < fdm.size*fdm.bandcount; ++i) {
	//	fdm.bands[i] = -fdm.bands[i];
	//}

	//for (int_t i = 0; i < fdm.size; ++i) {
	//	int_t idx = fdm.size*(fdm.bandcount-1) + i;
	//	//fdm.bands[idx] += potential(res.X[i], res.Y[i], res.wavefunction[i]);
	//}	
	//PROFILE_END("gen fd");

	printf("Trying to solve the eigenvalue problem\n");
	//PROFILE_BEGIN("EV prob");
	//{
	//	real_t* outD = malloc(fdm.size*sizeof(real_t));
	//	real_t* outE = malloc((fdm.size-1) * sizeof(real_t));
	//	complex_t* outQ = malloc(fdm.order*sizeof(complex_t));

	//	{
	//		int res = 0;

	//		// Reduce to tridiag-form
	//		res = LAPACKE_zhbtrd(LAPACK_ROW_MAJOR, 'V', 'U', 
	//				fdm.size,
	//				fdm.bandcount-1, 
	//				fdm.bands, 
	//				fdm.size, 
	//				outD, outE, outQ, 
	//				fdm.size);

	//		assert(res == 0);

	//		// Solve eigenvalue problem
	//		res = LAPACKE_zsteqr(LAPACK_ROW_MAJOR, 'V', 
	//				fdm.size, 
	//				outD, outE, outQ, 
	//				fdm.size);

	//		assert(res == 0);
	//	}

	//}
	//PROFILE_END("EV prob");

	printf(
			"item results:\n"
			"\titerations: %d/%d\n"
			"\terror: %e\n"
			"\terror tol: %e\n",
		res.iterations,
		res.settings.max_iterations,
		res.error,
		res.settings.error_tol
	);

	//{
	//	FILE* fdwf = fopen("gss_item_wf", "w");
	//	FILE* fdpt = fopen("gss_item_pt", "w");
	//	FILE* fdev = fopen("gss_item_ev", "w");
	//	FILE* fdvv = fopen("gss_item_vv", "w");
	//	for (int_t i = 0; i < settings.g.total_pointcount; ++i) {
	//		fprintf(fdwf, "%lf\n", cabs(res.wavefunction[i]));
	//		//fprintf(fdpt, "%lf\n", potential(res.X[i], res.Y[i], res.wavefunction[i]));
	//		//fprintf(fdev, "%lf\n", outD[i]);
	//		//for (int_t j = 0; j < fdm.size; ++j) {
	//		//	fprintf(fdvv, "%lf", outQ[i*fdm.size + j]);
	//		//	if (j < fdm.size-1)
	//		//		fprintf(fdvv, "\t");
	//		//	else
	//		//		fprintf(fdvv, "\n");
	//		//}
	//	}
	//	FOREACH_GRIDPOINT(settings.g, i) {
	//		fprintf(fdpt, "%lf\n", potential(&settings.g.points[i], settings.g.dimensions, res.wavefunction[i]));
	//	}
	//	fclose(fdwf);
	//	fclose(fdpt);
	//	fclose(fdev);
	//	fclose(fdvv);
	//}
	//{
	//	FILE* fd = fopen("gss_item_error", "w");
	//	for (int_t i = 0; i < res.debug.current; ++i) {
	//		fprintf(fd, "%lf\n", res.debug.error[i]);
	//	}
	//	fclose(fd);
	//}
	
	gss_free_result(res);
	free_grid(g);

	profile_print_results();
	return 0;
}
