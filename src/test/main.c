#include <sbmf/groundstate_solver/groundstate_solver.h>
#include <sbmf/common/profile.h>
#include <sbmf/common/fdm.h>
#include <sbmf/common/grid.h>

#include <lapacke.h>

#include <stdio.h>
#include <assert.h>

real_t potential(real_t x, real_t y, complex_t u) {
	real_t cosx = cos(x);
	real_t cosy = cos(y);
	//return x*x + y*y + cabs(u)*cabs(u);
	//return 3*(cosx*cosx + cosy*cosy)*u + u*u*u;
	return 3*(cosx*cosx + cosy*cosy) + cabs(u)*cabs(u);
	//return x*x + y*y;
	//return x*x + y*y;
	//return 3*(x*x)*u;
	//return (fabs(x) < 15 && fabs(y) < 15) ? 0 : 10*u;
}

complex_t initial_guess(real_t x, real_t y) {
	//return (1.0/(10.0*M_PI*10.0*M_PI + 10.0*M_PI*10.0*M_PI))*exp(-(x*x + y*y));
	return 1.0/(10.0*M_PI*10.0*M_PI);
}

void aitem_test() {
	gss_settings settings = {
		.resolution = 256,
		.max_iterations = 500,
		.error_tol = 1e-10,
		.c = 0.7,
		.A = 1.0384,
		.dt = 1,
		.length_x = 10.0*M_PI,
		.length_y = 10.0*M_PI
	};

	PROFILE_BEGIN("aitem");

	gss_result res = aitem_execute(settings,
																 potential,
																 initial_guess);
	PROFILE_END("aitem");

	FILE* fdwf = fopen("gss_aitem_out", "w");
	for (int_t i = 0; i < res.rows*res.cols; ++i) {
		real_t cr = creal(res.wavefunction[i]);
		real_t ci = cimag(res.wavefunction[i]);
		fprintf(fdwf, "%lf\n", cabs(res.wavefunction[i]));
		//fprintf(fdwf, "%lf%s%lfi\n",
		//				cr,
		//				(ci < 0.0) ? "-" : "+",
		//				fabs(ci));
	}


	gss_free_result(res);
}

void item_test() {
	const int_t N = 256;
	grid g = generate_grid(2, (dimension_info[]){
			{-10.0*M_PI,10.0*M_PI,N},
			{-10.0*M_PI,10.0*M_PI,N},
			});

	{
		FILE* fd = fopen("debug_grid", "w");
		for (int_t i = 0; i < g.total_pointcount; ++i) {
			for (int_t j = 0; j < g.dimensions; ++j) {
				fprintf(fd, "%lf", g.points[g.dimensions*i + j]);
				if (j < g.dimensions-1)
					fprintf(fd, "\t");
				else
					fprintf(fd, "\n");
			}
		}
		fclose(fd);
	}


	gss_settings settings = {
		.g = g,
		.resolution = N,
		.max_iterations = 2000,
		.error_tol = 1e-10,
		.dt = 0.01,
		.length_x = 10.0*M_PI,
		.length_y = 10.0*M_PI,

		.measure_every = 0,
	};

	printf("Finding groundstate\n");
	PROFILE_BEGIN("item");

	gss_result res = item_execute(settings,
																potential,
																initial_guess);
	PROFILE_END("item");

	printf("Generating FD matrix\n");
	PROFILE_BEGIN("gen fd");
	bandmat fdm = generate_fd_matrix(settings.resolution, 4, 2, (real_t[]){
			settings.length_x/settings.resolution,
			settings.length_y/settings.resolution});

	for (int_t i = 0; i < fdm.size*fdm.bandcount; ++i) {
		fdm.bands[i] = -fdm.bands[i];
	}

	for (int_t i = 0; i < fdm.size; ++i) {
		int_t idx = fdm.size*(fdm.bandcount-1) + i;
		fdm.bands[idx] += potential(res.X[i], res.Y[i], res.wavefunction[i]);
	}	
	PROFILE_END("gen fd");

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

	{
		FILE* fdwf = fopen("gss_item_wf", "w");
		FILE* fdpt = fopen("gss_item_pt", "w");
		FILE* fdev = fopen("gss_item_ev", "w");
		FILE* fdvv = fopen("gss_item_vv", "w");
		for (int_t i = 0; i < res.rows*res.cols; ++i) {
			fprintf(fdwf, "%lf\n", cabs(res.wavefunction[i]));
			fprintf(fdpt, "%lf\n", potential(res.X[i], res.Y[i], res.wavefunction[i]));
			//fprintf(fdev, "%lf\n", outD[i]);
			//for (int_t j = 0; j < fdm.size; ++j) {
			//	fprintf(fdvv, "%lf", outQ[i*fdm.size + j]);
			//	if (j < fdm.size-1)
			//		fprintf(fdvv, "\t");
			//	else
			//		fprintf(fdvv, "\n");
			//}
		}
		fclose(fdwf);
		fclose(fdpt);
		fclose(fdev);
		fclose(fdvv);
	}
	//{
	//	FILE* fd = fopen("gss_item_error", "w");
	//	for (int_t i = 0; i < res.debug.current; ++i) {
	//		fprintf(fd, "%lf\n", res.debug.error[i]);
	//	}
	//	fclose(fd);
	//}
	
	gss_free_result(res);
	free_grid(g);
}

/*
void niter_test() {
	niter_settings settings = {
		.resolution = 20,
		.max_iterations = 500,
		.error_tol = 1e-10,
		.length_x = 10.0*M_PI,
		.length_y = 10.0*M_PI,
	};

	niter_result res = niter_execute(settings, potential, initial_guess);
	niter_free_result(res);
}
*/
int main(int argc, char** argv) {
	//aitem_test();
	item_test();

	printf("Timing results:\n");
	for (size_t i = 0; i < profile_used_data; ++i) {
		profile_data_t data = profile_data[i];

		long avg = 0;
		for (size_t j = 0; j < data.delta_num_samples; ++j) {
			avg += 1000000000*data.deltas[j].tv_sec + data.deltas[j].tv_nsec;
		}
		avg /= data.delta_num_samples;

		printf("%10s -- %ld (ms)\n", data.name, avg/1000);
	}

	return 0;
}
