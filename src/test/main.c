#include <sbmf/groundstate_solver/groundstate_solver.h>
#include <sbmf/common/profile.h>
#include <sbmf/common/fdm.h>
#include <sbmf/common/grid.h>
#include <sbmf/common/eigenproblem.h>

#include <lapacke.h>

#include <stdio.h>
#include <assert.h>
#include <unistd.h>

#include "plotting/plt.h"

real_t potential(real_t* v, int_t n, complex_t u) {
	//real_t temp = 0.0f;
	//for (int_t i = 0; i < n; ++i)
	//	temp += cos(v[i])*cos(v[i]);
	//return 5*temp + cabs(u)*cabs(u);

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

static inline void test_eigenvalue_solving() {
	printf("Testing eigenvalue solving for band matrices\n");

	complex_t bands[] = {
		0,0,0,1,2,3,4,5,6,7,
		0,9,8,7,6,5,4,3,2,1,
		1,1,1,1,1,1,1,1,1,1
	};

	bandmat bm = {
		.order = 100,
		.size = 10,
		.bandcount = 3,
		.bands = bands,
	};

	real_t eigvals[bm.size];
	complex_t eigvecs[bm.order];
	eigvp_bandmat(eigvals, eigvecs, bm);

	// produces by matlab, warning COLUMN MAJOR
	complex_t eigvecs_answer[] = {
		-0.4489,    0.0085,    0.2281,   -0.2686,    0.3344,    0.4139,    0.2168,   -0.4228,   -0.3533,    0.2015,
		0.6377 ,  -0.0089 ,  -0.2163 ,   0.1788 ,  -0.1397 ,   0.0469 ,   0.0917 ,  -0.3495 ,  -0.4866 ,   0.3542 ,
		-0.5494,   -0.0178,   -0.0802,    0.2295,   -0.3348,   -0.4370,   -0.1732,    0.1008,   -0.3448,    0.4203,
		0.2818 ,   0.1498 ,   0.4358 ,  -0.4890 ,   0.1939 ,  -0.1808 ,  -0.2166 ,   0.3980 ,  -0.0938 ,   0.4297 ,
		-0.0251,   -0.4047,   -0.3178,    0.3091,    0.5092,    0.2221,    0.0615,    0.3800,    0.1375,    0.4042,
		-0.0817,    0.3838,   -0.3454,   -0.2371,   -0.4336,    0.4985,   -0.0242,    0.1086,    0.3041,    0.3595,
		0.0342 ,   0.2572 ,   0.4964 ,   0.4522 ,  -0.0600 ,  -0.0767 ,   0.5003 ,  -0.0762 ,   0.3593 ,   0.2951 ,
		0.0375 ,  -0.6138 ,   0.2489 ,  -0.0935 ,  -0.2515 ,   0.0495 ,  -0.3502 ,  -0.3963 ,   0.3852 ,   0.2395 ,
		-0.0204,   -0.0829,   -0.3886,   -0.4527,    0.1130,   -0.5182,    0.4676,   -0.2221,    0.2552,    0.1495,
		-0.0189,    0.4651,   -0.1586,    0.1849,    0.4382,   -0.1679,   -0.5212,   -0.4027,    0.2381,    0.1154
	};

	real_t eigvals_answer[] = {
		-11.7856,
		-8.4171 ,
		-7.5371 ,
		-4.9893 ,
		-2.7593 ,
		2.27    ,
		4.866   ,
		8.448   ,
		13.3963 ,
		16.8239 
	};

	// check answer
	for (int i = 0; i < bm.size; ++i) {
		if (fabs(eigvals[i] - eigvals_answer[i]) > 0.01) {
			fprintf(stderr, "Eigenvalue calculation %d failed: got %lf; expected %lf\n", i, eigvals[i], eigvals_answer[i]);
		}
	}

	(void)eigvecs_answer;
	// Check eigenvectors
	for (int i = 0; i < bm.size; ++i) {
		for (int j = 0; j < bm.size; ++j) {
			complex_t got = eigvecs[j*bm.size + i]; // NOTE: output of lapacke is apparently COLUMN MAJOR
			complex_t expected = eigvecs_answer[j*bm.size + i];
			if (cabs(got)-cabs(expected) > 0.01) {
				fprintf(stderr, "Eigenvector calculation failed: got %lf + %lf i; expected: %lf + %lf i\n", creal(got), cimag(got), creal(expected), cimag(expected));
			}
		}
	}
}

int main(int argc, char** argv) {
	test_eigenvalue_solving();

	PlotState* state = plt_init();

	const int_t N = 128;

	grid g = generate_grid(1, 
			(real_t[]){-5, -5},
			(real_t[]){ 5,  5},
			 (int_t[]){ N,  	 	 N     }
			);

	gss_settings settings = {
		.g = g,
		.max_iterations = 2000,
		.error_tol = 1e-10,
		.dt = 0.01,

		.measure_every = 0,
		//.dbgcallback = dbgcallback,
	};

	printf("Finding groundstate\n");

	PROFILE_BEGIN("item");

	gss_result res = item_execute(settings,
																potential,
																initial_guess);
	PROFILE_END("item");

	float x[g.total_pointcount];
	float y[g.total_pointcount];
	float v[g.total_pointcount];
	for (int_t i = 0; i < g.total_pointcount; ++i) {
		x[i] = (float)g.points[i];
		y[i] = (float)creal(res.wavefunction[i]);
		v[i] = (float)potential(&g.points[i], 1, res.wavefunction[i]);
	}

	plt_1d(state, x, v, g.total_pointcount); 
	plt_1d(state, x, y, g.total_pointcount); 








	printf("Generating FD matrix\n");
	PROFILE_BEGIN("gen. fdm");
	// NOTE: g.pointcounts[0] forces it to be square!
	bandmat fdm = generate_fd_matrix(g.pointcounts[0], pow(2,g.dimensions), g.dimensions, g.deltas);

	for (int_t i = 0; i < fdm.size*fdm.bandcount; ++i) {
		fdm.bands[i] = -0.5*fdm.bands[i];
	}

	for (int_t i = 0; i < fdm.size; ++i) {
		int_t idx = fdm.size*(fdm.bandcount-1) + i;
		fdm.bands[idx] += (complex_t) potential(&g.points[i], 1, res.wavefunction[i]);
	}	
	PROFILE_END("gen. fdm");




	printf("Trying to solve the eigenvalue problem\n");
	PROFILE_BEGIN("e.v. prob.");
	{
		real_t* eigenvalues = malloc(fdm.size * sizeof(real_t));
		complex_t* eigenvectors = malloc(fdm.order * sizeof(complex_t));

		eigvp_bandmat(eigenvalues, eigenvectors, fdm);

		float x[fdm.size];
		float u[fdm.size];
		for (int i = 1; i < 6; ++i) {
			for (int j = 0; j < fdm.size; ++j) {
				x[j] = (float)g.points[j];
				u[j] = (float)eigenvalues[i] +  5*(float)creal(eigenvectors[j*fdm.size + i]);
			}
			plt_1d(state, x,u, fdm.size);
		}
	}
	PROFILE_END("e.v. prob.");
















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
	
	gss_free_result(res);
	free_grid(g);

	profile_print_results();

	plt_wait_on_join(state);
	plt_shutdown(state);

	return 0;
}
