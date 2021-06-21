#include <sbmf/sbmf.h>
#include <stdio.h>

#define USE_TF_GUESS 0
#define USE_GAUSSIAN_GUESS 0
#define USE_RANDOM_GUESS 0

/* Logging function for sbmf to use */
void log_callback(enum sbmf_log_level log_level, const char* msg) {
	printf("%s\n", msg);
}

void gaussian0(f64* out, f64* in, u32 len, void* data) {
	for (u32 i = 0; i < len; ++i) {
		//out[i] = (1 - 0.5 * in[i]*in[i]);
		out[i] = gaussian(in[i] + 1.0, 0.0, 0.1);// + gaussian(in[i] - 1.0, 0.0, 0.2);
	}
}
void gaussian1(f64* out, f64* in, u32 len, void* data) {
	for (u32 i = 0; i < len; ++i)
		out[i] = gaussian(in[i] - 1.0, 0.0, 0.1);
}



int main() {
	sbmf_set_log_callback(log_callback);

#if USE_GAUSSIAN_GUESS
	struct nlse_guess guesses[] = {
		[0] = {
			.type = SPATIAL_GUESS,
			.data.spatial_guess = gaussian0,
		},
		[1] = {
			.type = SPATIAL_GUESS,
			.data.spatial_guess = gaussian1,
		},
	};
#elif USE_RANDOM_GUESS
	struct nlse_guess guesses[] = {
		[0] = {
			.type = RANDOM_GUESS,
		},
		[1] = {
			.type = RANDOM_GUESS,
		},
	};
#else
	struct nlse_guess* guesses = NULL;
#endif

	struct nlse_settings settings = {
		.max_iterations = 3000,
		.max_quadgk_iters = 1000,
		.abs_error_tol = 1e-14,

		.num_basis_funcs = 48,
		.basis = ho_basis,

		.zero_threshold = 1e-10,
		.hamiltonian_mixing = 0.97,

		.gk=gk20
    };

	/* This time we use two components */
	const u32 component_count = 2;

	/* Here are the particle counts per component to try */
	i64 Ns[] = {4,6, 8,10,15,20,25,30,35,40,45,50,55,60,65,70};
	/* These are the values of the oscillator strength of the
	 * harmonic potential. We use a weak trapping potential here. */
	f64 Os[] = {0.0050};
	/* These are the ratios between the inter- and intra-component
	 * interactions strengths we use. The intra-component one
	 * will be kept fix. */
	f64 gAB_factors[] = {-0.90,-0.95,-0.99};

	f64 lambda = 0.5;
	for (u32 k = 0; k < sizeof(gAB_factors)/sizeof(gAB_factors[0]); ++k) {
		f64 gAB_factor = gAB_factors[k];

		for (u32 j = 0; j < sizeof(Os)/sizeof(Os[0]); ++j) {
			/* Set the value of OMEGA */
			OMEGA = Os[j];

			/* Write header to file */
			char buf[128];
			snprintf(buf, 128, "out_%lf_gab_%lf", Os[j], gAB_factor);
			{
				FILE* fd = fopen(buf, "a");
				fprintf(fd, "# N\tE\tRS2\tRS3\tEN2\tEN3\n");
				fclose(fd);
			}

			for (u32 i = 0; i < sizeof(Ns)/sizeof(Ns[0]); ++i) {
				/* Particle counts per component */
				i64 N = Ns[i];
				i64 occupations[] = {N,N};

				/* Compute interaction strengths */
				const i64 max_N = 100;
				const f64 g = lambda/((f64)max_N - 1.0);
				f64 g0[] = {
					g,            gAB_factor*g,
					gAB_factor*g, g
				};

				sbmf_init();

				/* Solve the GP equation */
				struct nlse_result res = grosspitaevskii(settings, component_count, occupations, guesses, g0);
				if (!res.converged)
					break;

				/* Sample GP wave function at x = 0 so we can estimate the number
				 * density of the condensate later on */
				{
					f64 n = 0;
					f64 x_in = 0;
					ho_sample(res.coeff_count, res.coeff, 1, &n, &x_in);

					FILE* fd = fopen("out_n", "a");
					fprintf(fd, "%ld\t%.10e\n", N, n);
					fclose(fd);
				}

				/* Compute GP energy */
				f64 Efull = grosspitaevskii_energy(settings, res.coeff_count, component_count, res.coeff, occupations, g0);
				printf("\nfull energy: %lf\n", Efull);

				/* Compute perturbative energy shifts */
				struct pt_result rspt = rspt_2comp_cuda(&settings, res, 0, 1, g0[0], g0[1], occupations[0], occupations[1]);
				struct pt_result enpt = enpt_2comp_cuda(&settings, res, 0, 1, g0[0], g0[1], occupations[0], occupations[1]);

				/* Save results to file */
				{
					FILE* fd = fopen(buf, "a");
					fprintf(fd, "%ld\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\n",
							N,
							Efull,
							rspt.E0+rspt.E1+rspt.E2,
							rspt.E0+rspt.E1+rspt.E2+rspt.E3,
							enpt.E0+enpt.E1+enpt.E2,
							enpt.E0+enpt.E1+enpt.E2+enpt.E3
							);
					fclose(fd);
				}

				sbmf_shutdown();
			}
		}
	}



}
