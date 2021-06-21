#include <sbmf/sbmf.h>
#include <stdio.h>

/* Logging function for sbmf to use */
void log_callback(enum sbmf_log_level log_level, const char* msg) {
	printf("%s\n", msg);
}

int main() {
	sbmf_set_log_callback(log_callback);

	/* Use default guesses for both components */
	struct nlse_guess* default_guesses = NULL;

	/* Values of lambda = g(N-1) to try */
	f64 l0s[] = {-0.5, -1.0, 1.0, 0.5};
	/* Particle counts to try */
	i64 Ns[] = {
		4,
		6,
		10,
		20,
		30,
		40,
		60,
		80,
		100,
		120,
		140,
		160
	};

	struct nlse_settings settings = {
		.max_iterations = 1e5,
		.max_quadgk_iters = 1000,
		.abs_error_tol = 1e-14,

		.num_basis_funcs = 64,
		.basis = ho_basis,

		.zero_threshold = 1e-10,
		.gk=gk20
    };

	const u32 component_count = 1;

	for (u32 j = 0; j < sizeof(l0s)/sizeof(l0s[0]); ++j) {
		for (u32 k = 0; k < sizeof(Ns)/sizeof(Ns[0]); ++k) {
			sbmf_init();

			/* We compute the interaction strength such that lambda = g(N-1) is constant at
			 * the desired value when we increase N */
			i64 N = Ns[k];
			f64 g0 = l0s[j]/((f64)N-1.0);

			/* Solve GP equation and compute energy */
			struct nlse_result res = grosspitaevskii(settings, component_count, &N, default_guesses, &g0);
			f64 Egp = grosspitaevskii_energy(settings, res.coeff_count, component_count, res.coeff, &N, &g0);

			/* Compute petrubative energy shifts */
			struct pt_result rs_ptres = rspt_1comp_cuda(&settings, res, 0, g0, N);
			struct pt_result en_ptres = enpt_1comp_cuda(&settings, res, 0, g0, N);

			/* Save data to file */
			{
				char buf[256];
				snprintf(buf, 256, "out_l%.2lf", l0s[j]);

				FILE* fd = fopen(buf, "a");
				fprintf(fd, "%ld\t%.10f\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\n",
						N,
						Egp,
						rs_ptres.E0+rs_ptres.E1,
						rs_ptres.E0+rs_ptres.E1+rs_ptres.E2,
						rs_ptres.E0+rs_ptres.E1+rs_ptres.E2+rs_ptres.E3,
						en_ptres.E0+en_ptres.E1,
						en_ptres.E0+en_ptres.E1+en_ptres.E2,
						en_ptres.E0+en_ptres.E1+en_ptres.E2+en_ptres.E3
					   );
				fclose(fd);
			}

			sbmf_shutdown();
		}

	}

}
