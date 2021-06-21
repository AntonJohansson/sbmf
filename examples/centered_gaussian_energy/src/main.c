#include <sbmf/sbmf.h>
#include <math.h>
#include <stdio.h>

#define PERTURBATION(x)

/* Define the perturbation to use relative to trap in the used basis.
 * In this case as we are using the HO basis, this perturbation will add a centered
 * gaussian to create a double well */
void perturbation(const u32 len, f64 out[static len],
                                f64 in_x[static len], const u32 component_count,
                                f64 in_u[static len*component_count],
								void* userdata) {
    for (u32 i = 0; i < len; ++i) {
        out[i] = 2*gaussian(in_x[i], 0, 0.2);
	}
}

/* Logging function for sbmf to use */
void log_callback(enum sbmf_log_level log_level, const char* msg) {
	printf("%s\n", msg);
}

int main() {
	sbmf_set_log_callback(log_callback);
	sbmf_set_thread_storage_size(256*1024*1024);

	/* Use random guesses for both component */
	struct nlse_guess random_guesses[] = {
		[0] = {
			.type = RANDOM_GUESS,
		},
		[1] = {
			.type = RANDOM_GUESS,
		},
	};
	/* NULL means use default, that is eigenstates of the basis */
	struct nlse_guess* default_guesses = NULL;

	struct nlse_settings settings = {
        .spatial_pot_perturbation = perturbation,
		.max_iterations = 15000,
		.max_quadgk_iters = 500,
		.abs_error_tol = 1e-14,

		.num_basis_funcs = 64,
		.basis = ho_basis,
		.hamiltonian_mixing = 0.75,

		.zero_threshold = 1e-10,

		.gk=gk20,
    };

	const u32 component_count = 1;

	/* Particle counts to solve for */
	i64 Ns[] = {4, 6, 8, 12, 16, 20, 24, 28, 32};

	/* Interaction strenghts lambda = g(N-1) to use */
	f64 gs[] = {-1.0/(4.0-1.0)};

	/* We are using sbmf_init()/sbmf_shutdown() in the loop to guarantee a full reset of all state/memory in case of
	 * for instance a memory leak. */
	for (u32 i = 0; i < sizeof(Ns)/sizeof(Ns[0]); ++i) {

		i64 N = Ns[i];
		for (u32 j = 0; j < sizeof(gs)/sizeof(gs[0]); ++j) {
			f64 g0 = gs[j];

			FILE* fd = fopen("out", "a");
			fprintf(fd, "%ld\t", N);

			sbmf_init();
			{
				/* Solve the GP equation using the default guess */
				struct nlse_result gp_default_res = grosspitaevskii(settings, component_count, &N, default_guesses, &g0);
				if (!gp_default_res.converged)
					break;

				/* Compute GP energy */
				f64 Egp_default = grosspitaevskii_energy(settings, gp_default_res.coeff_count, component_count, gp_default_res.coeff, &N, &g0);

				/* Compute Best mean field fractional occupation */
				f64 bmf_default = bestmf_find_fractional_occupation(settings, N, g0, default_guesses);

				/* Apply peturbation theory using cuda */
				struct pt_result rs_default_ptres = rspt_1comp_cuda(&settings, gp_default_res, 0, g0, N);
				struct pt_result en_default_ptres = enpt_1comp_cuda(&settings, gp_default_res, 0, g0, N);

				/* Save results */
				fprintf(fd, "%.10f\t%.10lf\t", Egp_default, bmf_default);
				fprintf(fd, "%.10lf\t%.10lf\t%.10lf\t%.10lf\t",
						rs_default_ptres.E0 + rs_default_ptres.E1 + rs_default_ptres.E2,
						rs_default_ptres.E0 + rs_default_ptres.E1 + rs_default_ptres.E2 + rs_default_ptres.E3,
						en_default_ptres.E0 + en_default_ptres.E1 + en_default_ptres.E2,
						en_default_ptres.E0 + en_default_ptres.E1 + en_default_ptres.E2 + en_default_ptres.E3
						);
			}
			sbmf_shutdown();

			sbmf_init();
			{
				/* Solve the GP equation using the random guess */
				struct nlse_result gp_random_res = grosspitaevskii(settings, component_count, &N, random_guesses, &g0);
				if (!gp_random_res.converged)
					break;

				/* Compute GP energy */
				f64 Egp_random = grosspitaevskii_energy(settings, gp_random_res.coeff_count, component_count, gp_random_res.coeff, &N, &g0);

				/* Compute Best mean field fractional occupation */
				f64 bmf_random = bestmf_find_fractional_occupation(settings, N, g0, random_guesses);

				/* Apply peturbation theory using cuda */
				struct pt_result rs_random_ptres = rspt_1comp_cuda(&settings, gp_random_res, 0, g0, N);
				struct pt_result en_random_ptres = enpt_1comp_cuda(&settings, gp_random_res, 0, g0, N);

				/* Save reults */
				fprintf(fd, "%.10f\t%.10lf\t", Egp_random, bmf_random);
				fprintf(fd, "%.10lf\t%.10lf\t%.10lf\t%.10lf\n",
						rs_random_ptres.E0 + rs_random_ptres.E1 + rs_random_ptres.E2,
						rs_random_ptres.E0 + rs_random_ptres.E1 + rs_random_ptres.E2 + rs_random_ptres.E3,
						en_random_ptres.E0 + en_random_ptres.E1 + en_random_ptres.E2,
						en_random_ptres.E0 + en_random_ptres.E1 + en_random_ptres.E2 + en_random_ptres.E3
						);
			}
			sbmf_shutdown();

			fclose(fd);
		}
	}

	return 0;
}
