#include <sbmf/sbmf.h>

#include <stdio.h>

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

	/* We are only interested in the wave function for N = 4 particles and
	 * g0 = -1/3 interaction strength. As we will scaling in the mean-field limit
	 * the product g(N-1) will be constant and all wavefunctions will look the same.
	 */
	i64 Ns[] = {4};
	f64 g0 = -1.0/3.0;

	struct nlse_settings settings = {
        .spatial_pot_perturbation = perturbation,
		.max_iterations = 3000,
		.max_quadgk_iters = 500,
		.abs_error_tol = 1e-14,

        .num_basis_funcs = 64,
		.basis = ho_basis,
		.hamiltonian_mixing = 0.95,

		.zero_threshold = 1e-10,

		.gk=gk20,
    };

	const u32 component_count = 1;

	/* We are using sbmf_init()/sbmf_shutdown() in the loop to guarantee a full reset of all state/memory in case of
	 * for instance a memory leak. */
	for (u32 i = 0; i < sizeof(Ns)/sizeof(Ns[0]); ++i) {
		sbmf_init();
		i64 N = Ns[i];

		/* Solve the GP equations for both guesses */
		struct nlse_result gp_default_res = grosspitaevskii(settings, component_count, &N, default_guesses, &g0);
		struct nlse_result gp_random_res = grosspitaevskii(settings, component_count, &N, random_guesses, &g0);

		/* Output results to file */
		{
			char buf[256];
			snprintf(buf, 256, "out_W_N%ld", N);

			FILE* fd = fopen(buf, "a");
			fprintf(fd, "# x\tpotential\tdefault\trandom\n");

			for (f64 x = -3; x <= 3; x += 0.01) {
				/* Sampling one x at a time is inefficient, but it doesnt really matter here */
				f64 default_out = 0, random_out = 0;
				ho_sample(gp_default_res.coeff_count, &gp_default_res.coeff[0], 1, &default_out, &x);
				ho_sample(gp_random_res.coeff_count, &gp_random_res.coeff[0], 1, &random_out, &x);
				f64 pot = 0.5*x*x + 2*gaussian(x, 0, 0.2);
				fprintf(fd, "%lf\t%lf\t%lf\t%lf\n", x, pot,
						default_out*default_out,
						random_out*random_out);
			}

			fclose(fd);
		}

		sbmf_shutdown();
	}

	return 0;
}
