#include <sbmf/sbmf.h>
#include <stdio.h>


/* This will be called every SCF iteration */
void debug_callback(struct nlse_settings settings, struct nlse_result res) {
	/* Save current iterations eigenvalue stored in res.energy[0] to file */
	{
		FILE* fd = fopen("debug_out", "a");
		fprintf(fd, "%u\t%lf\n", res.iterations, res.energy[0]);
		fclose(fd);
	}
}

/* Logging function for sbmf to use */
void log_callback(enum sbmf_log_level log_level, const char* msg) {
	printf("%s\n", msg);
}

int main() {
	sbmf_set_log_callback(log_callback);
	sbmf_init();

	/* NULL means just use default guess for SCF calculations */
	struct nlse_guess* guesses = NULL;

	const u32 component_count = 1;

	/* Particle count */
	i64 N = 4;
	/* Large interaction strength, corresponds to g(N-1) = 4 */
	f64 g0 = 4.0/3.0;

	struct nlse_settings settings = {
		.max_iterations = 1000,
		.max_quadgk_iters = 500,
		.abs_error_tol = 1e-14,

        .num_basis_funcs = 16,
		.basis = ho_basis,

		.zero_threshold = 1e-10,

		/* Specify mixing parameters */
		//.orbital_mixing = 0.95,
		.hamiltonian_mixing = 0.3,

		/* Call debug_callback every iteration */
		.measure_every = 1,
		.debug_callback = debug_callback,

		.gk=gk20
    };

	/* Solve GP equation */
	struct nlse_result res = grosspitaevskii(settings, component_count, &N, guesses, &g0);

	/* Compute energy */
	f64 Efull = grosspitaevskii_energy(settings, res.coeff_count, component_count, res.coeff, &N, &g0);
	printf("\nfull energy: %lf\n", Efull);

	sbmf_shutdown();
}
