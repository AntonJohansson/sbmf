#include <sbmf/sbmf.h>
#include <stdio.h>

/* Logging function for sbmf to use */
void log_callback(enum sbmf_log_level log_level, const char* msg) {
	printf("%s\n", msg);
}

/* We later request spatial guesses for SCF and supply these functions as that
 * spatial guess, meaning our guesses will be gaussians centered at -1, +1.
 */
void gaussian0(f64* out, f64* in, u32 len, void* data) {
    for (u32 i = 0; i < len; ++i)
        out[i] = gaussian(in[i] + 1.0, 0.0, 0.1);
}
void gaussian1(f64* out, f64* in, u32 len, void* data) {
    for (u32 i = 0; i < len; ++i)
        out[i] = gaussian(in[i] - 1.0, 0.0, 0.1);
}

int main() {
	sbmf_set_log_callback(log_callback);

	/* Use gaussian guesses */
    struct nlse_guess gaussian_guesses[] = {
        [0] = {
            .type = SPATIAL_GUESS,
            .data.spatial_guess = gaussian0,
        },
        [1] = {
            .type = SPATIAL_GUESS,
            .data.spatial_guess = gaussian1,
        },
    };
	/* Use default guesses */
	struct nlse_guess* default_guesses = NULL;

	/* Particle count to use */
	i64 N = 100;
	/* Values of lambda = g(N-1) to use. As we keep N fix, this scales interaction strength */
	f64 l0s[] = {-2,-1.75,-1.5,-1.25,-1,-0.75,-0.5,-0.25,0.25,0.5,0.75,1,1.25,1.5,1.75,2};

	struct nlse_settings settings = {
		.max_iterations = 1e5,
		.max_quadgk_iters = 500,
		.abs_error_tol = 1e-14,

		.num_basis_funcs = 64,
		.basis = ho_basis,
		.hamiltonian_mixing = 0.7,

		.zero_threshold = 1e-10,
		.gk=gk20
    };

	const u32 component_count = 1;

	for (u32 j = 0; j < sizeof(l0s)/sizeof(l0s[0]); ++j) {
		sbmf_init();

		struct nlse_result res;
		f64 bmf_gaussian_res;
		f64 bmf_default_res;

		/* Compute interaction strength to use for given lambda and N */
		f64 g0 = l0s[j]/((f64)N-1.0);

		/* Solve GP equation and compute energy */
		res = grosspitaevskii(settings, component_count, &N, default_guesses, &g0);
		f64 Egp = grosspitaevskii_energy(settings, res.coeff_count, component_count, res.coeff, &N, &g0);

		/* Compute fractional occupations so we know the GP ground state is an
		 * accurate starting point. */
		bmf_gaussian_res = bestmf_find_fractional_occupation(settings, N, g0, gaussian_guesses);
		bmf_default_res = bestmf_find_fractional_occupation(settings, N, g0, default_guesses);

		/* Apply perturbation theory */
		struct pt_result rs_ptres = rspt_1comp_cuda(&settings, res, 0, g0, N);
		struct pt_result en_ptres = enpt_1comp_cuda(&settings, res, 0, g0, N);

		/* Save results to file */
		{
			FILE* fd = fopen("out", "a");
			fprintf(fd, "%.10lf\t%.10lf\t%.10f\t%.10f\t%.10f\t%.10f\t%.10lf\t%.10lf\t%.10lf\t%.10lf\n",
					l0s[j],
					Egp,
					bmf_gaussian_res,
					bmf_default_res,
					rs_ptres.E0+rs_ptres.E1														- Egp,
					rs_ptres.E0+rs_ptres.E1+rs_ptres.E2								- Egp,
					rs_ptres.E0+rs_ptres.E1+rs_ptres.E2+rs_ptres.E3		- Egp,
					en_ptres.E0+en_ptres.E1					   								- Egp,
					en_ptres.E0+en_ptres.E1+en_ptres.E2								- Egp,
					en_ptres.E0+en_ptres.E1+en_ptres.E2+en_ptres.E3		- Egp
				   );
			fclose(fd);
		}


		sbmf_shutdown();
	}

}
