#include <sbmf/sbmf.h>
#include <stdio.h>

//#define PERTURBATION(x) 2*gaussian(x, 0, 0.2)
#define PERTURBATION(x) 0.0

void perturbation(const u32 len, f64 out[static len],
                                f64 in_x[static len], const u32 component_count,
                                f64 in_u[static len*component_count],
                                void* userdata) {
    for (u32 i = 0; i < len; ++i) {
        out[i] = PERTURBATION(in_x[i]);
    }
}





void log_callback(enum sbmf_log_level log_level, const char* msg) {
	printf("%s\n", msg);
}


int main() {
	sbmf_set_log_callback(log_callback);

	struct nlse_guess* guesses = NULL;

	i64 occupations[] = {2};

	struct nlse_settings settings = {
		.spatial_pot_perturbation = perturbation,
		.max_iterations = 1000,
		.max_integration_evals = 1e5,
		.error_tol = 1e-9,

        .num_basis_funcs = 32,
		.basis = ho_basis,

		.hamiltonian_mixing = 0.6,

		.zero_threshold = 1e-10,
		.gk=gk15
    };

	const u32 component_count = 1;

	f64 g0s[] = {3.6, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0};
	//f64 g0s[] = {/*-10.0, -9.5, -9.0, -8.5, -8.0, -7.5, -7.0, -6.5, -6.0, -5.5, -5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5,*/ 4.0, 4.1};
	//f64 g0s[] = {-10.0, -9.5, -9.0, -8.5, -8.0, -7.5, -7.0, -6.5, -6.0, -5.5, -5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5};
	for (u32 i = 0; i < sizeof(g0s)/sizeof(f64); ++i) {
		f64* g0 = &g0s[i];
		sbmf_init();
		struct nlse_result res = grosspitaevskii(settings, component_count, occupations, guesses, g0);

		f64 Egp = 0;
		if (res.converged)
			Egp = full_energy(settings, res.coeff_count, component_count, res.coeff, occupations, g0);

		{
			FILE* fd = fopen("out", "a");
			fprintf(fd, "%lf\t%lf\n",
					g0s[i],
					Egp
					);
			fclose(fd);
		}

		sbmf_shutdown();
	}

}
