#include <sbmf/sbmf.h>
#include <sbmf/methods/best_meanfield.h>
#include <sbmf/math/harmonic_oscillator.h>
#include <sbmf/methods/quadgk_vec.h>

#include <stdio.h>

//#define PERTURBATION(x) gaussian(x, 0, 0.2)
#define PERTURBATION(x) (-1.5015*sqrt(x*x - 1.5*1.5 + 1.5015*1.5015));

void perturbation(const u32 len, f64 out[static len],
                                f64 in_x[static len], const u32 component_count,
                                f64 in_u[static len*component_count],
								void* userdata) {
    assert(component_count == 0);
    for (u32 i = 0; i < len; ++i) {
        out[i] = PERTURBATION(in_x[i]);
    }
}

struct inner_product_integrand_params {
	f64* coeff_a;
	f64* coeff_b;
	u32  coeff_count;
	struct basis basis;
};

void inner_product_integrand(f64* out, f64* in, u32 len, void* data) {
	struct inner_product_integrand_params* p = data;

	f64 sample_a[len];
	p->basis.sample(p->coeff_count, p->coeff_a, len, sample_a, in);

	f64 sample_b[len];
	p->basis.sample(p->coeff_count, p->coeff_b, len, sample_b, in);

	for (u32 i = 0; i < len; ++i) {
		out[i] = sample_a[i]*sample_b[i];
	}
}





void log_callback(enum sbmf_log_level log_level, const char* msg) {
	printf("%s\n", msg);
}

int main(int argc, char** argv) {
	u32 particle_count = 100;
	f64 g0 = (0.5)/(particle_count-1);

	sbmf_set_log_callback(log_callback);
	sbmf_init();

	struct nlse_settings settings = {
        .num_basis_funcs = 16,
        .max_iterations = 1e8,
        .error_tol = 1e-9,
        .spatial_pot_perturbation = perturbation,
        .gk = gk15,
		.basis = ho_basis,
		.zero_threshold = 1e-10,
    };

	struct nlse_result res = best_meanfield(settings, particle_count, g0);

	{
		struct inner_product_integrand_params p = {
			.coeff_a = &res.coeff[0],
			.coeff_b = &res.coeff[res.coeff_count],
			.coeff_count = res.coeff_count,
			.basis = settings.basis,
		};

		integration_settings int_settings = {
			.gk = gk15,
			.abs_error_tol = 1e-10,
			.max_evals = 1e5,
			.userdata = &p,
		};

		integration_result ires = quadgk_vec(inner_product_integrand, -INFINITY, INFINITY, int_settings);

		/*
		 * <p1|p2> = (n2-n1)/N = (N-n1-n1)/N
		 * 	=> n1 = (N-N<p1|p2>)/2
		 */
		f64 n1 = particle_count * (1.0 - ires.integral) / 2.0;
		printf("n1: %lf\n", n1);
	}

	sbmf_shutdown();
}
