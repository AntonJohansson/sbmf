#include <sbmf/methods/best_meanfield.h>

#include <sbmf/sbmf.h>
#include <sbmf/math/functions.h>

#include <sbmf/methods/quadgk_vec.h>
#include <sbmf/math/harmonic_oscillator.h>

struct bme_integrand_params {
	u32 coeff_count;
	f64* orbital_1_coeffs;
	f64* orbital_2_coeffs;
};

void bme_integrand(f64* out, f64* in, u32 len, void* data) {
	struct bme_integrand_params* params = data;

	f64 sample_1_out[len];
	ho_sample(params->coeff_count, params->orbital_1_coeffs, len, sample_1_out, in);

	f64 sample_2_out[len];
	ho_sample(params->coeff_count, params->orbital_2_coeffs, len, sample_2_out, in);

	for (u32 i = 0; i < len; ++i) {
		f64 abs1 = fabs(sample_1_out[i]);
		f64 abs2 = fabs(sample_2_out[i]);
		out[i] = abs1*abs1*abs2*abs2;
	}
}

f64 best_meanfield_energy( const u32 coeff_count,
                           f64 orbital_1_coeffs[static coeff_count],
                           f64 orbital_2_coeffs[static coeff_count],
                           const u32 occupation_1,
                           const u32 occupation_2,
						   const f64 interaction_strength) {
	struct bme_integrand_params params = {
		.coeff_count = coeff_count,
		.orbital_1_coeffs = orbital_1_coeffs,
		.orbital_2_coeffs = orbital_2_coeffs,
	};
	struct integration_settings int_settings = {
		.gk = gk7,
		.abs_error_tol = 1e-10,
		.rel_error_tol = 1e-10,
		.max_evals = 1e5,
		.userdata = &params,
	};

	struct integration_result res = quadgk_vec(bme_integrand, -INFINITY, INFINITY, int_settings);
	assert(res.converged);

	f64 E =
		occupation_1*gp_energy_per_particle(occupation_1, interaction_strength, coeff_count, orbital_1_coeffs) +
		occupation_2*gp_energy_per_particle(occupation_2, interaction_strength, coeff_count, orbital_2_coeffs) +
		0*
		2*interaction_strength*occupation_1*occupation_2*res.integral;
	return E;
}































f64 best_meanfield_energy_new(const u32 state_count,
							  const u32 occupation[static state_count],
							  const f64 coupling[static state_count*state_count],
							  const u32 coeff_count,
							  f64 state_coeff[static state_count*coeff_count]) {

	struct bme_integrand_params params = {
		.coeff_count = coeff_count,
		.orbital_1_coeffs = NULL,
		.orbital_2_coeffs = NULL,
	};
	struct integration_settings int_settings = {
		.gk = gk7,
		.abs_error_tol = 1e-10,
		.rel_error_tol = 1e-10,
		.max_evals = 1e5,
		.userdata = &params,
	};

	f64 E = 0.0;
	for (u32 i = 0; i < state_count; ++i) {
		E += occupation[i]*gp_energy_per_particle(occupation[i], coupling[i*state_count], coeff_count, &state_coeff[i*coeff_count]);
		for (u32 j = 0; j < state_count; ++j) {
			/* This case is included in energy of state i */
			if (i == j)
				continue;

			/* Calculate exchange integral */
			params.orbital_1_coeffs = &state_coeff[i*coeff_count];
			params.orbital_2_coeffs = &state_coeff[j*coeff_count];
			struct integration_result res = quadgk_vec(bme_integrand, -INFINITY, INFINITY, int_settings);
			assert(res.converged);

			f64 g = coupling[i*state_count + j];
			E += g*occupation[i]*occupation[j]*res.integral;
		}
	}

	return E;
}













struct energy_per_particle_integrand_params {
    u32 coeff_count;
    f64* coeff;
};

void energy_per_particle_integrand(f64* out, f64* in, u32 len, void* data) {
    struct energy_per_particle_integrand_params* params = data;

    f64 sample_out[len];
	ho_sample(params->coeff_count, params->coeff, len, sample_out, in);

    for (u32 i = 0; i < len; ++i) {
        f64 c = fabs(sample_out[i]);
        out[i] = c*c*c*c;
    }
}

f64 gp_energy_per_particle(const u32 particle_count,
						   const f64 interaction_strength,
						   const u32 coeff_count,
						   f64 coeff[static coeff_count]) {
    static struct integration_settings settings = {
        .abs_error_tol = 1e-10,
        .rel_error_tol = 1e-10,
        .max_evals = 1e5,
    };
    struct energy_per_particle_integrand_params params = {
        .coeff_count = coeff_count,
        .coeff = coeff
    };
    settings.userdata = &params;

    struct integration_result res = quadgk_vec(energy_per_particle_integrand, -INFINITY, INFINITY, settings);

    f64 E = 0.0;
    for (u32 i = 0; i < coeff_count; ++i) {
        f64 c = fabs(coeff[i]);
        E += c*c*ho_eigenval(i);
    }

    E += 0.5*interaction_strength*(particle_count-1)*res.integral;

    return E;
}
