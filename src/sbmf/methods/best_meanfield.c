#include <sbmf/methods/best_meanfield.h>

#include <sbmf/sbmf.h>
#include <sbmf/debug/log.h>
#include <sbmf/math/functions.h>
#include <sbmf/methods/gp2c.h>

#include <sbmf/methods/quadgk_vec.h>
#include <sbmf/math/harmonic_oscillator.h>

static f64 _lambda0;
static f64 _lambda;
static c64* _orbital_1_coeffs;
static c64* _orbital_2_coeffs;
static u32 _particle_count;
static u32 _n1;
static u32 _n2;

static void bestmf_perturbation(const u32 len, f64 out[static len],
                                f64 in_x[static len], const u32 component_count,
                                c64 in_u[static len*component_count]) {
	assert(component_count == 0);
	for (u32 i = 0; i < len; ++i) {
		out[i] = gaussian(in_x[i],0,0.2);
	}
}

static void op_1(const u32 len, f64 out[static len],
                                f64 in_x[static len], const u32 component_count,
                                c64 in_u[static len*component_count]) {
	for (u32 i = 0; i < len; ++i) {
		f64 abs1 = cabs(in_u[0*len + i]);
		f64 abs2 = cabs(in_u[1*len + i]);
		out[i] = 0.75 * _lambda * abs1*abs1 + 0.25 * _lambda * abs2*abs2;
	}
}

static void op_2(const u32 len, f64 out[static len],
                                f64 in_x[static len], const u32 component_count,
                                c64 in_u[static len*component_count]) {
	for (u32 i = 0; i < len; ++i) {
		f64 abs1 = cabs(in_u[0*len + i]);
		f64 abs2 = cabs(in_u[1*len + i]);
		out[i] = 0.75 * _lambda * abs2*abs2 + 0.25 * _lambda * abs1*abs1;
	}
}

static void guess_1(c64* out, u32 len) {
	for (u32 i = 0; i < len; ++i) {
		f64 c1 = sqrt((f64)_n1/(f64)_particle_count);
		f64 c2 = sqrt((f64)_n2/(f64)_particle_count);
		out[i]  = c1 * _orbital_1_coeffs[i] + c2 * _orbital_2_coeffs[i];
	}
}

static void guess_2(c64* out, u32 len) {
	for (u32 i = 0; i < len; ++i) {
		f64 c1 = sqrt((f64)_n1/(f64)_particle_count);
		f64 c2 = sqrt((f64)_n2/(f64)_particle_count);
		out[i]  = c2 * _orbital_2_coeffs[i] - c1 * _orbital_1_coeffs[i];
	}
}

static void ensure_structure_of_func(c64* a, c64* b, u32 len) {
	c64 sum = 0.0;
	for (u32 i = 0; i < len; ++i) {
		sum += conj(a[i])*b[i];
	}

	f64 n2_minus_n1 = (f64)_particle_count * creal(sum);
	/*log_info("n2 - n1 = %lf", n2_minus_n1);*/
	/* n2 - n1 = (N-n1)-n1 = N-2n1 => N-2n1 - (n2-n1) = 0
	 * => n1 = (N-(n2-n1))/2
	 */
	//c64 n1 = ((f64)_particle_count - n2_minus_n1)/2.0;
	//c64 n2 = (f64)_particle_count - n1;
	u32 n1 = (u32)(((f64)_particle_count - n2_minus_n1)/2.0);
	u32 n2 = _particle_count - n1;
	log_info("n1: %u -- n2: %u", n1,n2);
	//log_info("(%.2lf+%.2lfi) + (%.2lf+%.2lfi)", CCOMP(n1), CCOMP(n2));

	for (u32 i = 0; i < len; ++i) {
		f64 c1 = sqrt((f64)(n1)/ (f64)_particle_count);
		f64 c2 = sqrt((f64)(n2)/ (f64)_particle_count);

		a[i] = c1 * _orbital_1_coeffs[i] + c2 * _orbital_2_coeffs[i];
		b[i] = c2 * _orbital_2_coeffs[i] - c1 * _orbital_1_coeffs[i];
	}

	_n1 = n1;
	_n2 = n2;
}

struct best_meanfield_results find_best_meanfield_occupations(const u32 particle_count,
									 const f64 g,
                                     const u32 coeff_count,
                                     c64 orbital_1_coeffs[static coeff_count],
                                     c64 orbital_2_coeffs[static coeff_count],
                                     const f64 energy_1,
									 const f64 energy_2) {
	_n1 = particle_count;
	_n2 = particle_count - _n1;
	_particle_count = particle_count;
	_lambda0 = g;
	_lambda = g*(_particle_count - 1);
	_orbital_1_coeffs = orbital_1_coeffs;
	_orbital_2_coeffs = orbital_2_coeffs;

	struct gp2c_settings settings = {
		.num_basis_functions = coeff_count,
		.max_iterations = 1e6,
		.error_tol = 1e-8,
		.post_normalize_callback = ensure_structure_of_func,
		.ho_potential_perturbation = bestmf_perturbation,
		.gk = gk15,
	};
	struct gp2c_component component[2] = {
		[0] = {
			.guess = guess_1,
			.op = op_1,
		},
		[1] = {
			.guess = guess_2,
			.op = op_2,
		},
	};
	struct gp2c_result res = gp2c(settings, 2, component);

	const f64 energy = best_meanfield_energy(coeff_count, orbital_1_coeffs, orbital_2_coeffs, _n1, _n2, g);

	best_meanfield_energy(coeff_count, orbital_1_coeffs, orbital_2_coeffs, 0.6*particle_count, 0.4*particle_count, g);

	return (struct best_meanfield_results) {
		.occupation_1 = _n1,
		.occupation_2 = _n2,
		.energy = energy
	};
}











struct bme_integrand_params {
	u32 coeff_count;
	c64* orbital_1_coeffs;
	c64* orbital_2_coeffs;
};

void bme_integrand(f64* out, f64* in, u32 len, void* data) {
	struct bme_integrand_params* params = data;

	c64 sample_1_out[len];
	ho_sample(params->coeff_count, params->orbital_1_coeffs, len, sample_1_out, in);

	c64 sample_2_out[len];
	ho_sample(params->coeff_count, params->orbital_2_coeffs, len, sample_2_out, in);

	for (u32 i = 0; i < len; ++i) {
		f64 abs1 = cabs(sample_1_out[i]);
		f64 abs2 = cabs(sample_2_out[i]);
		out[i] = abs1*abs1*abs2*abs2;
	}
}

f64 best_meanfield_energy( const u32 coeff_count,
                           c64 orbital_1_coeffs[static coeff_count],
                           c64 orbital_2_coeffs[static coeff_count],
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
							  c64 state_coeff[static state_count*coeff_count]) {

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
    c64* coeff;
};

void energy_per_particle_integrand(f64* out, f64* in, u32 len, void* data) {
    struct energy_per_particle_integrand_params* params = data;

    c64 sample_out[len];
	ho_sample(params->coeff_count, params->coeff, len, sample_out, in);

    for (u32 i = 0; i < len; ++i) {
        f64 c = cabs(sample_out[i]);
        out[i] = c*c*c*c;
    }
}

f64 gp_energy_per_particle(const u32 particle_count,
						   const f64 interaction_strength,
						   const u32 coeff_count,
						   c64 coeff[static coeff_count]) {
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
        f64 c = cabs(coeff[i]);
        E += c*c*ho_eigenval(i);
    }

    E += 0.5*interaction_strength*(particle_count-1)*res.integral;

    return E;
}
