#include "best_meanfield.h"

#include <sbmf/sbmf.h>
#include <sbmf/debug/log.h>
#include <sbmf/math/functions.h>
#include "find_groundstate.h"

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

	log_info("sum: %lf + %lfi", CCOMP(sum));
	f64 n2_minus_n1 = (f64)_particle_count * creal(sum);
	log_info("n2 - n1 = %lf", n2_minus_n1);
	/* n2 - n1 = (N-n1)-n1 = N-2n1 => N-2n1 - (n2-n1) = 0
	 * => n1 = (N-(n2-n1))/2
	 */
	//u32 n1 = lroundl(((f64)_particle_count - n2_minus_n1)/2.0);
	u32 n1 = (u32)(((f64)_particle_count - n2_minus_n1)/2.0);
	u32 n2 = _particle_count - n1;
	log_info("%u,%u", n1,n2);

	for (u32 i = 0; i < len; ++i) {
		f64 c1 = sqrt((f64)n1/(f64)_particle_count);
		f64 c2 = sqrt((f64)n2/(f64)_particle_count);

		a[i] = c1 * _orbital_1_coeffs[i] + c2 * _orbital_2_coeffs[i];
		b[i] = c2 * _orbital_2_coeffs[i] - c1 * _orbital_1_coeffs[i];
	}

	_n1 = n1;
	_n2 = n2;
}

void find_best_meanfield_occupations(const u32 particle_count,
									 const f64 g,
                                     const u32 coeff_count,
                                     c64 orbital_1_coeffs[static coeff_count],
                                     c64 orbital_2_coeffs[static coeff_count],
                                     const f64 energy_1,
									 const f64 energy_2) {
	_n1 = particle_count;
	_n2 = 0;
	_lambda0 = g;
	_lambda = g*(_particle_count - 1);
	_particle_count = particle_count;
	_orbital_1_coeffs = orbital_1_coeffs;
	_orbital_2_coeffs = orbital_2_coeffs;

	struct gp2c_settings settings = {
		.num_basis_functions = coeff_count,
		.max_iterations = 1e7,
		.error_tol = 1e-8,
		.post_normalize_callback = ensure_structure_of_func,
		.ho_potential_perturbation = bestmf_perturbation,
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

	best_meanfield_energy(particle_count, coeff_count, orbital_1_coeffs, orbital_2_coeffs, energy_1, energy_2, _n1, _n2);
	best_meanfield_energy(particle_count, coeff_count, orbital_1_coeffs, orbital_2_coeffs, energy_1, energy_2, 6, 4);
	best_meanfield_energy(particle_count, coeff_count, orbital_1_coeffs, orbital_2_coeffs, energy_1, energy_2, 7, 3);
}

struct bme_integrand_params {
	u32 coeff_count;
	c64* orbital_1_coeffs;
	c64* orbital_2_coeffs;
};

void bme_integrand(f64* out, f64* in, u32 len, void* data) {
	struct bme_integrand_params* params = data;

	c64 sample_1_out[params->coeff_count];
	hob_sample_vec(params->orbital_1_coeffs, params->coeff_count, sample_1_out, in, len);

	c64 sample_2_out[params->coeff_count];
	hob_sample_vec(params->orbital_2_coeffs, params->coeff_count, sample_2_out, in, len);

	for (u32 i = 0; i < len; ++i) {
		f64 abs1 = cabs(sample_1_out[i]);
		f64 abs2 = cabs(sample_2_out[i]);
		out[i] = abs1*abs1*abs2*abs2;
	}
}

void best_meanfield_energy(const u32 particle_count,
                           const u32 coeff_count,
                           c64 orbital_1_coeffs[static coeff_count],
                           c64 orbital_2_coeffs[static coeff_count],
						   const f64 energy_1,
						   const f64 energy_2,
                           const u32 occupation_1,
                           const u32 occupation_2) {

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

	f64 E = occupation_1*energy_1 + occupation_2*energy_2 + 2*_lambda0*occupation_1*occupation_2*res.integral;
	log_info("(%u,%u) Energy: %lf", occupation_1, occupation_2, E);
	log_info("(%u,0) Energy: %lf", occupation_1+occupation_2, (occupation_1+occupation_2)*energy_1);
	log_info("(0,%u) Energy: %lf", occupation_1+occupation_2, (occupation_1+occupation_2)*energy_2);
}
