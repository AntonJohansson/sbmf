#include "best_meanfield.h"

#include <sbmf/sbmf.h>
#include <sbmf/debug/log.h>
#include <sbmf/math/functions.h>
#include "find_groundstate.h"

static f64 lambda = -0.25 * (10 - 1);

static c64* _orbital_1_coeffs;
static c64* _orbital_2_coeffs;
static u32 _particle_count;
static u32 _n1;
static u32 _n2;

static void op1(f64* out, f64* in_x, c64* in_1, c64* in_2, u32 len) {
	for (u32 i = 0; i < len; ++i) {
		f64 abs1 = cabs(in_1[i]);
		f64 abs2 = cabs(in_2[i]);
		out[i] = gaussian(in_x[i],0,0.2) + 0.75 * lambda * abs1*abs1 + 0.25 * lambda * abs2*abs2;
	}
}

static void op2(f64* out, f64* in_x, c64* in_1, c64* in_2, u32 len) {
	for (u32 i = 0; i < len; ++i) {
		f64 abs1 = cabs(in_1[i]);
		f64 abs2 = cabs(in_2[i]);
		out[i] = gaussian(in_x[i],0,0.2) + 0.75 * lambda * abs2*abs2 + 0.25 * lambda * abs1*abs1;
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

	f64 n2_minus_n1 = (f64)_particle_count * cabs(sum);
	u32 n1 = lroundl(((f64)_particle_count - n2_minus_n1)/2.0);
	u32 n2 = _particle_count - n1;

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
									const u32 coeff_count,
									c64 orbital_1_coeffs[static coeff_count],
									c64 orbital_2_coeffs[static coeff_count]) {
	_n1 = particle_count/2;
	_n2 = particle_count - _n1;
	_particle_count = particle_count;
	_orbital_1_coeffs = orbital_1_coeffs;
	_orbital_2_coeffs = orbital_2_coeffs;

	struct gp2c_settings settings = {
		.num_basis_functions = coeff_count,
		.max_iterations = 1e7,
		.error_tol = 1e-10,
		.guess_a = guess_1,
		.guess_b = guess_2,
		.post_normalize_callback = ensure_structure_of_func,
	};
	struct gp2c_result res = gp2c(settings, op1, op2);

	log_info("n1: %u -- n2: %u", _n1, _n2);

	/* Need to copy gp2c and force scaling relation between bases */
}
