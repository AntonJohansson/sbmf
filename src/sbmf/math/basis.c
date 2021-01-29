//#define OMEGA (1.0)
f64 OMEGA = 1.0;

/* Currently doesnt handle 2d/3d/... case */
void ho_eigenfunc(const u32 n, const u32 len, f64 out[static len], f64 in[static len]) {
	/*
	 * HN = 2xH_{N-1} - 2(N-1)H_{N-2}
	 * psiN = 1/sqrt(2^N * N!) * (1/pi^(1/4)) * exp(-x^2/2) * HN
	 */
	assert(n < 270);

	const f64 pi_factor = pow(OMEGA/M_PI,0.25);
	const f64 normalization_factor = pi_factor / sqrt(pow(2,n) * factorial_128(n));

	const f64 sqrt_omega = sqrt(OMEGA);
	const f64 half_omega = OMEGA/2.0;

	for (u32 i = 0; i < len; ++i) {
		f64 x = in[i];

		f64 init_value = exp(-half_omega*x*x) * normalization_factor;
		f64 H[3] = {
			init_value,
			init_value,
			init_value,
		};

		for (u32 j = 1; j <= n; ++j) {
			H[2] = 2*(sqrt_omega*x*H[1] - (j-1)*H[0]);
			H[0] = H[1];
			H[1] = H[2];
		}

		out[i] = H[2];
	}
}

/* energy eigenvalue */
f64 ho_eigenval(const u32 n) {
	/*
	f64 sum = 0.0;
	for (i32 i = 0; i < dims; ++i)
		sum += states[i];

	return sum + dims/2.0;
	*/
	return OMEGA*((f64)n + 0.5);
}

/* Currently doesnt handle 2d/3d/... case */
void ho_sample(const u32 coeff_count,
		f64 coeffs[static coeff_count],
		const u32 len,
		f64 out[static len],
		f64 in[static len]) {
	for (u32 i = 0; i < len; ++i)
		out[i] = 0;

	f64 eigenfunc_out[len];
	for (u32 i = 0; i < coeff_count; ++i) {
		ho_eigenfunc(i, len, eigenfunc_out, in);
		for (u32 j = 0; j < len; ++j) {
			out[j] += coeffs[i]*eigenfunc_out[j];
		}
	}
}

struct basis ho_basis = {
	.eigenfunc = ho_eigenfunc,
	.eigenval  = ho_eigenval,
	.sample    = ho_sample
};









#if 0
#include <gsl/gsl_sf_hermite.h>

static void ho_gsl_eigenfunc(const u32 n, const u32 len, f64 out[static len], f64 in[static len]) {
	for (u32 i = 0; i < len; ++i) {
		out[i] = gsl_sf_hermite_func(n, in[i]);
	}
}

static f64 ho_gsl_eigenval(const u32 n) {
	return n + 0.5;
}

static void ho_gsl_sample(const u32 coeff_count,
		f64 coeffs[static coeff_count],
		const u32 len,
		f64 out[static len],
		f64 in[static len]) {
	for (u32 i = 0; i < len; ++i)
		out[i] = 0;

	f64 eigenfunc_out[len];
	for (u32 i = 0; i < coeff_count; ++i) {
		ho_gsl_eigenfunc(i, len, eigenfunc_out, in);
		for (u32 j = 0; j < len; ++j) {
			out[j] += coeffs[i]*eigenfunc_out[j];
		}
	}
}


static struct basis ho_gsl_basis = {
	.eigenfunc = ho_gsl_eigenfunc,
	.eigenval  = ho_gsl_eigenval,
	.sample    = ho_gsl_sample
};
#endif




























f64 ho_potential(f64* v, i32 n, c64 u) {
	SBMF_UNUSED(u);

	f64 temp = 0.0;
	for (i32 i = 0; i < n; ++i) {
		temp += v[i]*v[i];
	}
	return (OMEGA*OMEGA)*temp*0.5;
}

/* Currently doesnt handle 2d/3d/... case */
void ho_potential_vec(f64* out, f64* in, u32 len) {
	for (u32 i = 0; i < len; ++i) {
		out[i] = 0.5*in[i]*in[i];
	}
}
