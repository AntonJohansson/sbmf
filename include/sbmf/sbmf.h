#pragma once

#include <stdint.h>
#include <complex.h>
#include <math.h>
#include <stdbool.h>

/* Redefinition of numerical types for personal reasons */
typedef int8_t 				i8;
typedef uint8_t 			u8;
typedef int16_t 			i16;
typedef uint16_t 			u16;
typedef int32_t 			i32;
typedef uint32_t 			u32;
typedef int64_t 			i64;
typedef uint64_t 			u64;
typedef __int128_t  		i128;
typedef __uint128_t 		u128;
typedef float 				f32;
typedef double 				f64;
typedef long double 		f128;
typedef float complex 		c32;
typedef double complex 		c64;
typedef long double complex c128;

/*
 * Initialization
 */

void sbmf_init();
void sbmf_shutdown();

/* Logging */
enum sbmf_log_level {
	SBMF_LOG_LEVEL_INFO    = 0,
	SBMF_LOG_LEVEL_WARNING = 1,
	SBMF_LOG_LEVEL_ERROR   = 2,
	SBMF_LOG_LEVEL_PANIC   = 3,
};

typedef void sbmf_log_callback_func(enum sbmf_log_level, const char*);

void sbmf_set_log_callback(sbmf_log_callback_func* func);

/*
 * Math functions
 */

static inline f64 gaussian(f64 x, f64 mu, f64 sigma) {
	return 1.0/(sigma*sqrt(2.0*M_PI)) * exp(-(x-mu)*(x-mu)/(2*sigma*sigma));
}

static inline void f64_normalize(f64* out, f64* data, u32 size) {
	f64 sum = 0.0;

#pragma omp parallel for shared(data) reduction(+: sum)
	for (u32 i = 0; i < size; ++i) {
		f64 absval = fabs(data[i]);
		sum += absval*absval;
	}

	f64 scaling = 1.0/sqrt(sum);
#pragma omp parallel for
	for (u32 i = 0; i < size; ++i) {
		out[i] = data[i] * scaling;
	}
}

/*
 * Matrix
 */

struct hermitian_bandmat {
	f64* data;
	u32 bandcount; /* rows */
	u32 size; /* cols */
};

struct hermitian_bandmat hermitian_bandmat_new(u32 bandcount, u32 size);
struct hermitian_bandmat hermitian_bandmat_new_zero(u32 bandcount, u32 size);

#define U32MIN(a,b) \
	((a < b) ? a : b)

#define HERMITIAN_BANDMAT_FOREACH(bm, r,c) 						\
	for (u32 r = 0; r < bm.size; ++r)									\
		for (u32 c = r; c < U32MIN(bm.size, r+bm.bandcount); ++c)

static inline u32 hermitian_bandmat_index(struct hermitian_bandmat bm, u32 row, u32 col) {
	return bm.size * (bm.bandcount - 1 + (row - col)) + col;
}

void hermitian_bandmat_mulv(f64* ans_vec, struct hermitian_bandmat bm, f64* vec);

enum which_eigenpairs {
	EV_LARGEST_MAG 		= 0,
	EV_SMALLEST_MAG 	= 1,
	EV_LARGEST_RE		= 2,
	EV_SMALLEST_RE 		= 3,
	EV_LARGEST_IM 		= 4,
	EV_SMALLEST_IM		= 5,
	EV_LARGEST			= 6,
	EV_SMALLEST 		= 7,
	EV_BOTH				= 8,
};

struct eigen_result_real {
	f64* eigenvalues;
	f64* eigenvectors;
	u32 num_eigenpairs;
	u32 points_per_eigenvector;
};

/* Find _all_ eigenpairs for a dense, symmetric,
 * upper tridiagonal matrix.
 */
struct eigen_result_real find_eigenpairs_full_real(struct hermitian_bandmat bm);

/* Find _some_ eigenpairs (specified by the enum which_eigenpairs)
 * for a dense, upper tridiagonal matrix.
 */
struct eigen_result_real find_eigenpairs_sparse_real(struct hermitian_bandmat bm, u32 num_eigenvalues, enum which_eigenpairs which);

/*
 * QUADGK
 */

#define MAX_GAUSS_POINTS 20
struct gk_data {
	f64 kronod_nodes[MAX_GAUSS_POINTS+1];
	f64 kronod_weights[MAX_GAUSS_POINTS+1];
	f64 gauss_weights[(MAX_GAUSS_POINTS+1)/2];
	u32 kronod_size;
	u32 gauss_size;
};

extern struct gk_data gk7;
extern struct gk_data gk10;
extern struct gk_data gk15;
extern struct gk_data gk20;

typedef f64 integrand(f64,void*);
typedef struct integration_settings {
	struct gk_data gk;

	// Aboslute and relative error tolarences.
	f64 abs_error_tol;
	f64 rel_error_tol;

	// The maximum allowed function evaluations of
	// the supplied integrand.
	i32 max_evals;

	void* userdata;
} integration_settings;

typedef struct integration_result {
	f64 integral;
	f64 error;
	i32 performed_evals;
	bool converged;
} integration_result;

typedef void integrand_vec(f64*,f64*,u32,void*);

integration_result quadgk_vec(integrand_vec* f, f64 start, f64 end, integration_settings settings);

/*
 * Basis
 */

/* A brief example to show how this structure will be used.
 * Consider the polynomial basis (in position representation)
 *
 * 		<x|0> = 1,
 * 		<x|1> = x,
 * 		<x|2> = x^2,
 * 		...,
 *
 * a set of coefficients such as (1, 2, 3, 4, 5) expressed
 * in this basis would correspond to
 *
 *	1<x|0> + 2<x|1> + 3<x|2> + 4<x|3> + 5<x|4>
 *	 = 1 + 2x + 3x^2 + 4x^3 + 5x^4,
 *
 * the conversion of a set of coefficients to a position
 * representation in some basis is exacly what the "sample"
 * function below does. It evaluates a set of coeffs. at
 * a particular x point.
 */

/* assuming 1D */
typedef void basis_eigenfunc_func(const u32 n, const u32 len,
		f64 out[static len], f64 in[static len]);
typedef f64  basis_energy_eigenval_func(const u32 n);
typedef void basis_sample_func(
		const u32 coeff_count,
		f64 coeffs[static coeff_count],
		const u32 len,
		f64 out[static len],
		f64 in[static len]);

struct basis {
	basis_eigenfunc_func*       eigenfunc;
	basis_energy_eigenval_func* eigenval;
	basis_sample_func*          sample;
};

/*
 * Harmonic oscillator stuff
 */

void ho_eigenfunc(const u32 n, const u32 len, f64 out[static len], f64 in[static len]);
f64 ho_eigenval(const u32 n);
void ho_sample(const u32 coeff_count,
		f64 coeffs[static coeff_count],
		const u32 len,
		f64 out[static len],
		f64 in[static len]);

f64 ho_potential(f64* v, i32 n, c64 u);
void ho_potential_vec(f64* out, f64* in, u32 len);

extern struct basis ho_basis;

/*
 * NLSE solving
 */

/*
 * Solves a general, non-linear system of
 * Schrödinger equations by the means of
 * iteration until self-consistency in
 * a specified basis.
 */

typedef void nlse_operator_func(const u32 len, f64 out[static len],
		f64 in_x[static len], const u32 component_count,
		f64 in_u[static len*component_count],
		void* userdata);
typedef void nlse_callback(c64* a, c64* b, u32 len);

struct nlse_settings;
struct nlse_result;
typedef void nlse_debug_callback(struct nlse_settings, struct nlse_result);

struct nlse_settings {
	u32 max_iterations;
	f64 error_tol;

	/* Separating out the spatial potentiential
	 * allows for optimizations */
	nlse_operator_func* spatial_pot_perturbation;

	u32 measure_every;
	nlse_debug_callback* debug_callback;

	void* post_normalize_userdata;
	nlse_debug_callback* post_normalize_callback;

	/* Choice of GK rule used for numerical integration internally */
	struct gk_data gk;

	/* Choice of basis to solve to problem in */
	u32 num_basis_funcs;
	struct basis basis;

	/* Everything below the zero_threshold is considered
	 * 0 in the hamiltonian. */
	f64 zero_threshold;
};

/* Initial Guess */
typedef void nlse_coeff_guess_func(f64* out, u32 len, u32 component);
typedef void nlse_spatial_guess_func(f64* out, f64* in, u32 len, void* data);

struct nlse_guess {
	enum {
		DEFAULT_GUESS 	= 0,
		SPATIAL_GUESS 	= 1,
		COEFF_GUESS   	= 2
	} type;
	union {
		nlse_coeff_guess_func* 	 coeff_guess;
		nlse_spatial_guess_func* spatial_guess;
	} data;
};

struct nlse_component {
	struct nlse_guess guess;
	nlse_operator_func* op;
	void* userdata;
};

struct nlse_result {
	u32 iterations;
	u32 component_count;
	u32 coeff_count;
	f64* coeff;
	f64* error;
	f64* energy;
	struct hermitian_bandmat* hamiltonian;
	bool converged;
};

struct nlse_result nlse_solver(struct nlse_settings settings, const u32 component_count, struct nlse_component components[static component_count]);

/*
 * Gross-pitaevskii solving
 */

struct gp_settings {
	nlse_operator_func* pot;

	u32 num_basis_funcs;
	struct basis basis;

	u32 component_count;
	u32* occupations; 				/* [static component_count] */
	struct nlse_guess* guesses; 	/* [static component_count] */
	f64* g0; 						/* [static component_count*component_count] */

	nlse_debug_callback* debug_callback;
	u32 measure_every;

	f64 zero_threshold;
};

struct nlse_result grosspitaevskii(struct nlse_settings settings,
		const u32 comp_count,
		u32 occupations[static comp_count],
		struct nlse_guess guesses[static comp_count],
		f64 g0[static comp_count*comp_count]);

f64 full_energy_naked(struct nlse_settings settings,
		const u32 coeff_count, const u32 comp_count,
		f64 coeff[static coeff_count*comp_count],
		u32 occupations[static comp_count],
		f64 g0[static comp_count*comp_count]
		);

/*
 * Best mean-field
 */

struct bestmf_result {
	f64 energy;
	u32 coeff_count;
	u32 comp_count;
	f64* coeff;
	f64 n1;
	f64 n2;
};

struct bestmf_result best_meanfield(struct nlse_settings settings,
		const u32 particle_count, f64 g0, struct nlse_guess* guesses);

struct bestmf_2comp_result {
	f64 energy;
	u32 coeff_count;
	u32 comp_count;
	f64* coeff;
	f64 n1;
	f64 n2;
	f64 n3;
	f64 n4;
};
struct bestmf_2comp_result best_meanfield_2comp(struct nlse_settings settings,
		const u32 particle_count,
		f64 g0[static 2*2],
		struct nlse_guess* guesses);

/*
 * Perturbation theory
 */

struct pt_result {
	f64 E0, E1, E2;
};

struct pt_result rayleigh_schroedinger_pt(struct nlse_result res, f64* g0, u32* particle_count);
