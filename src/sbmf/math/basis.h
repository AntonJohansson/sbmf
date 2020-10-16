#pragma once

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
typedef void basis_sample_func(u32 len, c64 coeffs[static len], f64 sample_x);
typedef void basis_eigenfunc(u32 state_index, u32 len, f64 out[static len], f64 in[static len]);
typedef f64  basis_energy_eigenval_func(u32 state_index);

struct basis {
	basis_sample_func* sample;
};
