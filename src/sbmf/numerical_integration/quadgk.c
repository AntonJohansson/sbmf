#include "quadgk.h"
#include <sbmf/common/prioqueue.h>
#include <string.h> // memcpy
#include <stdio.h> // printf
#include <math.h> // INFINITY and isinf

const f64 xd7[8] = {
	-9.9145537112081263920685469752598e-01,
	-9.4910791234275852452618968404809e-01,
	-8.6486442335976907278971278864098e-01,
	-7.415311855993944398638647732811e-01,
	-5.8608723546769113029414483825842e-01,
	-4.0584515137739716690660641207707e-01,
	-2.0778495500789846760068940377309e-01,
	0.0
};

const f64 wd7[8] = {
	2.2935322010529224963732008059913e-02,
	6.3092092629978553290700663189093e-02,
	1.0479001032225018383987632254189e-01,
	1.4065325971552591874518959051021e-01,
	1.6900472663926790282658342659795e-01,
	1.9035057806478540991325640242055e-01,
	2.0443294007529889241416199923466e-01,
	2.0948214108472782801299917489173e-01
};

const f64 gwd7[4] = {
	1.2948496616886969327061143267787e-01,
	2.797053914892766679014677714229e-01,
	3.8183005050511894495036977548818e-01,
	4.1795918367346938775510204081658e-01
};

typedef struct {
	f64 integral;
	f64 error;

	f64 start;
	f64 end;
} segment;

static inline f64 norm(f64 r) {
	return r*r;
}

segment evaluate_rule(integrand* f, f64 start, f64 end, integration_settings settings) {
	static i32 xsize = 8;
	static const f64* x = xd7;
	static const f64* w = wd7;

	static i32 gsize = 4;
	static const f64* gw = gwd7;

	f64 midpoint = 0.5 * (end - start);


	// n1 is:
	// 	0 -- if even order
	// 	1 -- if odd order
	i32 sample_count = xsize;
	i32 n1 = 1 - (sample_count & 1);

	f64 fg, fk;
	f64 Ig = 0, Ik = 0;

	i32 size = gsize - n1;
	for (i32 i = 0; i < size; ++i) {
		fg = f(start + (1+x[2*i+1])*midpoint, settings.userdata) + f(start + (1-x[2*i+1])*midpoint, settings.userdata);
		fk = f(start + (1+x[2*i-1+1])*midpoint, settings.userdata) + f(start + (1-x[2*i-1+1])*midpoint, settings.userdata);
		Ig += fg * gw[i];
		Ik += fg * w[2*i+1] + fk * w[2*i-1+1];
	}

	// apparently last point has to be handled differently depending on
	// wether or not the order of integration is odd/even.
	if (n1 == 0) {
		Ik += f(start + midpoint, settings.userdata) * w[7]; // 7 == end
	} else {
		f64 f0 = f(start + midpoint, settings.userdata);
		Ig += f0 * gw[3]; // 3 == end
		Ik += f0 * w[7] + (f(start + (1+x[7-1])*midpoint, settings.userdata) + f(start + (1-x[7-1])*midpoint, settings.userdata)) * w[7-1]; // 7-1 == end-1
	}

	f64 Ik_s = Ik*midpoint;
	f64 Ig_s = Ig*midpoint;
	f64 error = norm(Ik_s - Ig_s);

	// add nan/inf check on error

	return (segment){Ik_s, error, start, end};
}

// Sort segments by error in descending order.
static bool compare_segments(void* a, void* b) {
	return (((segment*)a)->error > ((segment*)b)->error);
}

static inline bool should_exit(integration_result result, integration_settings settings) {
	return (result.error <= settings.abs_error_tol ||
					result.error <= settings.rel_error_tol*norm(result.integral) ||
					result.performed_evals >= settings.max_evals);
}











typedef struct {
	integrand* original_f;
	f64 original_start;
	f64 original_end;
	void* userdata;
} coordinate_transform_data;

static inline f64 integrand_both_endpoints_inf(f64 t, void* data) {
	coordinate_transform_data* tdata = data;
	// Both endpoints infinity
	// Do the coordinate transformation:
	// 		x = t/(1-t^2) 	<=> (1-t^2)x = t
	// 										<=> t^2*x + t - x = 0
	// 										<=> t^2 + t/x - 1 = 0
	// 										<=> (t + 1/(2x))^2 - 1/(4x^2) -1 = 0
	// 										<=> t = -1/(2x) +- sqrt(1 + 1/(4x^2))
	// 		x -> +-inf => t -> +- 1
	// 		dx/dt 	= (1*(1-t^2) - t(-2*t))/(1-t^2)^2
	// 						= (1-t^2 + 2*t^2)/(1-t^2)^2
	// 						= (1+t^2)/(1-t^2)^2
	// NOTE: to handle -inf -> inf and inf -> -inf, sign of endpoints need
	// to be checked.
	f64 t2 = t*t;
	f64 one_minus_t2 = (1.0-t2);
	return tdata->original_f(t/one_minus_t2, tdata->userdata) * (1.0+t2)/(one_minus_t2*one_minus_t2);
}
static inline f64 integrand_end_endpoint_inf(f64 t, void* data) {
	coordinate_transform_data* tdata = data;
	// x = a + t/(1-t) 	<=> (x-a)(1-t) = t
	// 									<=> (x-a) = t + t(x-a)
	// 									<=> (x-a) = t(1 + (x-a))
	// 									<=> t = (x-a)/(1 + (x-a))
	// dx/dt 	= d/dt (t/(1-t))
	// 				= (1*(1-t) - t*(-1))/(1-t)^2
	// 				= 1/(1-t)^2
	f64 one_minus_t = (1.0 - t);
	return tdata->original_f(tdata->original_start + t/one_minus_t, tdata->userdata) * (1.0/(one_minus_t*one_minus_t));
}
static inline f64 integrand_start_endpoint_inf(f64 t, void* data) {
	coordinate_transform_data* tdata = data;
	// x = a - t/(1-t)	<=> (x-a)(1-t) = -t
	// 									<=> (x-a) = t(x-a) - t
	// 									<=> (x-a) = t((x-a) - 1)
	// 									<=> t = (x-a)/((x-a) - 1)
	// dx/dt 	= - d/dt (t/(1-t))
	// 				= - 1/(1-t)^2
	// NOTE: when integrating -inf -> b, switch to b -> -inf and then apply
	// transformation to get 0->1
	f64 one_minus_t = (1.0 - t);
	return tdata->original_f(tdata->original_start - t/one_minus_t, tdata->userdata) * (-1.0/(one_minus_t*one_minus_t));
}

static inline integration_result hadapt(integrand* f, f64 start, f64 end, integration_settings settings) {
	integration_result result = {
		.integral = 0,
		.error = 0,
		.performed_evals = 0,
		.converged = false,
	};

	// Perform one evaluation of the rule and see if we can
	// exit early without having to initalize extra memory and
	// so on.
	segment s = evaluate_rule(f, start, end, settings);
	result.performed_evals += 4*settings.order + 2;
	result.integral = s.integral;
	result.error = s.error;

	if (should_exit(result, settings)) {
		if (result.performed_evals <= settings.max_evals)
			result.converged = true;
		return result;
	}

	// If we get to this point we were not able to exit early and have do to more subdivisions
	// to get desired error tolerance. yay.
	prioqueue* pq = prioqueue_new(2*settings.max_evals, sizeof(segment), compare_segments);
	prioqueue_push(pq, &s);

	while (!should_exit(result, settings)) {
		segment* largest_error_seg = (segment*)prioqueue_top(pq);
		prioqueue_pop(pq);

		f64 midpoint = 0.5 * (largest_error_seg->start + largest_error_seg->end);
		segment left_seg  = evaluate_rule(f, largest_error_seg->start, midpoint, settings);
		segment right_seg = evaluate_rule(f, midpoint, largest_error_seg->end, settings);

		result.performed_evals += 4*settings.order+2;
		result.integral = (result.integral - largest_error_seg->integral) + left_seg.integral + right_seg.integral;
		result.error = (result.error - largest_error_seg->error) + left_seg.error + right_seg.error;

		prioqueue_push(pq, &left_seg);
		prioqueue_push(pq, &right_seg);
	}

	prioqueue_free(pq);

	if (result.performed_evals <= settings.max_evals)
		result.converged = true;

	return result;
}

integration_result quadgk(integrand* f, f64 start, f64 end, integration_settings settings) {
	// Check for infinities in integration interval
	i32 start_isinf = isinf(start);
	i32 end_isinf = isinf(end);
	if (start_isinf || end_isinf) {
		coordinate_transform_data tdata = {
			.original_f = f,
			.original_start = start,
			.original_end = end,
			.userdata = settings.userdata,
		};

		settings.userdata = &tdata;

		if (start_isinf && end_isinf) {
			return hadapt(integrand_both_endpoints_inf, copysign(-1, start), copysign(1, end), settings);
		} else if (end_isinf) {
			// endpoint infinity
			return hadapt(integrand_end_endpoint_inf, 0, copysign(1, end), settings);
		} else /*if (start_isinf) */ {
			// startoint infinity
			return hadapt(integrand_start_endpoint_inf, 0, copysign(1, end), settings);
		}
	} else {
		return hadapt(f, start, end, settings);
	}
}
