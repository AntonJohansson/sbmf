#include "quadgk.h"
#include <sbmf/memory/prioqueue.h>
#include <sbmf/debug/log.h>
#include <sbmf/debug/profile.h>
#include <math.h> // INFINITY and isinf

struct gk_data gk7 = {
	.kronod_nodes = {
		9.9145537112081263920685469752598e-01,
		9.4910791234275852452618968404809e-01,
		8.6486442335976907278971278864098e-01,
		7.415311855993944398638647732811e-01,
		5.8608723546769113029414483825842e-01,
		4.0584515137739716690660641207707e-01,
		2.0778495500789846760068940377309e-01,
		0.0
	},
	.kronod_weights = {
		2.2935322010529224963732008059913e-02,
		6.3092092629978553290700663189093e-02,
		1.0479001032225018383987632254189e-01,
		1.4065325971552591874518959051021e-01,
		1.6900472663926790282658342659795e-01,
		1.9035057806478540991325640242055e-01,
		2.0443294007529889241416199923466e-01,
		2.0948214108472782801299917489173e-01
	},
	.gauss_weights = {
		1.2948496616886969327061143267787e-01,
		2.797053914892766679014677714229e-01,
		3.8183005050511894495036977548818e-01,
		4.1795918367346938775510204081658e-01
	},
	.kronod_size = 8,
	.gauss_size = 4,
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

typedef struct {
	segment seg;
	u32 func_evals;
	bool valid;
} eval_result;

eval_result evaluate_rule(integrand* f, f64 start, f64 end, integration_settings settings) {
	struct gk_data* gk = &settings.gk;
	f64 mid = 0.5 * (end - start);

	// 	0 -- if even size
	// 	1 -- if odd size
	i32 is_order_odd = 1 - (gk->kronod_size & 1);

	// fg/k - gauss/kronod function eval
	// Ig/k - gauss/kronod integral estimate
	// xg/k - gauss/kronod sample nodes
	f64 Ig = 0, Ik = 0;

	// Size of the gauss weights array (other than the last element if odd order)
	u32 size = gk->gauss_size - is_order_odd;
	f64 xvals[4*size + 3];
	for (u32 i = 0; i < size; ++i) {
		f64 xg = gk->kronod_nodes[2*i+1];
		f64 xk = gk->kronod_nodes[2*i];
		xvals[4*i + 0] = start + (1-xg)*mid;
		xvals[4*i + 1] = start + (1+xg)*mid;
		xvals[4*i + 2] = start + (1-xk)*mid;
		xvals[4*i + 3] = start + (1+xk)*mid;
	}
	xvals[4*size + 0] = start+mid;
	xvals[4*size + 1] = start + (1-gk->kronod_nodes[gk->kronod_size-2])*mid;
	xvals[4*size + 2] = start + (1+gk->kronod_nodes[gk->kronod_size-2])*mid;

	f64 out[4*size+3];
	f(out, xvals, 4*size+3, settings.userdata);

	for (u32 i = 0; i < size; ++i) {
		f64 fg = out[4*i + 0] + out[4*i + 1];
		f64 fk = out[4*i + 2] + out[4*i + 3];
		Ig += gk->gauss_weights[i] * fg;
		Ik += gk->kronod_weights[2*i+1] * fg + gk->kronod_weights[2*i] * fk;
	}

	// In the odd-order case, the last point has to be handled
	// differently
	if (is_order_odd == 0) {
		Ik += gk->kronod_weights[gk->kronod_size-1] * out[4*size + 0];
	} else {
		Ig += gk->gauss_weights[gk->gauss_size-1] * out[4*size + 0];
		Ik += gk->kronod_weights[gk->kronod_size-1] * out[4*size + 0]
			+ gk->kronod_weights[gk->kronod_size-2] * (out[4*size+1] + out[4*size+2]);
	}

	f64 Ik_s = Ik*mid;
	f64 Ig_s = Ig*mid;
	f64 error = norm(Ik_s - Ig_s);

	bool valid = (!isinf(error) && !isnan(error));
	if (!valid) {
		log_error("evaluation of gk rule resulted in [%s]", (isinf(error) == 1) ? "inf" : "nan");
	}

	return (eval_result){
		.seg = {Ik_s, error, start, end},
		.func_evals = 4*gk->gauss_size + 1,
		.valid = valid
	};
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
	//return tdata->original_f(t/one_minus_t2, tdata->userdata) * (1.0+t2)/(one_minus_t2*one_minus_t2);
	return 0;
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
	//return tdata->original_f(tdata->original_start + t/one_minus_t, tdata->userdata) * (1.0/(one_minus_t*one_minus_t));
	return 0;
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
	//return tdata->original_f(tdata->original_end - t/one_minus_t, tdata->userdata) * (-1.0/(one_minus_t*one_minus_t));
	return 0;
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
	eval_result eval_res = evaluate_rule(f, start, end, settings);
	result.performed_evals += eval_res.func_evals;

	segment s = eval_res.seg;
	result.integral = s.integral;
	result.error = s.error;

	if (!eval_res.valid || should_exit(result, settings)) {
		if (result.performed_evals <= settings.max_evals)
			result.converged = true;
		return result;
	}

	u32 memory_marker = sbmf_stack_marker();

	// If we get to this point we were not able to exit early and have do to more subdivisions
	// to get desired error tolerance.

	u32 memory_marker = sbmf_stack_marker();

	prioqueue* pq = prioqueue_new(2*settings.max_evals, sizeof(segment), compare_segments);
	prioqueue_push(pq, &s);

	while (!should_exit(result, settings)) {
		PROFILE_BEGIN("quadgk -- hadapt iter");
			PROFILE_BEGIN("quadgk -- hadapt pq");
				segment* largest_error_seg = (segment*)prioqueue_top(pq);
				prioqueue_pop(pq);
			PROFILE_END("quadgk -- hadapt pq");

			f64 midpoint = 0.5 * (largest_error_seg->start + largest_error_seg->end);
			PROFILE_BEGIN("quadgk -- hadapt eval");
				eval_result left_eval_res = evaluate_rule(f, largest_error_seg->start, midpoint, settings);
				eval_result right_eval_res = evaluate_rule(f, midpoint, largest_error_seg->end, settings);
			PROFILE_END("quadgk -- hadapt eval");
			if (!left_eval_res.valid || !right_eval_res.valid)
				return result;

			segment left_seg = left_eval_res.seg;
			segment right_seg = right_eval_res.seg;

			result.performed_evals += left_eval_res.func_evals + right_eval_res.func_evals;
			result.integral = (result.integral - largest_error_seg->integral) + left_seg.integral + right_seg.integral;
			result.error = (result.error - largest_error_seg->error) + left_seg.error + right_seg.error;

			prioqueue_push(pq, &left_seg);
			prioqueue_push(pq, &right_seg);
		PROFILE_END("quadgk -- hadapt iter");
	}

	sbmf_stack_free_to_marker(memory_marker);

	if (result.performed_evals <= settings.max_evals)
		result.converged = true;

	return result;
}

integration_result quadgk(integrand* f, f64 start, f64 end, integration_settings settings) {
	// Make sure the endpoints are ordered start < end
	f64 output_factor = 1.0;
	if (start > end) {
		f64 temp = start;
		start = end;
		end = temp;
		output_factor = -1.0;
	}

	// Check for infinities in integration interval
	i32 start_isinf = isinf(start);
	i32 end_isinf = isinf(end);

	integration_result res;

	if (start_isinf || end_isinf) {
		coordinate_transform_data tdata = {
			.original_f = f,
			.original_start = start,
			.original_end = end,
			.userdata = settings.userdata,
		};

		settings.userdata = &tdata;

		if (start_isinf && end_isinf) {
			//res = hadapt(integrand_both_endpoints_inf, -1, 1, settings);
		} else if (end_isinf) {
			//res = hadapt(integrand_end_endpoint_inf, 0, 1, settings);
		} else if (start_isinf) {
			//res = hadapt(integrand_start_endpoint_inf, 1, 0, settings);
		}
	} else {
		res = hadapt(f, start, end, settings);
	}

	res.integral = output_factor*res.integral;
	return res;
}
