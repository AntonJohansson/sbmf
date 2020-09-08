#include "quadgk_vec.h"

#include <sbmf/memory/prioqueue.h>
#include <sbmf/debug/log.h>
#include <math.h> // INFINITY and isinf

typedef f64 coord_transform_input_func(f64 x, f64 a, f64 b);
typedef f64 coord_transform_output_func(f64 x);

struct coordinate_transform {
	coord_transform_input_func* input;
	coord_transform_output_func* output;
	f64 original_start;
	f64 original_end;
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

static inline void calculate_sample_points_for_interval(integration_settings settings, u32 size, f64 out[size], f64 start, f64 end) {
	struct gk_data* gk = &settings.gk;
	f64 mid = 0.5 * (end - start);

	u32 reduced_size = (size-3)/4;
	for (u32 i = 0; i < reduced_size; ++i) {
		f64 xg = gk->kronod_nodes[2*i+1];
		f64 xk = gk->kronod_nodes[2*i];

		out[4*i + 0] = start + (1-xg)*mid;
		out[4*i + 1] = start + (1+xg)*mid;
		out[4*i + 2] = start + (1-xk)*mid;
		out[4*i + 3] = start + (1+xk)*mid;
	}
	out[4*reduced_size + 0] = start + mid;
	out[4*reduced_size + 1] = start + (1-gk->kronod_nodes[gk->kronod_size-2])*mid;
	out[4*reduced_size + 2] = start + (1+gk->kronod_nodes[gk->kronod_size-2])*mid;
};

static eval_result evaluate_rule(
		integration_settings settings,
		u32 sample_count,
		f64 sampled_func_vals[sample_count],
		f64 start, f64 end,
		u32 is_order_odd
		) {
	struct gk_data* gk = &settings.gk;

	f64 mid = 0.5*(end - start);

	// fg/k - gauss/kronod function eval
	// Ig/k - gauss/kronod integral estimate
	// xg/k - gauss/kronod sample nodes

	u32 reduced_size = (sample_count - 3)/4;

	/* Compute the integral estimate */
	f64 Ig = 0, Ik = 0;
	for (u32 i = 0; i < reduced_size; ++i) {
		f64 fg = sampled_func_vals[4*i + 0] + sampled_func_vals[4*i + 1];
		f64 fk = sampled_func_vals[4*i + 2] + sampled_func_vals[4*i + 3];

		Ig += gk->gauss_weights[i] * fg;
		Ik += gk->kronod_weights[2*i+1] * fg + gk->kronod_weights[2*i] * fk;
	}

	/*
	 * In the odd-order case, the last point has to be handled
	 * differently
	 * */
	if (is_order_odd  == 0) {
		Ik += gk->kronod_weights[gk->kronod_size-1] * sampled_func_vals[4*reduced_size + 0];
	} else {
		Ig += gk->gauss_weights[gk->gauss_size-1] * sampled_func_vals[4*reduced_size + 0];
		Ik += gk->kronod_weights[gk->kronod_size-1] * sampled_func_vals[4*reduced_size + 0]
			+ gk->kronod_weights[gk->kronod_size-2] * (sampled_func_vals[4*reduced_size + 1] + sampled_func_vals[4*reduced_size + 2]);
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

/* Sort segments by error in descending order. */
static bool compare_segments(void* a, void* b) {
	return (((segment*)a)->error > ((segment*)b)->error);
}

static inline bool should_exit(integration_result result, integration_settings settings) {
	return (result.error <= settings.abs_error_tol ||
			result.error <= settings.rel_error_tol*norm(result.integral) ||
			result.performed_evals >= settings.max_evals);
}


/*
 * f(x)
 * x = g(t)
 * dx/dt = g'(t)
 * => f(x)*dx = f(g(t))*g'(t)*dt
 */


static f64 coord_transform_input_identity(f64 t, f64 a, f64 b) {
	return t;
}
static f64 coord_transform_output_identity(f64 t) {
	return 1.0;
}

static f64 coord_transform_input_both_inf(f64 t, f64 a, f64 b) {
	return t/(1.0 - t*t);
}
static f64 coord_transform_output_both_inf(f64 t) {
	f64 one_minus_t2 = 1.0 - t*t;
	return (1.0 + t*t)/(one_minus_t2*one_minus_t2);
}

static f64 coord_transform_input_start_inf(f64 t, f64 a, f64 b) {
	return b - t/(1.0 - t);
}
static f64 coord_transform_output_start_inf(f64 t) {
	f64 one_minus_t = 1.0 - t;
	return -1.0/(one_minus_t*one_minus_t);
}

static f64 coord_transform_input_end_inf(f64 t, f64 a, f64 b) {
	return a + t/(1.0 - t);
}
static f64 coord_transform_output_end_inf(f64 t) {
	f64 one_minus_t = 1.0 - t;
	return 1.0/(one_minus_t*one_minus_t);
}





















static inline integration_result hadapt(integrand_vec* f, f64 start, f64 end,
		struct coordinate_transform transform,
		integration_settings settings) {
	integration_result result = {
		.integral = 0,
		.error = 0,
		.performed_evals = 0,
		.converged = false,
	};

	u32 is_order_odd = 1 - (settings.gk.kronod_size & 1);
	u32 sample_size = 4*(settings.gk.gauss_size - is_order_odd) + 3;
	f64 sample_points[2*sample_size];
	f64 transformed_sample_points[2*sample_size];
	f64 sampled_func_vals[2*sample_size];

	calculate_sample_points_for_interval(settings, sample_size, sample_points, start, end);

	{
		/* Apply transform scaling to input */
		for (u32 i = 0; i < sample_size; ++i) {
			transformed_sample_points[i] = transform.input(sample_points[i], transform.original_start, transform.original_end);
		}

		/* Sample the function */
		f(sampled_func_vals, transformed_sample_points, sample_size, settings.userdata);

		/* Apply transform output scaling */
		for (u32 i = 0; i < sample_size; ++i) {
			sampled_func_vals[i] *= transform.output(sample_points[i]);
		}
	}

	/*
	 * Perform one evaluation of the rule and see if we can
	 * exit early without having to initalize extra memory and
	 * so on.
	 */
	//eval_result eval_res = evaluate_rule(f, start, end, transform, settings);
	eval_result eval_res = evaluate_rule(settings, sample_size, sampled_func_vals, start, end, is_order_odd);
	result.performed_evals += eval_res.func_evals;

	segment s = eval_res.seg;
	result.integral = s.integral;
	result.error = s.error;

	if (!eval_res.valid || should_exit(result, settings)) {
		if (result.performed_evals <= settings.max_evals)
			result.converged = true;
		return result;
	}

	/*
	 * If we get to this point we were not able to exit early and have do to more subdivisions
	 * to get desired error tolerance.
	 */

	u32 memory_marker = sbmf_stack_marker();

	prioqueue* pq = prioqueue_new(2*settings.max_evals, sizeof(segment), compare_segments);
	prioqueue_push(pq, &s);

	while (!should_exit(result, settings)) {
		segment* largest_error_seg = (segment*)prioqueue_top(pq);
		prioqueue_pop(pq);

		f64 midpoint = 0.5 * (largest_error_seg->start + largest_error_seg->end);
		calculate_sample_points_for_interval(settings, sample_size, sample_points, largest_error_seg->start, midpoint);
		calculate_sample_points_for_interval(settings, sample_size, &sample_points[sample_size], midpoint, largest_error_seg->end);

		{
			/* Apply transform scaling to input */
			for (u32 i = 0; i < 2*sample_size; ++i) {
				transformed_sample_points[i] = transform.input(sample_points[i], transform.original_start, transform.original_end);
			}

			/* Sample the function */
			f(sampled_func_vals, transformed_sample_points, 2*sample_size, settings.userdata);

			/* Apply transform output scaling */
			for (u32 i = 0; i < 2*sample_size; ++i) {
				sampled_func_vals[i] *= transform.output(sample_points[i]);
			}
		}

		eval_result left_eval_res = evaluate_rule(settings, sample_size, sampled_func_vals, largest_error_seg->start, midpoint, is_order_odd);
		eval_result right_eval_res = evaluate_rule(settings, sample_size, &sampled_func_vals[sample_size], midpoint, largest_error_seg->end, is_order_odd);
		if (!left_eval_res.valid || !right_eval_res.valid)
			return result;

		segment left_seg = left_eval_res.seg;
		segment right_seg = right_eval_res.seg;

		result.performed_evals += left_eval_res.func_evals + right_eval_res.func_evals;
		result.integral = (result.integral - largest_error_seg->integral) + left_seg.integral + right_seg.integral;
		result.error = (result.error - largest_error_seg->error) + left_seg.error + right_seg.error;

		prioqueue_push(pq, &left_seg);
		prioqueue_push(pq, &right_seg);
	}

	sbmf_stack_free_to_marker(memory_marker);

	if (result.performed_evals <= settings.max_evals)
		result.converged = true;

	return result;
}

integration_result quadgk_vec_inl(integrand_vec* f, f64 start, f64 end, integration_settings settings) {
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
	struct coordinate_transform transform = {
		.input = coord_transform_input_identity,
		.output = coord_transform_output_identity,
		.original_start = start,
		.original_end = end,
	};

	if (start_isinf || end_isinf) {
		if (start_isinf && end_isinf) {
			transform.input  = coord_transform_input_both_inf;
			transform.output = coord_transform_output_both_inf;
			res = hadapt(f, -1, 1, transform, settings);
		} else if (end_isinf) {
			transform.input  = coord_transform_input_end_inf;
			transform.output = coord_transform_output_end_inf;
			res = hadapt(f, 0, 1, transform, settings);
		} else if (start_isinf) {
			transform.input  = coord_transform_input_start_inf;
			transform.output = coord_transform_output_start_inf;
			res = hadapt(f, 1, 0, transform, settings);
		}
	} else {
		res = hadapt(f, start, end, transform, settings);
	}

	res.integral = output_factor*res.integral;
	return res;
}
