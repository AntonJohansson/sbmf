static void bestmf_operator(const u32 len, f64* out, f64* in, const u32 component_count, f64* wf, void* userdata) {
	SBMF_UNUSED(in);
	for (u32 i = 0; i < len; ++i) {
		out[i] = 0.0;
		for (u32 j = 0; j < component_count; ++j) {
			f64* g = userdata;
			f64 c = fabs(wf[j*len + i]);
			out[i] += g[j]*c*c;
		}
	}
}












struct inner_product_integrand_params {
    f64* coeff_a;
    f64* coeff_b;
    u32  coeff_count;
    struct basis basis;
};

void inner_product_integrand(f64* out, f64* in, u32 len, void* data) {
    struct inner_product_integrand_params* p = data;

    f64 sample_a[len];
    p->basis.sample(p->coeff_count, p->coeff_a, len, sample_a, in);

    f64 sample_b[len];
    p->basis.sample(p->coeff_count, p->coeff_b, len, sample_b, in);

    for (u32 i = 0; i < len; ++i) {
        out[i] = sample_a[i]*sample_b[i];
    }
}

static void compute_occupations(struct nlse_settings settings, const i64 particle_count, const u32 coeff_count, f64* coeff, i64* n1, i64* n2) {
	/* Calculates innner product to get n1,n2 */
	sbmf_log_info("best_meanfield: finding occupations");
	{
		struct inner_product_integrand_params p = {
			.coeff_a = &coeff[0],
			.coeff_b = &coeff[coeff_count],
			.coeff_count = coeff_count,
			.basis = settings.basis,
		};

		struct quadgk_settings int_settings = {
			.gk = gk15,
			.abs_error_tol = 1e-10,
			.max_iters = settings.max_quadgk_iters,
			.userdata = &p,
		};

		u8 quadgk_memory[quadgk_required_memory_size(&int_settings)];

		struct quadgk_result ires;
		quadgk_infinite_interval(inner_product_integrand, &int_settings, quadgk_memory, &ires);

		/*
		 * <p1|p2> = (n2-n1)/N = (N-n1-n1)/N
		 *  => n1 = (N-N<p1|p2>)/2
		 */
		*n1 = lround(particle_count * (1.0 - ires.integral) / 2.0);
		*n2 = lround(particle_count - *n1);
	}
}


static void compute_wavefunctions(struct nlse_settings settings, const u32 particle_count, const u32 coeff_count, f64* coeff, const f64 n1, const f64 n2, f64* out) {
		SBMF_UNUSED(settings);
		/* P1 = sqrt(n1/N)p1 + sqrt(n2/N)p2
		 * P2 = sqrt(n2/N)p2 - sqrt(n1/N)p1
		 *
		 * sqrt(n1/N)p1 = P1 - sqrt(n2/N)p2
		 * P2 = sqrt(n2/N)p2 - P1 + sqrt(n2/N)p2
		 * p2 = 0.5*(P1+P2)*sqrt(N/n2)
		 *
		 * sqrt(n2/N)p2 = P2 + sqrt(n1/N)p1
		 * P1 = sqrt(n1/N)p1 + P2 + sqrt(n1/N)p1
		 * p1 = 0.5*(P1-P2)*sqrt(N/n1)
		 */

		f64 scaling_p1 = 0.5*sqrt((f64)particle_count/n1);
		f64 scaling_p2 = 0.5*sqrt((f64)particle_count/n2);
		for (u32 i = 0; i < coeff_count; ++i) {
			if (!isinf(scaling_p1) && !isnan(scaling_p1))
				out[i]               = scaling_p1 * (coeff[i] - coeff[coeff_count+i]);
			if (!isinf(scaling_p2) && !isnan(scaling_p2))
				out[coeff_count + i] = scaling_p2 * (coeff[i] + coeff[coeff_count+i]);
		}

		//f64_normalize(out, out, coeff_count);
		//f64_normalize(&out[coeff_count], &out[coeff_count], coeff_count);
}






//void ensure_structure(struct nlse_settings settings, struct nlse_result res) {
//	u32 particle_count = *(u32*)settings.post_normalize_userdata;
//
//	u32 n1, n2;
//	compute_occupations(settings, particle_count, res.coeff_count, res.coeff, &n1, &n2);
//
//	f64 p[2*res.coeff_count];
//	compute_wavefunctions(settings, particle_count, res.coeff_count, res.coeff, n1, n2, p);
//
//	for (u32 i = 0; i < res.coeff_count; ++i) {
//		res.coeff[i] =
//			  sqrt(n1/particle_count)*p[i]
//			+ sqrt(n2/particle_count)*p[res.coeff_count + i];
//		res.coeff[res.coeff_count + i] =
//			  sqrt(n2/particle_count)*p[res.coeff_count + i]
//			- sqrt(n1/particle_count)*p[i];
//	}
//
//	f64_normalize(res.coeff, res.coeff, res.coeff_count);
//	f64_normalize(&res.coeff[res.coeff_count], &res.coeff[res.coeff_count], res.coeff_count);
//}






f64 best_meanfield_energy(struct nlse_settings settings, const u32 coeff_count, f64* coeff, i64 n1, i64 n2, f64 g0) {
	return grosspitaevskii_energy(settings, coeff_count, 2, coeff, (i64[]){n1, n2}, (f64[]){g0,g0,g0,g0});
}








f64 bestmf_find_fractional_occupation(struct nlse_settings settings, const i64 particle_count, f64 g0, struct nlse_guess* guesses) {
	settings.post_normalize_userdata = &particle_count;
	//settings.post_normalize_callback = ensure_structure;

	/* Construct new g matrix including occupations factors */
	f64 g[2*2] = {
		0.75 * g0 * (particle_count - 1), 0.25 * g0 * (particle_count - 1),
		0.25 * g0 * (particle_count - 1), 0.75 * g0 * (particle_count - 1),
	};

	struct nlse_component comps[2];
	for (u32 i = 0; i < 2; ++i) {
		if (guesses)
			comps[i].guess = guesses[i];
		else
			comps[i].guess.type = DEFAULT_GUESS;
		comps[i].op = bestmf_operator;
		comps[i].userdata = &g[i*2];
	}

	struct nlse_result res = nlse_solver(settings, 2, comps);

	f64 n1_over_N;
	{
		/* Computes inner product <1|2> */
		struct quadgk_result ires;
		{
			struct inner_product_integrand_params p = {
				.coeff_a = &res.coeff[0],
				.coeff_b = &res.coeff[res.coeff_count],
				.coeff_count = res.coeff_count,
				.basis = settings.basis,
			};

			struct quadgk_settings int_settings = {
				.gk = gk15,
				.abs_error_tol = 1e-10,
				.max_iters = settings.max_quadgk_iters,
				.userdata = &p,
			};

			u8 quadgk_memory[quadgk_required_memory_size(&int_settings)];

			quadgk_infinite_interval(inner_product_integrand, &int_settings, quadgk_memory, &ires);
		}

		/*
		 * <1|2> = (n2-n1)/N = (N-n1-n1)/N = 1 - 2n1/N
		 *  => n1/N = (1 - <1|2>)/2
		 */
		n1_over_N =  (1.0 - ires.integral)/2.0;
	}

	return n1_over_N;
}
