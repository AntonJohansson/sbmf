static void bestmf_operator(const u32 len, f64 out[static len],
			  f64 in[static len], const u32 component_count,
			  f64 wf[static len*component_count],
			  void* userdata) {
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

static void compute_occupations(struct nlse_settings settings, const u32 particle_count,
		const u32 coeff_count, f64 coeff[static coeff_count],
		u32* n1, u32* n2) {
	/* Calculates innner product to get n1,n2 */
	sbmf_log_info("best_meanfield: finding occupations");
	{
		struct inner_product_integrand_params p = {
			.coeff_a = &coeff[0],
			.coeff_b = &coeff[coeff_count],
			.coeff_count = coeff_count,
			.basis = settings.basis,
		};

		integration_settings int_settings = {
			.gk = gk15,
			.abs_error_tol = 1e-10,
			.max_evals = 1e5,
			.userdata = &p,
		};

		integration_result ires = quadgk(inner_product_integrand, -INFINITY, INFINITY, int_settings);

		/*
		 * <p1|p2> = (n2-n1)/N = (N-n1-n1)/N
		 *  => n1 = (N-N<p1|p2>)/2
		 */
		*n1 = particle_count * (1.0 - ires.integral) / 2.0;
		*n2 = particle_count - *n1;
	}
}



static void compute_wavefunctions(struct nlse_settings settings, const u32 particle_count,
						   const u32 coeff_count, f64 coeff[static coeff_count],
						   const f64 n1, const f64 n2, f64* out) {
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






void ensure_structure(struct nlse_settings settings, struct nlse_result res) {
	u32 particle_count = *(u32*)settings.post_normalize_userdata;

	u32 n1, n2;
	compute_occupations(settings, particle_count, res.coeff_count, res.coeff, &n1, &n2);

	f64 p[2*res.coeff_count];
	compute_wavefunctions(settings, particle_count, res.coeff_count, res.coeff, n1, n2, p);

	for (u32 i = 0; i < res.coeff_count; ++i) {
		res.coeff[i] =
			  sqrt(n1/particle_count)*p[i]
			+ sqrt(n2/particle_count)*p[res.coeff_count + i];
		res.coeff[res.coeff_count + i] =
			  sqrt(n2/particle_count)*p[res.coeff_count + i]
			- sqrt(n1/particle_count)*p[i];
	}

	f64_normalize(res.coeff, res.coeff, res.coeff_count);
	f64_normalize(&res.coeff[res.coeff_count], &res.coeff[res.coeff_count], res.coeff_count);
}














struct bestmf_result best_meanfield(struct nlse_settings settings, const u32 particle_count, f64 g0, struct nlse_guess* guesses) {
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

	u32 n1, n2;
	compute_occupations(settings, particle_count, res.coeff_count, res.coeff, &n1, &n2);


	/* Find the wavefunctions */
	sbmf_log_info("best_meanfield: finding wavefunctions");
	f64* p = (f64*)sbmf_stack_push(2*res.coeff_count*sizeof(f64));
	memset(p, 0, 2*res.coeff_count*sizeof(f64));
	compute_wavefunctions(settings, particle_count, res.coeff_count, res.coeff, n1, n2, p);

	sbmf_log_info("best_meanfield: finding energy");
	f64 E = full_energy(settings,
			res.coeff_count, res.component_count,
			p, (u32[]){(u32)n1, particle_count - (u32)n1},
			(f64[]){g0,g0,g0,g0});

	return (struct bestmf_result) {
		.energy = E,
		.coeff_count = res.coeff_count,
		.comp_count = res.component_count,
		.coeff = p,
		.n1 = n1,
		.n2 = n2,
	};
}

struct bestmf_2comp_result best_meanfield_2comp(struct nlse_settings settings,
		const u32 particle_count,
		f64 g0[static 2*2],
		struct nlse_guess* guesses) {
	f64 g[4*4] = {
		0.75*g0[0]*(particle_count - 1), 0.25*g0[0]*(particle_count - 1), 0.50*g0[1]*particle_count,     0.50*g0[1]*particle_count,
		0.25*g0[0]*(particle_count - 1), 0.75*g0[0]*(particle_count - 1), 0.50*g0[1]*particle_count,     0.50*g0[1]*particle_count,
		0.50*g0[2]*particle_count,       0.50*g0[2]*particle_count,       0.75*g0[3]*(particle_count-1), 0.25*g0[3]*(particle_count-1),
		0.50*g0[2]*particle_count,       0.50*g0[2]*particle_count,       0.25*g0[3]*(particle_count-1), 0.75*g0[3]*(particle_count-1),
	};

	struct nlse_component comps[4];
	for (u32 i = 0; i < 4; ++i) {
		if (guesses)
			comps[i].guess = guesses[i];
		else
			comps[i].guess.type = DEFAULT_GUESS;
		comps[i].op = bestmf_operator;
		comps[i].userdata = &g[i*4];
	}

	struct nlse_result res = nlse_solver(settings, 4, comps);

	u32 n1, n2;
	compute_occupations(settings, particle_count, res.coeff_count, res.coeff, &n1, &n2);

	u32 n3, n4;
	compute_occupations(settings, particle_count, res.coeff_count, &res.coeff[2*res.coeff_count], &n3, &n4);

	f64* p = sbmf_stack_push(4*res.coeff_count*sizeof(f64));
	memset(p, 0, 4*res.coeff_count*sizeof(f64));
	compute_wavefunctions(settings, particle_count, res.coeff_count, res.coeff, n1, n2, p);
	compute_wavefunctions(settings, particle_count, res.coeff_count, &res.coeff[2*res.coeff_count], n3, n4, &p[2*res.coeff_count]);

	f64 E = full_energy(settings,
			res.coeff_count, res.component_count,
			p,
			(u32[]){
				(u32)n1, particle_count - (u32)n1,
				(u32)n3, particle_count - (u32)n3
				},
			(f64[]){
				g0[0],g0[0],g0[1],g0[1],
				g0[0],g0[0],g0[1],g0[1],
				g0[2],g0[2],g0[3],g0[3],
				g0[2],g0[2],g0[3],g0[3]
				});

	return (struct bestmf_2comp_result) {
		.energy = E,
		.coeff_count = res.coeff_count,
		.comp_count = res.component_count,
		.coeff = p,
		.n1 = n1,
		.n2 = n2,
		.n3 = n3,
		.n4 = n4,
	};
}
