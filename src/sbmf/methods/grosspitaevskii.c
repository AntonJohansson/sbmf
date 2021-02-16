static void gp_operator(const u32 len, f64* out, f64* in, const u32 component_count, f64* wf, void* userdata) {
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

struct nlse_result grosspitaevskii(struct nlse_settings settings, const u32 comp_count, i64* occupations, struct nlse_guess* guesses, f64* g0) {

	/* Construct new g matrix including occupations factors */
	f64 g[comp_count * comp_count];
	for (u32 i = 0; i < comp_count; ++i) {
		for (u32 j = 0; j < comp_count; ++j) {
			f64 factor = 1.0;
			if (i == j) {
				factor = (occupations[i]-1);
			} else {
				factor = occupations[j];
			}
			g[i*comp_count + j] = factor*g0[i*comp_count + j];
		}
	}

	struct nlse_component comps[comp_count];
	for (u32 i = 0; i < comp_count; ++i) {
		if (guesses)
			comps[i].guess = guesses[i];
		else
			comps[i].guess.type = DEFAULT_GUESS;
		comps[i].op = gp_operator;
		comps[i].userdata = &g[i*comp_count];
	}

	struct nlse_result res = nlse_solver(settings, comp_count, comps);

	return res;
}











struct full_energy_integrand_params {
	u32 coeff_count;
	f64* coeff_a;
	f64* coeff_b;
	struct basis basis;
};

struct full_energy_integrand_pot_params {
	u32 coeff_count;
	f64* coeff_a;
	nlse_operator_func* V;
	struct basis basis;
};

/* Integrand of the form a*V*a */
void full_energy_integrand_pot(f64* out, f64* in, u32 len, void* data) {
	struct full_energy_integrand_pot_params* p = data;

	f64 sample[len];
	p->basis.sample(p->coeff_count, p->coeff_a, len, sample, in);

	f64 pot[len];
	p->V(len, pot, in, 0, NULL, NULL);

	for (u32 i = 0; i < len; ++i) {
		out[i] = sample[i]*sample[i]*pot[i];
	}
}

/* Integrand of the form |a|^2|b|^2 */
void full_energy_integrand(f64* out, f64* in, u32 len, void* data) {
	struct full_energy_integrand_params* p = data;

	f64 sample_a[len];
	p->basis.sample(p->coeff_count, p->coeff_a, len, sample_a, in);

	f64 sample_b[len];
	p->basis.sample(p->coeff_count, p->coeff_b, len, sample_b, in);

	for (u32 i = 0; i < len; ++i) {
		out[i] = sample_a[i]*sample_a[i]*sample_b[i]*sample_b[i];
	}
}

f64 grosspitaevskii_energy(struct nlse_settings settings, const u32 coeff_count, const u32 comp_count, f64* coeff, i64* occupations, f64* g0) {
	f64 E = 0.0;

	for (u32 i = 0; i < comp_count;++i) {
		for (u32 j = 0; j < coeff_count; ++j) {
			f64 c = fabs(coeff[i*coeff_count + j]);
			E += occupations[i]*settings.basis.eigenval(j)*c*c;
		}
	}

	struct full_energy_integrand_params p = {
		.coeff_count = coeff_count,
		.basis = settings.basis,
	};

	struct full_energy_integrand_pot_params ppot = {
		.coeff_count = coeff_count,
		.V = settings.spatial_pot_perturbation,
		.basis = settings.basis,
	};

	struct quadgk_settings int_settings = {
		.max_iters = settings.max_quadgk_iters,
		.abs_error_tol = 1e-15,
		.userdata = &ppot,
		.gk = gk20
	};

	u8 quadgk_memory[quadgk_required_memory_size(&int_settings)];


	/* pot terms */
	if (settings.spatial_pot_perturbation) {
		struct quadgk_result ires;
		for (u32 i = 0; i < comp_count; ++i) {
			ppot.coeff_a = &coeff[i*coeff_count];
			quadgk_infinite_interval(full_energy_integrand_pot, &int_settings, quadgk_memory, &ires);
			E += occupations[i]*ires.integral;
		}
	}

	/* |a|^4 terms within comp */
	int_settings.userdata = &p;
	struct quadgk_result ires;
	for (u32 i = 0; i < comp_count; ++i) {
		p.coeff_a = &coeff[i*coeff_count];
		p.coeff_b = &coeff[i*coeff_count];
		quadgk_infinite_interval(full_energy_integrand, &int_settings, quadgk_memory, &ires);
		E += 0.5 * g0[i*comp_count + i] * occupations[i] * (occupations[i]-1) * ires.integral;
	}

	/* |a|^2|b|^2 terms between comps */
	for (u32 i = 0; i < comp_count; ++i) {
		for (u32 j = i+1; j < comp_count; ++j) {
			p.coeff_a = &coeff[i*coeff_count];
			p.coeff_b = &coeff[j*coeff_count];
			quadgk_infinite_interval(full_energy_integrand, &int_settings, quadgk_memory, &ires);
			E += g0[i*comp_count + j] * occupations[i] * occupations[j] * ires.integral;
		}
	}

	return E;
}
