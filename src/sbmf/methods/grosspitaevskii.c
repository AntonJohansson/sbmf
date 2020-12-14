static void gp_operator(const u32 len, f64 out[static len],
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

struct nlse_result grosspitaevskii(struct nlse_settings settings,
		const u32 comp_count,
		u32 occupations[static comp_count],
		struct nlse_guess guesses[static comp_count],
		f64 g0[static comp_count*comp_count]) {

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

f64 full_energy(struct nlse_settings settings,
		const u32 coeff_count, const u32 comp_count,
		f64 coeff[static coeff_count*comp_count],
		u32 occupations[static comp_count],
		f64 g0[static comp_count*comp_count]
		) {
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

	integration_settings int_settings = {
		.max_evals = 1e5,
		.abs_error_tol = 1e-10,
		.userdata = &ppot,
		.gk = gk15
	};

	/* pot terms */
	if (settings.spatial_pot_perturbation) {
		for (u32 i = 0; i < comp_count; ++i) {
			ppot.coeff_a = &coeff[i*coeff_count];
			integration_result ires = quadgk(full_energy_integrand_pot, -INFINITY, INFINITY, int_settings);
			E += occupations[i]*ires.integral;
		}
	}

	/* |a|^2|b|^2 terms */
	int_settings.userdata = &p;
	for (u32 i = 0; i < comp_count; ++i) {
		for (u32 j = 0; j < comp_count; ++j) {
			f64 factor = 0.0;
			if (i == j)
				factor = 0.5*(occupations[i]-1);
			else
				factor = occupations[j];


			p.coeff_a = &coeff[i*coeff_count];
			p.coeff_b = &coeff[j*coeff_count];
			integration_result ires = quadgk(full_energy_integrand, -INFINITY, INFINITY, int_settings);
			E += g0[i*comp_count + j] * occupations[i] * factor * ires.integral;
		}
	}

	return E;
}
