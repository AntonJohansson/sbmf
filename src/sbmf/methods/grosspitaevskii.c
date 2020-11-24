#include <sbmf/methods/grosspitaevskii.h>

void operator(const u32 len, f64 out[static len],
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

struct nlse_result grosspitaevskii(struct gp_settings settings) {
	struct nlse_settings nlse_settings = {
		.max_iterations = 1e7,
		.error_tol = 1e-10,

		.spatial_pot_perturbation = settings.pot,

		.gk = gk15,

		.num_basis_funcs = settings.num_basis_funcs,
		.basis = settings.basis,

		.zero_threshold = settings.zero_threshold,

		.debug_callback = settings.debug_callback,
		.measure_every = settings.measure_every,
	};


	/* Construct new g matrix including occupations factors */
	f64 g[settings.component_count * settings.component_count];
	for (u32 i = 0; i < settings.component_count; ++i) {
		for (u32 j = 0; j < settings.component_count; ++j) {
			f64 factor = 1.0;
			if (i == j) {
				factor = (settings.occupations[i]-1);
			} else {
				factor = settings.occupations[j];
			}
			g[i*settings.component_count + j] = factor*settings.g0[i*settings.component_count + j];
		}
	}

	struct nlse_component comps[settings.component_count];
	for (u32 i = 0; i < settings.component_count; ++i) {
		if (settings.guesses)
			comps[i].guess = settings.guesses[i];
		else
			comps[i].guess.type = DEFAULT_GUESS;
		comps[i].op = operator;
		comps[i].userdata = &g[i*settings.component_count];
	}

	struct nlse_result res = nlse_solver(nlse_settings, settings.component_count, comps);

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

f64 full_energy(struct gp_settings settings, struct nlse_result res) {
	f64 E = 0.0;

	for (u32 i = 0; i < settings.component_count; ++i) {
		for (u32 j = 0; j < res.coeff_count; ++j) {
			f64 c = fabs(res.coeff[i*res.coeff_count + j]);
			E += settings.occupations[i]*settings.basis.eigenval(j)*c*c;
		}
	}

	struct full_energy_integrand_params p = {
		.coeff_count = res.coeff_count,
		.basis = settings.basis,
	};

	struct full_energy_integrand_pot_params ppot = {
		.coeff_count = res.coeff_count,
		.V = settings.pot,
		.basis = settings.basis,
	};

	integration_settings int_settings = {
		.max_evals = 1e5,
		.abs_error_tol = 1e-10,
		.userdata = &ppot,
	};

	/* pot terms */
	if (settings.pot) {
		for (u32 i = 0; i < res.component_count; ++i) {
			ppot.coeff_a = &res.coeff[i*res.coeff_count];
			integration_result ires = quadgk_vec(full_energy_integrand_pot, -INFINITY, INFINITY, int_settings);
			E += settings.occupations[i]*ires.integral;
		}
	}

	/* |a|^2|b|^2 terms */
	int_settings.userdata = &p;
	for (u32 i = 0; i < res.component_count; ++i) {
		for (u32 j = 0; j < res.component_count; ++j) {
			f64 factor = 0.0;
			if (i == j)
				factor = 0.5*(settings.occupations[i]-1);
			else
				factor = settings.occupations[j];


			p.coeff_a = &res.coeff[i*res.coeff_count];
			p.coeff_b = &res.coeff[j*res.coeff_count];
			integration_result ires = quadgk_vec(full_energy_integrand, -INFINITY, INFINITY, int_settings);
			E += settings.g0[i*res.component_count + j] * settings.occupations[i] * factor * ires.integral;
		}
	}

	return E;
}
