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

		.zero_threshold = 1e-10,
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
	struct nlse_result* res;
	struct gp_settings* settings;
};


void full_energy_integrand(f64* out, f64* in, u32 len, void* data) {
	struct full_energy_integrand_params* p = data;

	const u32 comp_count = p->settings->component_count;

	f64 pot[len];
	p->settings->pot(len, pot, in, 0, NULL, NULL);

	f64 samples[comp_count * len];
	for (u32 i = 0; i < comp_count; ++i) {
		p->settings->basis.sample(p->res->coeff_count, &p->res->coeff[i*p->res->coeff_count], len, &samples[i*len], in);
	}

	for (u32 i = 0; i < comp_count; ++i) {
		for (u32 j = 0; j < len; ++j) {
			f64 c = fabs(samples[i*len + j]);
			out[i] += p->settings->occupations[i]*pot[j]*c*c;
		}

		for (u32 j = 0; j < comp_count; ++j) {
			f64 factor;
			if (i == j)
				factor = 0.5*p->settings->occupations[i]*(p->settings->occupations[j]-1);
			else
				factor = p->settings->occupations[i]*p->settings->occupations[j];

			for (u32 k = 0; k < len; ++k) {
				f64 c = fabs(samples[j*len + k]);
				out[k] += p->settings->g0[i*comp_count + j] * factor * c*c*c*c;
			}
		}
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
		.settings = &settings,
		.res = &res,
	};

	integration_settings int_settings = {
		.max_evals = 1e5,
		.abs_error_tol = 1e-10,
		.userdata = &p,
	};

	integration_result ires = quadgk_vec(full_energy_integrand, -INFINITY, INFINITY, int_settings);
	E += ires.integral;

	return E;
}
