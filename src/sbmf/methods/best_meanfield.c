#include <sbmf/methods/best_meanfield.h>
#include <sbmf/math/functions.h>

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



void gaussian0(f64* out, f64* in, u32 len, void* data) {
    for (u32 i = 0; i < len; ++i)
        out[i] = gaussian(in[i] + 1.0, 0.0, 0.2);
}
void gaussian1(f64* out, f64* in, u32 len, void* data) {
    for (u32 i = 0; i < len; ++i)
        out[i] = gaussian(in[i] - 1.0, 0.0, 0.2);
}


struct nlse_result best_meanfield(struct nlse_settings settings, const u32 particle_count, f64 g0) {
	/* Construct new g matrix including occupations factors */
	f64 g[2*2] = {
		0.75 * g0 * (particle_count - 1), 0.25 * g0 * (particle_count - 1),
		0.25 * g0 * (particle_count - 1), 0.75 * g0 * (particle_count - 1),
	};

	struct nlse_component comps[2];
	for (u32 i = 0; i < 2; ++i) {
		comps[i].guess.type = SPATIAL_GUESS;
		comps[i].guess.data.spatial_guess = (i == 0) ? gaussian0 : gaussian1;
		comps[i].op = operator;
		comps[i].userdata = &g[i*2];
	}

	struct nlse_result res = nlse_solver(settings, 2, comps);

	return res;

}
