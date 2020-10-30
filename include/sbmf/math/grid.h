#pragma once

#include <sbmf/sbmf.h>

#include <string.h> // memcpy, memset

struct grid {
	i32 dimensions;
	i32 total_pointcount;

	f64* lens;
	f64* mins;
	f64* maxs;
	f64* deltas;
	i32* pointcounts;

	f64* points;
};

static inline struct grid generate_grid(i32 dimensions, f64 mins[], f64 maxs[], i32 pointcounts[]) {
	i32 total_pointcount = 1;
	for (i32 i = 0; i < dimensions; ++i)
		total_pointcount *= pointcounts[i];

	struct grid g = {
		.dimensions = dimensions,
		.total_pointcount = total_pointcount,
		.lens        = (f64*)sbmf_stack_push(dimensions*sizeof(f64)),
		.mins        = (f64*)sbmf_stack_push(dimensions*sizeof(f64)),
		.maxs        = (f64*)sbmf_stack_push(dimensions*sizeof(f64)),
		.deltas      = (f64*)sbmf_stack_push(dimensions*sizeof(f64)),
		.pointcounts = (i32*)sbmf_stack_push(dimensions*sizeof(i32)),
		.points      = (f64*)sbmf_stack_push(dimensions*total_pointcount*sizeof(f64)),
	};

	for (i32 i = 0; i < dimensions; ++i) {
		g.lens[i] = maxs[i] - mins[i];
		g.mins[i] = mins[i];
		g.maxs[i] = maxs[i];
		g.deltas[i] = g.lens[i]/(pointcounts[i]-1);
		g.pointcounts[i] = pointcounts[i];
	}

	{
		i32 indices[g.dimensions];
		memset(indices, 0, g.dimensions*sizeof(i32));
		for (i32 index = 0; index < total_pointcount; ++index) {
			// k = (index / (1))   		% l1
			// j = (index / (l1))  		% l2
			// i = (index / (l1*l2))  % l3
			// ...

			i32 prodlen = 1;
			// TODO: does it matter if xyzxyzxyz... vs. zyxzyxzyx... ?
			for (i32 n = g.dimensions-1; n >= 0; --n) {
				indices[n] = fmod((index / prodlen), g.pointcounts[n]);
				indices[n] = (index/prodlen) % g.pointcounts[n];
				prodlen *= g.pointcounts[n];
			}

			for (i32 n = 0; n < g.dimensions; ++n) {
				g.points[g.dimensions*index + n] = g.mins[n] + indices[n]*g.deltas[n];
			}
		}
	}

	return g;
}

static inline struct grid mimic_grid(struct grid base) {

	struct grid g = {
		.dimensions = base.dimensions,
		.total_pointcount = base.total_pointcount,
		.lens        = (f64*)sbmf_stack_push(base.dimensions*sizeof(f64)),
		.mins        = (f64*)sbmf_stack_push(base.dimensions*sizeof(f64)),
		.maxs        = (f64*)sbmf_stack_push(base.dimensions*sizeof(f64)),
		.deltas      = (f64*)sbmf_stack_push(base.dimensions*sizeof(f64)),
		.pointcounts = (i32*)sbmf_stack_push(base.dimensions*sizeof(i32)),
		.points      = (f64*)sbmf_stack_push(base.dimensions*base.total_pointcount*sizeof(f64)),
	};

	memcpy(g.lens, base.lens, base.dimensions*(4*sizeof(f64) + sizeof(i32)) + base.dimensions*base.total_pointcount*sizeof(f64));

	return g;
}
