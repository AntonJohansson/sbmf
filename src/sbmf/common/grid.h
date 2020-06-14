#pragma once


#include <sbmf/common/common.h>
#include <sbmf/memory/stack_allocator.h>

#include <string.h> // memcpy, memset

// 2D
//  0  1  2  3
//  4  5  6  7
//  8  9 10 11
// 12 13 14 15

// (0,0)  (0,1)  (0,2)  (0,3)
// (1,0)  (1,1)  (1,2)  (1,3)
// (2,0)  (2,1)  (2,2)  (2,3)
// (3,0)  (3,1)  (3,2)  (3,3)
//
// (i,j) = (, index%4)

// If we have a 30x20x10 array the corresponding index would be for (i,j,k) would be
//
// 		index = i*(20*10) + j*10 + k
//
// we also have that
//
// 		k = index % 10
// 		j = (index/10) % 20
// 		i = (index/(20*10)) % 30
//
// however
//
// 		index = q * d + r <=> q = (index - r)/d = (index - index % d)/d
//
// i.e for d = 10
//
// 		q = (index - index % 10)/10
//
// all in all
//
// 		k = ((index - 0)/1)  % 10
// 		j = ((index - k)/10) % 20
// 		i = ((index - j)/20) % 30

// For example, consider the block of lengths 1x2x3 sampled with 128,64,32 points in each dimension, this
// corresponds to the grid with
//
//  	dimension = 3,
//  	lengths = {1, 2, 3},
//  	points_per_length = {128,64,32},
//  	deltas = {1/128, 2/64, 3/32},
//  	points = {...} (128 x 64 x 32 array)

typedef struct {
	i32 dimensions;
	i32 total_pointcount;

	f64* lens;
	f64* mins;
	f64* maxs;
	f64* deltas;
	i32* pointcounts;

	f64* points;
} grid;

static inline grid generate_grid(i32 dimensions, f64 mins[], f64 maxs[], i32 pointcounts[]) {
	i32 total_pointcount = 1;
	for (i32 i = 0; i < dimensions; ++i)
		total_pointcount *= pointcounts[i];

	grid g = {
		.dimensions = dimensions,
		.total_pointcount = total_pointcount,
		.lens        = (f64*)sa_push(_sbmf.main_stack, dimensions*sizeof(f64)),
		.mins        = (f64*)sa_push(_sbmf.main_stack, dimensions*sizeof(f64)),
		.maxs        = (f64*)sa_push(_sbmf.main_stack, dimensions*sizeof(f64)),
		.deltas      = (f64*)sa_push(_sbmf.main_stack, dimensions*sizeof(f64)),
		.pointcounts = (i32*)sa_push(_sbmf.main_stack, dimensions*sizeof(i32)),
		.points      = (f64*)sa_push(_sbmf.main_stack, dimensions*total_pointcount*sizeof(f64)),
	};

	for (i32 i = 0; i < dimensions; ++i) {
		g.lens[i] = maxs[i] - mins[i];
		g.mins[i] = mins[i];
		g.maxs[i] = maxs[i];
		g.pointcounts[i] = pointcounts[i];
		g.deltas[i] = g.lens[i]/pointcounts[i];
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
			for (i32 n = g.dimensions-1; n >= 0; --n) {
				indices[n] = fmod((index / prodlen), g.pointcounts[n]);
				prodlen *= g.pointcounts[n];
			}

			for (i32 n = 0; n < g.dimensions; ++n) {
				g.points[(g.dimensions)*index + n] = g.mins[n] + indices[n]*g.deltas[n];
				//printf("%d\t", indices[n]);
			}
			//printf("\n");
		}
	}

	return g;
}

static inline grid mimic_grid(grid base) {

	grid g = {
		.dimensions = base.dimensions,
		.total_pointcount = base.total_pointcount,
		.lens        = (f64*)sa_push(_sbmf.main_stack, base.dimensions*sizeof(f64)),
		.mins        = (f64*)sa_push(_sbmf.main_stack, base.dimensions*sizeof(f64)),
		.maxs        = (f64*)sa_push(_sbmf.main_stack, base.dimensions*sizeof(f64)),
		.deltas      = (f64*)sa_push(_sbmf.main_stack, base.dimensions*sizeof(f64)),
		.pointcounts = (i32*)sa_push(_sbmf.main_stack, base.dimensions*sizeof(i32)),
		.points      = (f64*)sa_push(_sbmf.main_stack, base.dimensions*base.total_pointcount*sizeof(f64)),
	};

	memcpy(g.lens, base.lens, base.dimensions*(4*sizeof(f64) + sizeof(i32)) + base.dimensions*base.total_pointcount*sizeof(f64));

	return g;
}
