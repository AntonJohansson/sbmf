#pragma once

#include "common.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

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

	void* memory;
	f64* lens;
	f64* mins;
	f64* maxs;
	f64* deltas;
	i32* pointcounts;

	f64* points;
} grid;

static inline grid generate_grid(int_t dimensions, real_t mins[], real_t maxs[], int_t pointcounts[]) {
	int_t total_pointcount = 1;
	for (int_t i = 0; i < dimensions; ++i)
		total_pointcount *= pointcounts[i];

	int_t size_points  = dimensions*total_pointcount*sizeof(real_t);

	void* mem = malloc(dimensions*(4*sizeof(real_t) + sizeof(int_t)) + size_points);
	grid g = {
		.total_pointcount = total_pointcount,
		.dimensions = dimensions,
		.memory = mem,
		.lens					= (real_t*)(mem),
		.mins 				= (real_t*)((char*)mem + 		dimensions*sizeof(real_t)),
		.maxs 				= (real_t*)((char*)mem +  2*dimensions*sizeof(real_t)),
		.deltas 			= (real_t*)((char*)mem +  3*dimensions*sizeof(real_t)),
		.pointcounts 	=  (int_t*)((char*)mem + 	4*dimensions*sizeof(real_t)),
		.points 			= (real_t*)((char*)mem +  4*dimensions*sizeof(real_t) + dimensions*sizeof(int_t))
	};

	for (int_t i = 0; i < dimensions; ++i) {
		g.lens[i] = maxs[i] - mins[i];
		g.mins[i] = mins[i];
		g.maxs[i] = maxs[i];
		g.pointcounts[i] = pointcounts[i];
		g.deltas[i] = g.lens[i]/pointcounts[i];
	}

	int_t indices[dimensions];
	memset(indices, 0, dimensions*sizeof(int_t));

	//for (int_t i = 0; i < total_pointcount; ++i) {
	//	for (int_t j = 0; j < dimensions; ++j) {
	//		g.points[dimensions*i + j] = mins[j] + indices[j]*g.deltas[j];
	//	}

	//	indices[dimensions-1] = (indices[dimensions-1]+1) % pointcounts[dimensions-1];
	//	for (int_t j = dimensions-1; j > 0; --j) {
	//		if (indices[j] % pointcounts[j] == 0)
	//			indices[j-1] = (indices[j-1]+1) % pointcounts[j-1];
	//		else
	//			break;
	//	}
	//}

	{
		int_t indices[g.dimensions];
		memset(indices, 0, g.dimensions*sizeof(int_t));
		for (int_t index = 0; index < total_pointcount; ++index) {
			// k = (index / (1))   		% l1
			// j = (index / (l1))  		% l2
			// i = (index / (l1*l2))  % l3
			// ...

			int_t prodlen = 1;
			for (int_t n = g.dimensions-1; n >= 0; --n) {
				indices[n] = fmod((index / prodlen), g.pointcounts[n]);
				prodlen *= g.pointcounts[n];
			}

			for (int_t n = 0; n < g.dimensions; ++n) {
				g.points[dimensions*index + n] = g.mins[n] + indices[n]*g.deltas[n];
				//printf("%d\t", indices[n]);
			}
			//printf("\n");
		}
	}

	return g;
}

static inline grid mimic_grid(grid base) {
	int_t size_points  = base.dimensions*base.total_pointcount*sizeof(real_t);
	int_t total_size = base.dimensions*(4*sizeof(real_t) + sizeof(int_t)) + size_points;

	void* mem = malloc(total_size);
	grid g = {
		.dimensions = base.dimensions,
		.total_pointcount = base.total_pointcount,
		.memory = mem,
		.lens					= (real_t*)(mem),
		.mins 				= (real_t*)((char*)mem +   base.dimensions*sizeof(real_t)),
		.maxs 				= (real_t*)((char*)mem + 2*base.dimensions*sizeof(real_t)),
		.deltas 			= (real_t*)((char*)mem + 3*base.dimensions*2*sizeof(real_t)),
		.pointcounts 	=  (int_t*)((char*)mem + 4*base.dimensions*sizeof(real_t)),
		.points 			= (real_t*)((char*)mem + 4*base.dimensions*sizeof(real_t) + base.dimensions*sizeof(int_t)),
	};

	memcpy(g.memory, base.memory, total_size);

	return g;
}

static inline void free_grid(grid g) {
	free(g.memory);
}

#define FOREACH_GRIDPOINT(grid, i) \
		for (int_t i = 0; i < grid.total_pointcount; i += grid.dimensions)
