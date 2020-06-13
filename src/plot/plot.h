#pragma once

#include "common.h"

struct plotstate;
typedef struct plotstate plotstate;

extern plotstate* plot_init();
extern void plot_shutdown(plotstate* state);
extern void plot_update(plotstate* state);
extern void plot_wait_on_join(plotstate* state);
extern void plot_toggle_active(plotstate* state, i32 idx);
//extern void plt_update_until_closed(plotstate* state);

extern void plot_clear(plotstate* state);
extern void plot_1d(plotstate* state, f32* x, f32* y, u32 len);
extern void plot_2d(plotstate* state, f32* x, f32* y, f32* z, u32 len);

typedef struct {
	u32 size;
	f32* vbo_mem;
	u32* ebo_mem;
} plotsurface;

extern plotsurface* plot_make_surface(plotstate* state, u32 size);
extern void plot_free_surface(plotstate* state, plotsurface* surf);
extern void plot_surface(plotstate* state, plotsurface* surf, f32* z);
