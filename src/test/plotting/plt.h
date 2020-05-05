#pragma once

#include <sbmf/common/common.h>

struct PlotState;
typedef struct PlotState PlotState;

extern PlotState* plt_init();
extern void plt_shutdown(PlotState* state);
extern void plt_update(PlotState* state);
extern void plt_wait_on_join(PlotState* state);
extern void plt_toggle_active(PlotState* state, i32 idx);
//extern void plt_update_until_closed(PlotState* state);

extern void plt_clear(PlotState* state);
extern void plt_1d(PlotState* state, f32* x, f32* y, u32 len);
extern void plt_2d(PlotState* state, f32* x, f32* y, f32* z, u32 len);
