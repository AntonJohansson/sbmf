#pragma once

struct PlotState;
typedef struct PlotState PlotState;

extern PlotState* plt_init();
extern void plt_shutdown(PlotState* state);
extern void plt_update(PlotState* state);
extern void plt_wait_on_join(PlotState* state);
//extern void plt_update_until_closed(PlotState* state);

extern void plt_clear(PlotState* state);
extern void plt_1d(PlotState* state, unsigned int id, float* x, float* y, unsigned int len);
extern void plt_2d(PlotState* state, unsigned int id, float* points, unsigned int len);
extern void plt_3d(PlotState* state, unsigned int id, float* points, unsigned int len);
