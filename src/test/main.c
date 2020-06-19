#define SNOW_ENABLED
#include "snow.h"

#include <sbmf/groundstate_solver/groundstate_solver.h>
#include <sbmf/common/fdm.h>
#include <sbmf/common/grid.h>
#include <sbmf/common/common.h>

#include <sbmf/common/eigenproblem.h>

#include <plot/plot.h>

#define PLOT_SPARSE_1D_HO_FDM 		0
#define PLOT_SPARSE_2D_HO_FDM 		0
#define PLOT_SPARSE_1D_PB_FDM 		0
#define PLOT_SPARSE_2D_PB_FDM 		0
#define PLOT_SPARSE_1D_PB_HOBASIS 	1

#include "basic.c"
#include "integration.c"
#include "fdm_solving.c"
#include "ho_basis_func.c"

snow_main();
