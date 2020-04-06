#include "groundstate_solver.h"

#include <stdlib.h>

void gss_free_result(gss_result res) {
	free(res.wavefunction);
}
