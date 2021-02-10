#include <sbmf/sbmf.h>
#include <assert.h>

#define sbmf_stack_push(size_in_bytes) \
    sbmf_stack_push_impl(size_in_bytes, __LINE__, __FILE__, __func__)

__host__
void* sbmf_stack_push_impl(u32 size_in_bytes, const u32 linenumber, const char file[], const char func[]);

#include "functions_cuda.cu"
#include "basis_cuda.cu"
#include "indices_cuda.cu"
#include "hermite_integrals_cuda.cu"
#include "perturbation_theory_cuda.cu"
