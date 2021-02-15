extern "C" {
#include <sbmf/sbmf.h>
#include <assert.h>
#include <stdio.h>

#define sbmf_stack_push(size_in_bytes) \
    sbmf_stack_push_impl(size_in_bytes, __LINE__, __FILE__, __func__)

__host__
void* sbmf_stack_push_impl(u32 size_in_bytes, const u32 linenumber, const char file[], const char func[]);

u32 sbmf_stack_marker();
void sbmf_stack_free_to_marker(u32);

f64 hermite_integral_4(u32 i, u32 j, u32 k, u32 l);

#include "functions_cuda.cu"
//#include "basis_cuda.cu"
#include "indices_cuda.cu"
//#include "hermite_integrals_cuda.cu"
#include "perturbation_theory_cuda.cu"
}
