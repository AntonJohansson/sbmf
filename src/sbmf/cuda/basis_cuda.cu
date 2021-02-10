#define OMEGA_cuda 1.0

__host__
static inline f64 ho_K_cuda(const u32 n) {
	const f64 pi_factor = pow(OMEGA_cuda/M_PI,0.25);
	const f64 normalization_factor = pi_factor / sqrt(pow(2.0,n) * (f64)factorial_128_cuda(n));

	return normalization_factor;
}
