f64 gaussian(f64 x, f64 mu, f64 sigma) {
	return 1.0/(sigma*sqrt(2.0*M_PI)) * exp(-(x-mu)*(x-mu)/(2*sigma*sigma));
}

void f64_normalize(f64* out, f64* data, u32 size) {
	f64 sum = 0.0;

#pragma omp parallel for shared(data) reduction(+: sum)
	for (u32 i = 0; i < size; ++i) {
		f64 absval = fabs(data[i]);
		sum += absval*absval;
	}

	f64 scaling = 1.0/sqrt(sum);
#pragma omp parallel for
	for (u32 i = 0; i < size; ++i) {
		out[i] = data[i] * scaling;
	}
}
