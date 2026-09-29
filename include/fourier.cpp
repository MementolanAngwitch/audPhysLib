#pragma once

#include "fourier.hpp"

//X[k] = sum_[n=0,N-1] x[n] * (cos(2pi/N * kn) - jsin(2pi/N * kn))
cvec dft(const cvec& x) {
	int N = x.size();
	cvec output(N,0.0);

	for (int k = 0; k < N; ++k) {
		output[k] = dft_helper(x,k, N);
	}

	return output;
}

std::complex<double> dft_helper(const cvec& x, int k, int N) {
	std::complex<double> out = {0.0,0.0};
	for(int n = 0; n < N-1; ++n) {
		std::complex<double> entry = {std::cos((2*pi/N) * k*n ),-1 * std::sin((2*pi/N) * k*n )}; 
		entry = entry * x[n];
		out += entry;
	}
	return out;
}	