#pragma once

#include "fourier.hpp"

namespace audphys{

namespace {
std::complex<double> dft_helper(const cvec& x, size_t k) {
	const std::size_t N = x.size();
	std::complex<double> sum{0.0,0.0};
	for(std::size_t n = 0; n < N; ++n) {
		double angle = 2.0 * std::numbers::pi * static_cast<double>(((k*n)%N)/static_cast<double>(N));
		sum += x[n] * std::polar(1.0,-angle);
	}
	return sum;
}	

}

//X[k] = sum_[n=0,N-1] x[n] * (cos(2pi/N * kn) - jsin(2pi/N * kn))
cvec dft(const cvec& x) {
	std::size_t N = x.size();
	cvec output(N);

	for (std::size_t k = 0; k < N; ++k) {
		output[k] = dft_helper(x,k);
	}

	return output;
}

}
