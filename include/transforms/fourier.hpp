#pragma once
#include <vector>
#include <complex>
#include <cmath>
#include <numbers>

namespace audphys::transforms {

using cvec = std::vector<std::complex<double>>;

inline std::complex<double> dft_helper(const cvec& x, size_t k) {
	const std::size_t N = x.size();
	std::complex<double> sum{0.0,0.0};
	for(std::size_t n = 0; n < N; ++n) {
		double angle = 2.0 * std::numbers::pi * static_cast<double>(((k*n)%N)/static_cast<double>(N));
		sum += x[n] * std::polar(1.0,-angle);
	}
	return sum;
}	

//X[k] = sum_[n=0,N-1] x[n] * (cos(2pi/N * kn) - jsin(2pi/N * kn))
inline cvec dft(const cvec& x) {
	std::size_t N = x.size();
	cvec output(N);

	for (std::size_t k = 0; k < N; ++k) {
		output[k] = dft_helper(x,k);
	}

	return output;
}

cvec idft(const cvec& X);    // inverse, includes the 1/N factor
cvec fft(const cvec& x);     // radix-2; throws if size is not a power of 2
cvec ifft(const cvec& X);

bool is_power_of_two(std::size_t n);                 // n && !(n & (n - 1))
std::size_t bit_reverse(std::size_t i, int bits);    // for the in-place FFT
double bin_frequency(std::size_t k, std::size_t N, double sample_rate);  // k * fs / N  [Hz]

}