#pragma once
#include <vector>
#include <complex>
#include <cmath>
#include <numbers>

namespace audphys {

using cvec = std::vector<std::complex<double>>;
auto pi = std::numbers::pi;

cvec dft_helper(const cvec& x, int k, int N);
cvec dft(const cvec& x);     // reference: direct from the definition, O(N^2)
cvec idft(const cvec& X);    // inverse, includes the 1/N factor
cvec fft(const cvec& x);     // radix-2; throws if size is not a power of 2
cvec ifft(const cvec& X);

bool is_power_of_two(std::size_t n);                 // n && !(n & (n - 1))
std::size_t bit_reverse(std::size_t i, int bits);    // for the in-place FFT
double bin_frequency(std::size_t k, std::size_t N, double sample_rate);  // k * fs / N  [Hz]

}