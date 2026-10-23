#pragma once

#include <cmath>
#include <vector>
#include <algorithm>
#include "brent_method.hpp"

namespace audphys {
// Using Miller's Backward Recurrence, safe for X < 200
// J_{n-1} (x) = (2n/x) J_n(x) - J_{n+1} (x)
// 1. Choose N such that N >> m and x set J_{N+1} = 0 and J_N = small and recur down to J_0 
//    values are correct up to unknown scale factor
// 2. Find the unknown factor by using identity J_0 + 2(J_2 + J_4 + J_6 + ...) = 1
// 3. Divide the J_m recorded by that sum
inline double first_kind_bessel(std::size_t m, double x) {
	if (x == 0.0) return m == 0 ? 1.0 : 0.0;
	const double ax = std::abs(x);

	std::size_t N = std::static_cast<std::size_t>(x) * m * 40;
	double Jnext = 0.0;
	double JN = 1e-20;
	double Jprev;
	double Jm;

	double sum = 0.0;

	while (N > 0) {
		Jprev = (2.0 * N)/x * JN - Jnext;
		N -= 1;
		if (N == m) Jm = Jprev;
		if (N % 2 == 0 && N != 0) sum += 2.0 * Jprev; 
		if (N == 0) sum += Jprev;
	}

	const double val = Jm / sum;
	if (x < 0 && (m%2)) val = -val; // J_m(-x) = (-1)^m J_m(x)
	return val;



}

// I_n(z) = sum_{k=0, inf} (x/2)^(2k+m) / k!(m+k)!
// a_k = a_k-1 * (x/2)^2 * 1/k * 1/(n+k) 
// Fails when I_m overflows at x ~= 700
inline double modified_first_kind_bessel(std::size_t m, double x) {
	const double constant = (x/2.0) * (x/2.0);
	double term = std::pow(x/2.0,m) / std::tgamma(m+1.0);
	double sum = term;

	for (std::size_t k = 1; k < 500; ++k) {
		term *= constant * (1.0/k) * (1.0/(m+k));
		sum += term;
		if (std::abs(term) <= 1e-17 * std::abs(sum)) break;
	} 
	return sum;
}

inline std::vector<double> bessel_zero(std::size_t m, size_t N) {
	const double h = 1.0;
	double a = m;
	auto f =[m](double x){ return first_kind_bessel(m,x); };
	double fa = f(a);

	std::vector<double> zeros;

	while (zeros.size() < N) {
		const double b = a + h, fb = f(b);
		if (fb == 0.0) zeros.push_back(b);
		if (fa*fb < 0.0) zeros.push_back(brent(f,a,b));
		a = b; fa = fb;
	}
	return zeros;
}

}