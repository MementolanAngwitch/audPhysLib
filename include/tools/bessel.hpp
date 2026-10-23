#pragma once

#include <cmath>
#include <vector>
#include "brent_method.hpp"

namespace audphys {
// fails at x = 20 
// a_k = a_k-1 * (x/2)^2 * 1/k * 1/(n+k) * -1

inline double first_kind_bessel(std::size_t m, double x) {
	const double constant = -1.0 * (x/2.0) * (x/2.0);
	double term = std::pow(x/2.0,m) / std::tgamma(m+1.0);
	double sum = term;

	for (std::size_t k = 1; k < 500; ++k) {
		term *= constant * (1.0 / k) * (1.0 / (m+k));
		sum += term;
		if (std::abs(term) <= 1e-17 * std::abs(sum)) break;
	}
	return sum;
}
//TODO: Current implementation fails for large x
inline double modified_first_kind_bessel(std::size_t m, double x) {
	return 0.0; //TODO: Modified first kind bessel function
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