#pragma once

#include "membrane.hpp"
#include "2dwave.hpp"
#include <complex>
#include <cmath>

namespace audphys {
struct CircularMembrane{
	Membrane membrane;
	double a; // radius
}

// Bessel's equation d^2 R/dr^2 + 1/r dR/dd + (k^2 - m^2/r^2)R = 0 
// Solutions are bessel functions of order m of the first kind J_m(kr) and second kind Y_m(kr)
// R(r) = AJ_m(kr) + BY_m(kr)
// if the membrane extends across the origin, B =0
// Application of clamped conditions R(a) = 0 and J_m(ka) = 0 yields:
// y_mn (r,theta,t) = A_mn J_m (k_mn r) cos(mtheta + gamma_mn)e^(jw_mn t), k_mn a = j_mn


// natural frequenices
// f_mn = j_mn c/(2pi *a)


// Physical motion is real part of y_mn
// y_mn(r,theta,t) = A_mn J_m (k_mn r) cos(m*theta + gamma_mn) cos(w_mn t + phi_mn) where A_mn exp(j phi_mn)
// and gamma_mn depends on initial location of excitation


// m determines # of radial nodal lines and n determines # of nodal circles. n >= 1

// Bessel function of the first kind: j_n(x) = sum{k=[0,infty)} (-1)^k / (k!(k+n)!) (x/2)^(2k+n)
// a_k = a_k-1 * (x/2)^2 * 1/k * 1/(n+k)

inline double first_kind_bessel(double n, double k) {
	const double constant = (x / 2.0) * (x / 2.0);
	const double a0 = std::pow(x/2.0, n) / std::tgamma(n + 1.0);

	double sum = a0;

	for (std::size_t i = 1; i < 500; ++i) {
		term *= constant * (1/k) * (1/(n+k));
		sum += term;
		if (std::abs(term) <= 1e-17 * std::abs(sum)) break;
	}
	return sum;
}


}