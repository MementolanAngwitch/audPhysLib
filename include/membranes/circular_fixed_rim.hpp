#pragma once

// Bessel's equation d^2 R/dr^2 + 1/r dR/dd + (k^2 - m^2/r^2)R = 0 
// Solutions are bessel functions of order m of the first kind J_m(kr) and second kind Y_m(kr)
// R(r) = AJ_m(kr) + BY_m(kr)
// if the membrane extends across the origin, B =0
// Application of clamped conditions R(a) = 0 and J_m(ka) = 0 yields:
// y_mn (r,theta,t) = A_mn J_m (k_mn r) cos(mtheta + gamma_mn)e^(jw_mn t), 
// k_mn a = j_mn where j_mn is the argument that causes j_m to be zero



// natural frequenices
// f_mn = j_mn c/(2pi *a)


// Physical motion is real part of y_mn
// y_mn(r,theta,t) = A_mn J_m (k_mn r) cos(m*theta + gamma_mn) cos(w_mn t + phi_mn) where A_mn exp(j phi_mn)
// and gamma_mn depends on initial location of excitation


// m determines # of radial nodal lines and n determines # of nodal circles. n >= 1

// Bessel function of the first kind: j_n(x) = sum{k=[0,infty)} (-1)^k / (k!(k+n)!) (x/2)^(2k+n)
// a_k = a_k-1 * (x/2)^2 * 1/k * 1/(n+k)


#include "membranes.hpp"
#include "2dwave.hpp"
#include "tools\brent_method.hpp"
#include <complex>
#include <cmath>
#include <vector>

namespace audphys::membrane {
struct CircularMembrane{
	Membrane membrane;
	double a; // radius
};

inline double first_kind_bessel(std::size_t m, double x) {
	const double constant = (x / 2.0) * (x / 2.0);
	const double a0 = std::pow(x/2.0, m) / std::tgamma(m + 1.0);
	double term = a0;
	double sum = a0;

	for (std::size_t i = 1; i < 500; ++i) {
		term *= constant * (1/i) * (1/(m+i));
		sum += term;
		if (std::abs(term) <= 1e-17 * std::abs(sum)) break;
	}
	return sum;
}

inline std::vector<double> bessel_zero(std::size_t m, std::size_t N) {
	auto f = [m](double x){ return first_kind_bessel(m,x); };
	const double h = 1.0; //gap between is always 3 or greater
	std::vector<double> zeros;
	zeros.reserve(N);
	double a = m, fa = f(a); // first zero lies beyond x = m
	while (static_cast<int>(zeros.size()) < N) {
		const double b = a + h, fb = f(b);
		if (fb == 0.0) zeros.push_back(b);
		else if (fa*fb < 0) zeros.push_back(brent(f,a,b)); // zero between a and b
		a = b; fa = fb;
	}
	return zeros;
}

inline circular_wave_number(const CircularMembrane& membrane, std::size_t m, std::size_t n) {
	return bessel_zero(m,n)[n] / membrane.a;
}

inline std::complex<double> complex_circular_displacement(const CircularMembrane& membrane, double phi, double gamma, std::size_t m, size_t n, 
							double r, double theta, double t, double w) {
	const double Amn = 1.0; //TODO: Need amplitude, phi and gamma equation
	const double kmn = circular_wave_number(membrane, m,n);
	return Amn * first_kind_bessel(m,kmn*r) * std::cos(m*theta + gamma) * std::polar(1.0, w*t);
}

inline std::complex<double> real_circular_displacement(const CircularMembrane& membrane, double phi, double gamma, std::size_t m, size_t n, 
							double r, double theta, double t, double w) {
	const double Amn = 1.0;
	const double kmn = circular_wave_number(membrane,m,n);
	return Amn * first_kind_bessel(m,kmn*r)*std::cos(m*theta+gamma)*std::cos(w*t + phi);
}

}