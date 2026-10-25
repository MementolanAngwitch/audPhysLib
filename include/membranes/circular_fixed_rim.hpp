#pragma once 

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <vector>
#include "2dwave.hpp"
#include "membrane_modes.hpp"
#include "tools/bessel.hpp"

namespace audphys::membrane{

struct CircularMode {
	Membrane membrane;
	std::size_t m,n; // nodal diamters, nodal circles
	double a; // radius
	double k; // j_mn/a
	double w; // natural angular frequency
	double gamma = 0.0; // orientation of nodal features
	double A = 1.0, phi = 0.0, beta = 0.0;
};

// w = c·jmn/a

inline CircularMode circular_mode(const Membrane& mem, double a, std::size_t m, std::size_t n) { // don't use, inefficeint 
	const double k = bessel_zero(m,n)[n-1] / a;
	const double c = stretched_phase_speed(mem);
	return{mem,m,n,a,k,c*k};
}

inline std::vector<CircularMode> circular_modes(const Membrane& mem, double a, std::size_t m_max, std::size_t n_max) {
	std::vector<CircularMode> out;
	out.reserve(n_max * (m_max+1));
	const double c = stretched_phase_speed(mem);
	for (std::size_t m = 0; m <= m_max; ++m) {
		const std::vector<double> zeros = bessel_zero(m,n_max);
		for (std::size_t n = 1; n <= n_max; ++n) {
			const double k = bessel_zero(m,n)[n-1] / a;
			out.push_back({mem,m, n, a, k, c*k});
		}
	}
	std::sort(out.begin(), out.end(), [](const CircularMode& p, const CircularMode& q) { return p.w < q.w; });
	return out;
}

// Polar form. Zero outside the membrane.
inline double mode_shape_polar(const CircularMode& md, double r, double theta) {
	if (r < 0 || r > md.a) return 0.0;
	return first_kind_bessel(md.m, md.k * r) * std::cos(md.m * theta + md.gamma);
}

// Cartesian form origin at the centre.
inline double mode_shape(const CircularMode& md, double x, double z) {
	return mode_shape_polar(md, std::hypot(x, z), std::atan2(z, x));
}

// Interior nodal circles r = j_mq / k, q = 1..n-1 (the rim, q = n, is left out).
inline std::vector<double> nodal_circles(const CircularMode& md) {
	std::vector<double> out;
	if (md.n < 2) return out;
	const std::vector<double> zeros = bessel_zero(md.m, md.n - 1);
	for (double j : zeros) out.push_back(j / md.k);
	return out;
}

// Nodal diameters: angles where cos(m theta + gamma) = 0, theta = ((q + 1/2) pi - gamma) / m, q = 0..m-1.
inline std::vector<double> nodal_diameters(const CircularMode& md) {
	const double pi = std::numbers::pi;
	std::vector<double> out;
	for (std::size_t q = 0; q < md.m; ++q) {
		double theta = std::fmod(((q + 0.5) * pi - md.gamma) / md.m, pi);
		if (theta < 0) theta += pi;
		out.push_back(theta);
	}
	std::sort(out.begin(), out.end());
	return out;
}

inline double average_displacement_amplitude(const CircularMode& md) {
	if (md.m != 0) return 0.0;
	const double ka = md.k * md.a;
	return md.A * 2.0 * first_kind_bessel(1, ka) / ka;
}

inline void strike(std::vector<CircularMode>& mds, double I, double r0, double theta0) {
	// check that r0,theta0 is on the membrane
	for (CircularMode& md : mds) {
		const double V = I / md.membrane.ps;
		double epsilon = 2.0;
		if (md.m==0) epsilon = 1.0;

		md.A = V / (std::numbers::pi * md.a * md.a) * (epsilon / md.w) * 
		(first_kind_bessel(md.m, md.k*r0) / (first_kind_bessel(md.m+1,md.k*md.a)*first_kind_bessel(md.m+1,md.k*md.a)));
		md.gamma = -static_cast<double>(md.m) * theta0;
	}
}

inline double mode_displacement(const CircularMode& md, double r, double theta, double t) {
	return md.A * std::exp(-md.beta * t) * mode_shape_polar(md,r,theta) * std::sin(md.w * t + md.phi);
}

inline double total_mode_displacement(const std::vector<CircularMode>& mds, double r, double theta, double t) {
	double sum = 0.0;
	for (const CircularMode& md : mds) {
		sum += mode_displacement(md,r,theta,t);
	}
	return sum;
}

}