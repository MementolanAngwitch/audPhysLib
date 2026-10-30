#pragma once

#include <numbers>
#include <complex>
#include <cmath>
#include <stdexcept>

#include "tools/brent_method.hpp"
#include "tools/composite_simpson.hpp"
#include "bars.hpp"

namespace audphys::bars{

inline double strain(const Bar& b, double change) {
	return change/b.length;
}
//Hooke's Law
inline double stress(const Bar& b, double strain) {
	return -1.0 * b.youngMod * strain;
}

//equation of motion d^2 y / d^2 x = 1 /c^2 d^2 y/dt^2
inline double longitudinal_phase_speed(const Bar& b) {
	return std::sqrt(b.youngMod / b.density);
}

inline double wave_number(const Bar& b, double w) {
	return w / longitudinal_phase_speed(b);
}

//y(x,t) = A exp(j(wt-kx)) + B exp(j(wt+kx))
inline std::complex<double> position(const WaveAmplitude& amps, const Bar& b, double w, double x, double t) {
	const double k = wave_number(b,w);
	return amps.A * std::polar(1.0, w*t - k*x) + amps.B * std::polar(1.0, w*t + k*x);
}

// u(x,t) = y(x,t)/dt = jw y(x,t)
inline std::complex<double> velocity(const WaveAmplitude& amps, const Bar& b, double w, double x, double t) {
	return std::complex<double>(0.0,w) * position(amps,b,w,x,t);
}

//a(x,t) = u(x,t)/dt = -1w^2 y(x,t)
inline std::complex<double> acceleration(const WaveAmplitude& amps, const Bar& b, double w, double x, double t) {
	return -1.0 * w * w * position(amps,b,w,x,t);
}

// rect cross section kappa = t/sqrt(12)
// circular cross section kappa = a/2

inline CrossSection circular_section(double a) {
	return { std::numbers::pi * a * a, a/2.0};
}

inline CrossSection rectangular_section(double w, double t) {
	return { w*t, t/std::sqrt(12.0)};
}

// general equation 1/S int(r^2 dS) = 1/S int(0,x) int (0,y) y^2 dx dy = 1/S int y^2 w(y)dy 

template <class W> 
inline CrossSection general_cross_section(W width, double y_min, double y_max, int n=1000) {
	double area = simpson([&](double y) { return width(y); }, y_min, y_max, n);
	double ybar = simpson([&](double y) { return y*width(y); }, y_min, y_max, n) / area;  //neutral axis
	double I = simpson([&](double y) { return (y-ybar)*(y-ybar)*width(y); }, y_min, y_max, n);

	return { area, std::sqrt(I / area) };
}

inline std::complex<double> impedance(const End& L, std::complex<double> w) {
	if (L.fixed) throw std::domain_error("impedance: a fixed end has infinite impedance");
	const std::complex<double> j{ 0.0,1.0 };
	std::complex<double> Z = L.R + j * w * L.m;
	if (L.s >0) Z += L.s / (j*w);
	return Z;
}

inline End fixed_end() {return {true}; }
inline End free_end() {return {}; }
inline End loaded_end(double R, double m, double s = 0.0) {return {false, R, m, s};}

// Reflection coefficient of an end: R = (1 - z) / (1 + z),  z = Z / (rho c S).
// Finite for every end: free R = +1, fixed R = -1, matched resistance (z = 1) R = 0.
inline std::complex<double> reflection(const End& e, double z0, std::complex<double> w) {
	if (e.fixed) return -1.0;
	const std::complex<double> z= impedance(e,w) / z0;
	return (1.0-z) / (1.0+z);
}

}