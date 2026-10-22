#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <concepts>
#include <numbers>
#include <vector>

namespace audphys::membrane{

template <class M>
concept MembraneMode = requires(const M& md, double x, double z) {
	{ mode_shape(md, x, z) } -> std::convertible_to<double>;
	{ md.w } -> std::convertible_to<double>;   // natural angular frequency, rad/s
	{ md.A }-> std::convertible_to<double>;   // amplitude
	{ md.phi} -> std::convertible_to<double>;   // phase
	{ md.beta } -> std::convertible_to<double>;   // decay rate, 1/s (0 = undamped)
};
// w_d = (w^2 - b^2)^1/2
template <MembraneMode M>
inline double damped_angular_frequency(const M& md) {
	return std::sqrt(md.w * md.w - md.beta * md.beta); 
}

// Displace of one mode at a ponit
template <MembraneMode M>
inline std::complex<double> displacement(const M& md, double x, double z, double t) {
	return md.A * mode_shape(md,x,z) * std::exp(-md.beta*t) * std::polar(1.0, damped_angular_frequency(md)*t+md.phi);
}

template <MembraneMode M>
inline double total_displacement(const std::vector<M>& modes, double x, double z, double t) {
	double out =0.0;
	for (const M& md : modes) out += displacement(md,x,z,t);
	return out;
}
}