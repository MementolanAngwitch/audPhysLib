#pragma once
#include <complex>
#include <cmath>
#include "equation_of_motion.hpp"

namespace audphys::oscillators {

// Complex velocity amplitude u_hat = F / Zm   [m/s]
inline std::complex<double> velocity_phasor(const DampedOscillator& o, const Drive& d) {
	return d.F / mech_impedance(o, d);
}

// Complex displacement amplitude x_hat = F / (jw Zm)   [m]
inline std::complex<double> displacement_phasor(const DampedOscillator& o, const Drive& d) {
	const std::complex<double> jw{ 0.0, d.w };
	return d.F / (jw * mech_impedance(o, d));
}

// Real steady-state displacement x = (F / (w|Zm|)) sin(wt - theta)   [m]
inline double steady_position(const DampedOscillator& o, const Drive& d, double t) {
	return d.F / (d.w * std::abs(mech_impedance(o, d))) * std::sin(d.w * t - phase_angle(o, d));
}

// Real steady-state velocity u = (F / |Zm|) cos(wt - theta)   [m/s]
inline double steady_velocity(const DampedOscillator& o, const Drive& d, double t) {
	return d.F / std::abs(mech_impedance(o, d)) * std::cos(d.w * t - phase_angle(o, d));
}

}
