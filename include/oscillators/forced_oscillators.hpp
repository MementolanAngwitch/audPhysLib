#pragma once
#include "damped_oscillators.hpp"
#include <complex>

namespace audphys {

// Forced oscillations: driving force f = F e^{jwt}   (w = omega, j = imaginary unit)
struct Drive { double F; double w; };	// force amplitude [N], angular frequency [rad/s]

// Reactance Xm = wm - s/w   [kg/s]
inline double mech_reactance(const DampedOscillator& o, const Drive& d) {
	return d.w * o.base.m - o.base.s / d.w;
}

// Impedance Zm = Rm + jXm   [kg/s]
// Magnitude: std::abs(Zm) = (Rm^2 + Xm^2)^(1/2)
inline std::complex<double> mech_impedance(const DampedOscillator& o, const Drive& d) {
	return { o.Rm, mech_reactance(o, d) };
}

// Phase angle theta = atan2(Xm, Rm)   [rad]   (same as std::arg(Zm))
inline double phase_angle(const DampedOscillator& o, const Drive& d) {
	return std::atan2(mech_reactance(o, d), o.Rm);
}

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