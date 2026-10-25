#pragma once
#include <cmath>
#include <numbers>
#include <cassert>
#include <stdexcept>
#include <complex>
#include "oscillators.hpp"

namespace audphys::oscillators {

// omega = 2pi * f
inline double nat_angular_freq(const SimpleOscillator& o) {
	return std::sqrt(o.s / o.m);
}

inline double angular_from_frequency(double nat_freq) {
	assert(nat_freq > 0);
	return nat_freq * 2 * std::numbers::pi;
}

inline double nat_frequency(const SimpleOscillator& o) {
	return std::sqrt(o.s/o.m) / (2 * std::numbers::pi);
}

inline double frequency_from_angular(double ang_freq) {
	return ang_freq / (2 * std::numbers::pi);
}

// T = 1/f
inline double period(double f0) {
	assert(f0 > 0);
	return 1.0/f0;
}

//omega_d = (omega0 ^ 2 - beta^2) ^ 1/2
inline double nat_damped_ang_freq(const DampedOscillator& d) {
	double omega0 = nat_angular_freq(d.base);
	assert(d.beta < omega0);
	return std::sqrt(omega0 * omega0 - d.beta * d.beta);
}

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

// Wave equation general solution: x = A1 cos (omega * t) + A2 cos (omega * t)
// A = (x^2 + (u/omega)^2)^1/2, phi = tan^-1 (-u/(x * omega))

inline AmplitudePhase to_amp_phase( const SimpleOscillator& o, InitialConditions ic) {
	double A2 = ic.u0 / nat_angular_freq(o);
	return {
		std::sqrt(ic.x0 * ic.x0 + A2 * A2),
		std::atan2(-A2, ic.x0)
	};
}

inline AmplitudePhase to_amp_phase(const DampedOscillator& d, InitialConditions ic) {
	double A2 = (ic.u0 + d.beta * ic.x0) / nat_damped_ang_freq(d);
	return { std::sqrt(ic.x0 * ic.x0 + A2 * A2), std::atan2(-A2, ic.x0) };
}

struct Motion {
    SimpleOscillator osc;
    AmplitudePhase ap;

    // From initial conditions: the conversion uses this oscillator, so they always match
    Motion(SimpleOscillator o, InitialConditions ic) : osc(o), ap(to_amp_phase(o, ic)) {}

    Motion(SimpleOscillator o, AmplitudePhase a) : osc(o), ap(a) {
        if (!(std::isfinite(ap.A) && ap.A >= 0)) throw std::invalid_argument("amplitude must be non-negative");
    }
};

struct DampedMotion {
	DampedOscillator osc;
    AmplitudePhase ap;

    // From initial conditions: the conversion uses this oscillator, so they always match
    DampedMotion(DampedOscillator o, InitialConditions ic) : osc(o), ap(to_amp_phase(o, ic)) {}

    DampedMotion(DampedOscillator o, AmplitudePhase a) : osc(o), ap(a) {
        if (!(std::isfinite(ap.A) && ap.A >= 0)) throw std::invalid_argument("amplitude must be non-negative");
    }
};

}
