//Damped oscilliations
#pragma once
#include "simple_oscillators.hpp"

namespace audphys {
struct DampedOscillator {
	SimpleOscillator base;	
	double Rm;				// kg/s
	double beta;            //B = Rm / 2m

	DampedOscillator(const SimpleOscillator& b, double Rm_) : base(b), Rm(Rm_), beta(Rm_ / (2 * b.m)) {
		if (!(std::isfinite(Rm) && Rm >= 0)) throw std::invalid_argument("DampedOscillator: Rm must be non-negative and finite");
	}
};

//omega_d = (omega0 ^ 2 - beta^2) ^ 1/2
inline double nat_damped_ang_freq(const DampedOscillator& d) {
	double omega0 = nat_angular_freq(d.base);
	assert(d.beta < omega0);
	return std::sqrt(omega0 * omega0 - d.beta * d.beta);
}

inline AmplitudePhase to_amp_phase(const DampedOscillator& d, InitialConditions ic) {
	double A2 = (ic.u0 + d.beta * ic.x0) / nat_damped_ang_freq(d);
	return { std::sqrt(ic.x0 * ic.x0 + A2 * A2), std::atan2(-A2, ic.x0) };
}

struct DampedMotion {
	DampedOscillator osc;
    AmplitudePhase ap;

    // From initial conditions: the conversion uses this oscillator, so they always match
    DampedMotion(DampedOscillator o, InitialConditions ic) : osc(o), ap(to_amp_phase(o.base, ic)) {}

    DampedMotion(DampedOscillator o, AmplitudePhase a) : osc(o), ap(a) {
        if (!(std::isfinite(ap.A) && ap.A >= 0)) throw std::invalid_argument("amplitude must be non-negative");
    }
};

// A = A * exp(-Beta * t)
inline double damped_amplitude(const DampedMotion& mo, double t) {
	return mo.ap.A * std::exp(-1.0 * mo.osc.beta * t);
}

//x = Ae^(-Beta*t) cos (omega_d * t + phi)
inline double position(const DampedMotion& mo, double t) {
	return damped_amplitude(mo, t) * std::cos(nat_damped_ang_freq(mo.osc) * t + mo.ap.phi);
}

inline double relaxation_time(const DampedMotion& mo) {
	return 1.0 / mo.osc.beta; 
}

inline double relaxation_time(const DampedOscillator& d) {
	return 1.0 / d.beta;
}
}