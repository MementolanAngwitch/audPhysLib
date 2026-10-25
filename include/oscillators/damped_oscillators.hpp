//Damped oscilliations
#pragma once
#include <cmath>
#include "equation_of_motion.hpp"

namespace audphys::oscillators {

inline double relaxation_time(const DampedOscillator& d) {
	return 1.0 / d.beta;
}

inline double relaxation_time(const DampedMotion& mo) {
	return 1.0 / mo.osc.beta; 
}

// A = A * exp(-Beta * t)
inline double damped_amplitude(const DampedMotion& mo, double t) {
	return mo.ap.A * std::exp(-1.0 * mo.osc.beta * t);
}

//x = Ae^(-Beta*t) cos (omega_d * t + phi)
inline double position(const DampedMotion& mo, double t) {
	return damped_amplitude(mo, t) * std::cos(nat_damped_ang_freq(mo.osc) * t + mo.ap.phi);
}
}
