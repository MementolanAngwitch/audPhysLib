#pragma once 
#include <cmath>
#include "equation_of_motion.hpp"

namespace audphys::oscillator {

// x = x0 cos(omega * t) + (u/omega) sin (omega * t)
inline double position(const SimpleOscillator& o, const InitialConditions& ic, double t) {
	double omega0 = nat_angular_freq(o);
	return ic.x0 * std::cos(omega0 * t) + (ic.u0 / omega0) * std::sin(omega0 * t);
}
// x = A cos (omega * t + phi)
inline double position(const Motion& mo, double t) {
    double w0 = nat_angular_freq(mo.osc);
    return mo.ap.A * std::cos(w0 * t + mo.ap.phi);
}
// U = omega * A
inline double speed_amplitude(const Motion& mo) {
    return nat_angular_freq(mo.osc) * mo.ap.A;
}
// u = -U sin(omega0 * t + phi)
inline double speed(const Motion& mo, double t) {
    double w0 = nat_angular_freq(mo.osc);
    return -speed_amplitude(mo) * std::sin(w0 * t + mo.ap.phi);
}

//a = -omega * U * cos(omega * t + phi)
inline double acceleration(const Motion& mo, double t) {
    return -(mo.osc.s / mo.osc.m) * position(mo, t);
}
//Energy
//Ep = 1/2 sx^2
inline double potential_energy(const SimpleOscillator& o, double x) {
	return 0.5 * o.s * x * x;
}

inline double potential_energy(const Motion& mo, double t) {
	return potential_energy(mo.osc, position(mo, t));
}
//Ek = 1/2 mu^2
inline double kinetic_energy(const SimpleOscillator& o, double u) {
	return 0.5 * o.m * u * u;
}

inline double kinetic_energy(const Motion& mo, double t) {
	return kinetic_energy(mo.osc, speed(mo, t));
}

// E = 1/2 m * omega ** 2 * A **2
inline double total_energy(const Motion& mo) {
	return 0.5 * mo.osc.s * mo.ap.A * mo.ap.A;
}
// E = 1/2 m u0^2 + 1/2 s x0^2  
inline double total_energy(const SimpleOscillator& o, InitialConditions ic) {
	return kinetic_energy(o, ic.u0) + potential_energy(o, ic.x0);
}

}
