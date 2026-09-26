#pragma once 
#include <cmath>
#include <numbers>
#include <cassert>
#include <stdexcept>
namespace audphys {

struct SimpleOscillator { // should probably be a class but whatever
    double m;   // kg
    double s;   // N/m
    SimpleOscillator(double m_, double s_) : m(m_), s(s_) {
        if (!(std::isfinite(m) && m > 0)) throw std::invalid_argument("mass must be positive and finite");
        if (!(std::isfinite(s) && s > 0)) throw std::invalid_argument("stiffness must be positive and finite");
    }
};
struct InitialConditions {  double x0; double u0; }; 
struct AmplitudePhase { double A; double phi; };

// omega = 2pi * f
inline double angular_freq(SimpleOscillator o) {
	return std::sqrt(o.s / o.m);
}

inline double angular_from_frequency(double nat_freq) {
	assert(nat_freq > 0);
	return nat_freq * 2 * std::numbers::pi;
}

inline double frequency(SimpleOscillator o) {
	return std::sqrt(o.s/o.m) / (2 * std::numbers::pi);
}

inline double frequency_from_angular(double nat_ang_freq) {
	return nat_ang_freq / (2 * std::numbers::pi);
}

// T = 1/f
inline double period(double f0) {
	assert(f0 > 0);
	return 1.0/f0;
}

// Wave equation general solution: x = A1 cos (omega * t) + A2 cos (omega * t)
// A = (x^2 + (u/omega)^2)^1/2, phi = tan^-1 (-u/(x * omega))

inline AmplitudePhase to_amp_phase( SimpleOscillator o, InitialConditions ic) {
	double A2 = ic.u0 / angular_freq(o);
	return {
		std::sqrt(ic.x0 * ic.x0 + A2 * A2),
		std::atan2(-A2, ic.x0)
	};
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

// x = x0 cos(omega * t) + (u/omega) sin (omega * t)
inline double position(SimpleOscillator o, InitialConditions ic, double t) {
	double omega0 = angular_freq(o);
	return ic.x0 * std::cos(omega0 * t) + (ic.u0 / omega0) * std::sin(omega0 * t);
}
// x = A cos (omega * t + phi)
inline double position(const Motion& mo, double t) {
    double w0 = angular_freq(mo.osc);
    return mo.ap.A * std::cos(w0 * t + mo.ap.phi);
}
// U = omega * A
inline double speedAmplitude(const Motion& mo) {
    return angular_freq(mo.osc) * mo.ap.A;
}
// u = -U sin(omega0 * t + phi)
inline double speed(const Motion& mo, double t) {
    double w0 = angular_freq(mo.osc);
    return -speedAmplitude(mo) * std::sin(w0 * t + mo.ap.phi);
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