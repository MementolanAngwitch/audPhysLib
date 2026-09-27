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

// Wave equation general solution: x = A1 cos (omega * t) + A2 cos (omega * t)
// A = (x^2 + (u/omega)^2)^1/2, phi = tan^-1 (-u/(x * omega))

inline AmplitudePhase to_amp_phase( const SimpleOscillator& o, InitialConditions ic) {
	double A2 = ic.u0 / nat_angular_freq(o);
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
inline double position(const SimpleOscillator& o, InitialConditions& ic, double t) {
	double omega0 = nat_angular_freq(o);
	return ic.x0 * std::cos(omega0 * t) + (ic.u0 / omega0) * std::sin(omega0 * t);
}
// x = A cos (omega * t + phi)
inline double position(const Motion& mo, double t) {
    double w0 = nat_angular_freq(mo.osc);
    return mo.ap.A * std::cos(w0 * t + mo.ap.phi);
}
// U = omega * A
inline double speedAmplitude(const Motion& mo) {
    return nat_angular_freq(mo.osc) * mo.ap.A;
}
// u = -U sin(omega0 * t + phi)
inline double speed(const Motion& mo, double t) {
    double w0 = nat_angular_freq(mo.osc);
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

//Lightly Damped oscilliations

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

struct DampedMotion {
	DampedOscillator osc;
    AmplitudePhase ap;

    // From initial conditions: the conversion uses this oscillator, so they always match
    DampedMotion(DampedOscillator o, InitialConditions ic) : osc(o), ap(to_amp_phase(o.base, ic)) {}

    DampedMotion(DampedOscillator o, AmplitudePhase a) : osc(o), ap(a) {
        if (!(std::isfinite(ap.A) && ap.A >= 0)) throw std::invalid_argument("amplitude must be non-negative");
    }
};

inline AmplitudePhase to_amp_phase(const DampedOscillator& d, InitialConditions ic) {
	double A2 = (ic.u0 + d.beta * ic.x0) / nat_damped_ang_freq(d);
	return { std::sqrt(ic.x0 * ic.x0 + A2 * A2), std::atan2(-A2, ic.x0) };
}

// A = A * exp(-Beta * t)
inline double dampedAmplitude(const DampedMotion& mo, double t) {
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

//Forced oscillations, (w = omega, j is the complex constant)
//complex displacement x = (1/jw) * Fe^{jwt}/(Rm + j(wm-s/w))

//complex speed u = (Fe^{jwt})/(Rm + j(wm - s/w))

//impedance Zm = Rm + jXm

//reactance Xm = wm - s/w

//impedance magnitude Zm = (Rm^2 + (wm -s/w)^2)^(1/2)

//phase angle theta = tan^-1 (Xm/Rm) = tan-1 ((wm-s/w)/Rm) use atan2

//Zm = Zm e^(j theta)

//Simplfiied impedance Zm = f/u, f = Fexp(jwt)

//Simplified equations of motion
// u = f/Zm
// x = f/(j*w*Zm)

//Actual (real) value
// x = (F/(w*Zm) sin(wt - theta))
// u = (F/Zm)cos(wt - theta)

// Power Relations
//instantanetous power (watts) PI_i  = (F^2/Zm) cos(wt) * cos(wt - theta)
//average power (watts)        PI = F^2/(2*Zm) * cos(theta) = (F^2Rm)/(2Zm^2) 

//Mechanical Resonance 
//w0 is where Xm vanishes and Zm = Rm
//Resonance values
// u_res = (F/Rm) cos(w0 t)
// x_res = (F/(w0*Rm)) sin(w0 * t)
// Quality factor Q = w0/(w_u - w_l)
// w_u and w_l are the two angular frequenices above and below resoance at which avg is half its resonance value
// Q = w0m /Rm
// Q = w0/(2beta)
// Q = 1/2 w0 t (t = relaxation time)

//Linear Combination
//Same angular frequencies
// x = A cos(wt+phi)
// A = [(SUM(An*cos(phi_n))^2)^2 + (SUM(A_n*sin(phi_n)))^2]^(1/2)
// tan(phi) = SUM(A_n * sin(phi_n)) / SUM(A_n * cos(phi_n))
// different angular frequenices
// w2 = w1 + delta(w)
// x = Ae^(j(w1*t + phi))
// A = [A1^2 +A2^2 + 2A1A2 * cos(phi1 - phi2 - delta*w*t)]^1/2
// tan(phi) = (A1 * sin(phi_1) + A2 * sin(phi_2 + delta*w*t))/(A1 * cos(phi_1) + A2 cos(phi_2 + delta*w*t))
} 