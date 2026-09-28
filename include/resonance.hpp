#pragma once
#include "forced_oscillators.hpp"

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

namespace audphys {

inline double resonant_speed(const DampedOscillator& o, const Drive& d, double t) {
	double w0 = nat_angular_freq(o.base);
	return (d.F / o.Rm) * std::cos(w0 * t);
}

inline double resonant_position(const DampedOscillator& o, const Drive& d, double t) {
	double w0 = nat_angular_freq(o.base);
	return (d.F / (w0 * o.Rm)) * std::sin(w0 * t);
}

inline double quality_factor(const DampedOscillator& o) {
	double w0 = nat_angular_freq(o.base);
	return w0*o.base.m / o.Rm;
}

}