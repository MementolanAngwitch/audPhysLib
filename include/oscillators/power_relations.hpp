#pragma once
#include "forced_oscillators.hpp"
#include <complex>
namespace audphys{
// Power Relations
//instantanetous power (watts) PI_i  = (F^2/Zm) cos(wt) * cos(wt - theta)
//average power (watts)        PI = F^2/(2*Zm) * cos(theta) = (F^2Rm)/(2Zm^2) 

inline std::complex<double> instant_power(const DampedOscillator& o, const Drive& d, double t) {
	return ((d.F * d.F) / mech_impedance(o,d)) * std::cos(d.w * t) * cos(d.w *t - phase_angle(o,d));
}
inline std::complex<double> average_power(const DampedOscillator& o, const Drive& d) {
	std::complex<double> Zm = mech_impedance(o,d);
	return std::complex<double>((d.F * d.F * o.Rm)/2.0) / Zm;
} 
}