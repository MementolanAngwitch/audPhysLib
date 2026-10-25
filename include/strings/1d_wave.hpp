#pragma once

#include <numbers>
#include <complex>
#include <cmath>
#include "strings.hpp"

namespace audphys::strings{

// c = sqrt(T/pl)
inline double wave_speed(const String& s) {
	return std::sqrt(s.tension / s.lin_density);
}
// rhoL * c
inline double characteristic_impedance(const String& s) {
	return std::sqrt (s.tension * s.lin_density);
}
// k = w/c
inline double wave_number(const String& s, double w) {
	return w / wave_speed(s);
}
// lambda= 2pi/k
inline double wavelength(const String& s, double w) {
	return (2.0 * std::numbers::pi) / wave_number(s,w);
}
// real value of a phasor Re{p e^{jwt}}
inline double real_at(std::complex<double> p, double w, double t) {
	return std::real(p * std::polar(1.0, w*t));
}

//Displacement phasor y(x) = A e^{-jkx} + B e^{jkx} 

inline std::complex<double> displacement_phasor(const String& s, const WaveAmplitudes& amps, double w, double x) {
	const double k = wave_number(s,w);
	return amps.A * std::polar(1.0,-1.0 * k *x) + amps.B * std::polar(1.0,k * x);
}

// Velocity phasor u(x) = j w y(x) 

inline std::complex<double> velocity_phasor(const String& s, const WaveAmplitudes& amps, double w, double x) {
	return std::complex<double>(0.0, w) * displacement_phasor(s, amps, w, x);
}

// y(x,t) = y1(ct-x) + y2(ct+x) with ideal ends fixed end inverts (-1.0) and free end doesn't (1.0)

// B.C Fixed: y(x,t) = y1(ct-x) - y1(ct+x)
// B.C Free:  y(x,t) = y1(ct-x) + y1(ct+x)

inline double reflection_sign(End e) {
	return e == End::Fixed ? -1.0 : 1.0;
}

// Load impedances: mass jwm, resistance rm

inline std::complex<double> mass_impedance(double w, double m) {
	return {0.0, w*m};
}

inline std::complex<double> resistance_impedance(double rm) {
	return {rm, 0.0};
}

inline std::complex<double> reflection_coefficient(const String& s, std::complex<double> ZL) {
	const double Z0 = characteristic_impedance(s);
	return (Z0 - ZL) / (Z0 + ZL);
} 

inline std::complex<double> reflection_coefficient(End e) { return reflection_sign(e); }

}