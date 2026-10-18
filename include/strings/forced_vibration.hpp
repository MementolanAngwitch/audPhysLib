#pragma once

#include <numbers>
#include <stdexcept>
#include <complex>
#include <cmath>
#include "tools/brent_method.hpp"
#include "1d_wave.hpp"

namespace audphys::string{

// ---- Forced, fixed string: resonances, nodes, antinodes ----

// Resonance (Z_m0 = 0): f = (2n - 1) c / (4L), n = 1, 2, ...

inline double forced_fixed_res_freq(const String& s, std::size_t n) {
	return ((2.0 * n - 1) * wave_speed(s)) / (4.0 * s.length);
}
 
// Anti-resonance (|Z_m0| infinite): f = n c / (2L), n = 1, 2, ...

inline double forced_fixed_anti_res_freq (const String& s, std::size_t n) {
	return (n*wave_speed(s)) / (2.0 * s.length);
}

// ---- Forced, mass-loaded string: resonances ----

// Resonance: tan kL = -(m / m_s) kL, m_s = rho_L * L. No formula exists, so it is solved numerically,
// written as sin(kL) + (m/m_s) kL cos(kL) = 0 to avoid tan's poles. The n-th root lies between
// kL = (n - 1/2) pi (an infinitely heavy mass, a fixed end) and kL = n pi (no mass, a free end).
// f_res = kL * c / (2pi L)

inline double mass_loaded_resonance_freq(const String& s, double m, double n) {
	const double pi = std::numbers::pi, mu = m / string_mass(s);
	auto g = [mu](double x) { return std::sin(x) + mu * x * std::cos(x); };
	const double kL = brent(g, (n-0.5) * pi, n * pi);
	return kL * wave_speed(s) / (2*pi*s.length);
} 

// ==== Forced vibration of an infinite string ==================================
// Only an outgoing wave: u(x) = (F / rho_L c) e^{-jkx},  input impedance Z_m0 = rho_L c
inline std::complex<double> outgoing_wave(const String& s, double F, double w, double x) {
	const double rho = s.lin_density;
	const double c = wave_speed(s);
	const double k = wave_number(s,w);

	return (F / (rho * c)) * std::polar(1.0, -1.0 * k * x );
}

// Average power input F^2 / (2 rho_L c) 

inline double infinite_string_power(const String& s, double f) {
	return (f * f) /  (2 * characteristic_impedance(s));
}

// right ward infinite string
// k = w/c
// lambda = 2pi/k
// c = lambda * f
// Zm0 = f/u(0,t) for infinite string Zm0 = pl * C, property of string

// Power input
// instantaneous power pi_i  = fu = (Fcos(wt))[(F/pl*C)cos(wt)]
// average power input pi = F^2/(2plc) = 1/2 pl*c*U^2_0, U_0 = |u(0,t)| = F/plc

// ==== Forced vibration of a string of finite length ===========================
// y(x,t) = A e^{j(wt - kx)} + B e^{j(wt + kx)}, with A, B fixed by the two ends:
//   at x = 0:  F = -T dy/dx                       (the driving force)
//   at x = L:  -T dy/dx = Z_L u                   (the load, impedance Z_L)
// Writing z = Z_L / (rho_L c) covers every termination in the chapter:
// Using reflection cofficient R = (Z0-ZL)/(Z0+ZL)

inline WaveAmplitudes forced_amplitudes(const String& s, double F, double w, std::complex<double> R) {
	if (!(w > 0)) throw std::invalid_argument("forced_amplitudes: w must be positive");
	const std::complex<double> j = {0.0,1.0};
	const double k = wave_number(s,w);
	const double L = s.length;
	const double T = s.tension;
	const std::complex<double> e = std::polar(1.0, -2.0 * k * L);
	const std::complex<double> A = F / (j*k*T * (1.0 - R*e));
	return {
		A, R * e * A
	};
} 

inline std::complex<double> input_impedance(const String& s, double w, std::complex<double> R) {
	const double z0 = characteristic_impedance(s);
	const double k = wave_number(s,w);
	const std::complex<double> e = std::polar(1.0, -2.0 * k * s.length);
	return z0 * ((1.0-R*e) / (1.0+R*e));
}

inline double input_power(const String& s, double F, double w, std::complex<double> R) {
	const WaveAmplitudes amps = forced_amplitudes(s,F,w,R);
	return 0.5 * characteristic_impedance(s) * w * w * std::norm(amps.A) * (1.0 - std::norm(R));
}

// Fixed forced finite string
// A = (Fe^jkl) / (2jkTcoskL)
// B = -(Fe^-jkl) / (2jkTcoskL)

// inline WaveAmplitudes loaded_amplitudes(const String& s, double F, double w) {
// 	const double k = wave_number(s,w);
// 	const double L = s.length;
// 	const double T = s.tension;
// 	const std::complex<double> num = F  * std::polar(1.0, k * L);
// 	const std::complex<double> num2 = -F * std::polar(1.0, -k*L);
// 	const std::complex<double> den = std::complex(0.0, 2.0 * k * T * cos(k*L));

// 	return {num / den , num2 / den};
// }

// // Fixed end: Z_m0 = -j rho_L c cot kL 

// inline std::complex<double> input_impedance_forced_fixed_string(const String& s,double w) {
// 	const double rho = s.lin_density; 
// 	const double c = wave_speed(s);
// 	const double k = wave_number(s,w);
// 	const double L = s.length;
// 	return { 0.0, -1.0 * rho * c * (std::cos(k*L) / std::sin(k*L)) };
// }


// // Average power delivered by the driver: (F^2 / 2) Re{Z_m0} / |Z_m0|^2 

// inline double input_power_fixed_forced(const String& s, double F, double w) {
// 	const std::complex<double> zm0 = input_impedance_forced_fixed_string(s,w);//Fix
// 	return (F*F / 2) * std::real(zm0) / std::abs(zm0 * zm0)
// }

// Nodes: x_q = L - q lambda/2, q = 0, 1, 2, ... while x_q >= 0 

inline double forced_fixed_node_pos(const String& s, std::size_t q, double w) {
	return s.length - q * wavelength(s,w) / 2;
}

// Antinodes: x_q = L - (2q - 1) lambda/4, q = 1, 2, ... while x_q >= 0 

inline double forced_fixed_anti_node_pos(const String& s, std::size_t q, double w) {
	return s.length - ((2.0 * q - 1) * wavelength(s,w)) / 4.0; 
}

//Forced, fixed string
// y(x,t) = F/(2jkiTcoskL) * [exp(j[wt+k(L-x)]) - exp(j[wt-k(L-x)])] 
// y(x,t) = F/kT * (sin[k(L-x)])/cos(kL) * e^(jwt)
// Node positions: xq = L - qlambda/2, q = 0,1,2,...,\leq 2L/lambda
// Anti-Node positions: xq = L - q lambda/4, q = 0,1,2,...
// Resonance frequency f_rn  = [(2n-1)/4](c/L)
// Anti-resonance freq f_an = (n/2)(c/L)
// Z_m0 = -j * pl * c * cot(kL)

//Forced, mass-loaded string 
// F = -pl * c^2 (-jkA + jkB)
// A = -(Fe^jkL / 2w*pl*c) ((1+(jwm/plc))/((wm/plc)coskL + sinkL))
// B = -(Fe^jkL / 2w*pl*c) ((1-(jwm/plc))/((wm/plc)coskL + sinkL))
// u(x,t) = -j (F/plc) ((cos[k(L-x)]-(wm/plc)sin[k(L-x)])/((wm/plc)coskL + sinkL)) e^jwt
// Zm0 = jpLc ((wm/plc) + tankL)/(1-(wm/plc)tankL)
// fres: tankL = -(m/ms)kL where ms = pl * L
// node position: tan[k(L-xq)] = plc/wm q = 0,1,2,... \leq 2L/lambda

//Forced, resistance-loaded string
// A = -(Fe^jkl)/(2wplc) * (1+(Rm/plc)/(1+(Rm/jplc)coskL + sinkl))
// B = -(Fe^jkl)/(2wplc) * (1+(Rm/plc)/(1-(Rm/jplc)coskL + sinkl))
// u(x,t) = F/Plc (cos[k(L-x)] + j(Rm/Plc)sin[k(L-x)])/((Rm/plc)coskL + jsinkL) e^jwt
// Zm0 = plc ((Rm/plc) + jtankl)/(1+j(Rm/plc)tankL)
// speed amp. F/plc num/den where num and deno from u(x,t)

}