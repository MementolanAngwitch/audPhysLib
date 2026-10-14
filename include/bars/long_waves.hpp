#pragma once

#include <stdexcept>
#include <numbers>
#include <complex>
#include <vector>
#include <cmath>

#include "tools\brent_method.hpp"

namespace audphys{

struct Bar {
	double length;
	double crossSectionArea;
	double youngMod;
	double density;

	Bar(double l_,double a_,double y_,double d_) : length(l_), crossSectionArea(a_), youngMod(y_), density(d_){
		if (!(std::isfinite(l_) && l_ > 0)) throw std::invalid_argument("Invalid bar length");
		if (!(std::isfinite(a_) && a_ > 0)) throw std::invalid_argument("Invalid bar cross sectional area");
		if (!(std::isfinite(y_) && y_ > 0)) throw std::invalid_argument("Invalid young's modulus");
		if (!(std::isfinite(d_) && d_ > 0)) throw std::invalid_argument("Invalid density");
	};
};
// using y to represent longitudinal displacement
inline double strain(const Bar& b,double change) {
	return change/b.length;
}
// hooke's law
inline double stress(const Bar& b, double strain) {
	return -1.0 * b.youngMod * strain;
}


// eq of mtion d^2 y/ d^2 x = 1/c^2 d^2 y/dt^2

inline double phase_speed(const Bar& b) {
	return std::sqrt(b.youngMod / b.density);
}

inline double wave_number(const Bar& b, double w) {
	return w / phase_speed(b);
}


struct WaveAmplitude {
	std::complex<double> A;
	std::complex<double> B;
};

struct WaveModeAmplitudes {
	std::vector<WaveAmplitude> amps;
};
// y(x,t) = A exp(j(wt-kx)) + B exp(j(wt+kx))
inline std::complex<double> position(const WaveAmplitude& amps, const Bar& b, double w, double x, double t) {
	const double k = wave_number(b,w);
	return amps.A * std::polar(1.0,w*t-k*x) + amps.B * std::polar(1.0,w*t+k*x);
}
// u(x,t) = y(t)/dt = jw y(x,t)
inline std::complex<double> velocity(const WaveAmplitude& amps, const Bar& b, double w, double x, double t) {
	return std::complex<double>(0.0, w) * position(amps, b, w,x,t);
}

inline std::complex<double> acceleration(const WaveAmplitude& amps, const Bar& b, double w, double x,double t) {
	return -1.0 * w * w * position(amps,b,w,x,t);
}

// Rigidly fixed at both ends

inline double mode_wave_number(const Bar& b, std::size_t n) {
	return (n * std::numbers::pi / b.length);
}

inline double mode_ang_freq(const Bar& b, std::size_t n) {
	return (n * std::numbers::pi * phase_speed(b)) / b.length;
}

inline std::complex<double> fixed_complex_displacement(const Bar& b, const WaveModeAmplitudes& AmpList, std::size_t n, double x, double t) {
	const std::complex<double> j{ 0.0, 1.0 };
	return -2.0 * j * AmpList.amps[n].A
	     * std::polar(1.0, mode_ang_freq(b, n) * t)		// e^{j w_n t}
	     * std::sin(mode_wave_number(b, n) * x);
}

inline double fixed_real_displacement(const Bar& b, const WaveModeAmplitudes& AmpList, std::size_t n, double x, double t) {
	const double wn = mode_ang_freq(b, n);
	return std::real(AmpList.amps[n].A * std::cos(wn * t) + AmpList.amps[n].B * std::sin(wn*t)) * std::sin(mode_wave_number(b,n)*x);
}

inline double total_fixed_real_displacement(const Bar& b, const WaveModeAmplitudes& AmpList, double x, double t) {
	std::complex<double> out = 0.0;
	for (std::size_t i = 0; i < AmpList.amps.size(); ++i) {
		out += fixed_real_displacement(b,AmpList,i,x,t);
	}
	return std::real(out);
}

//free, free bar

inline std::complex<double> free_complex_displacement(const Bar& b, const WaveModeAmplitudes& AmpList, std::size_t n, double x, double t) {
	const std::complex<double> j{0.0, 1.0};
	return 2.0 * AmpList.amps[n].A 
		   * std::polar(1.0, mode_ang_freq(b,n) * t)
		   * std::cos(mode_wave_number(b,n) * x);
}

inline double free_real_displacement(const Bar& b, const WaveModeAmplitudes& AmpList, std::size_t n, double x, double t) {
	const double wn = mode_ang_freq(b,n);
	return std::real(
		AmpList.amps[n].A * std::cos(wn*t) 
		+ AmpList.amps[n].B * std::sin(wn*t))
		* std::cos(mode_wave_number(b,n) * x);
}

inline double total_free_real_displacement(const Bar& b, const WaveModeAmplitudes& AmpList, double x, double t) {
	std::complex<double> out = 0.0;
	for (std::size_t i = 0; i < AmpList.amps.size(); ++i ) {
		out += free_real_displacement(b, AmpList, i, x,t);
	}
	return std::real(out);
}
// free-fixed bar used free-free displacment eqs.

inline double free_fixed_mode_wave_number(const Bar& b, std::size_t n) {
	return ((2.0*n-1.0) * std::numbers::pi)/(2.0 * b.length);
}

inline double free_fixed_mode_ang_freq(const Bar& b, std::size_t n) {
	return ((2.0*n-1.0)/2.0) * phase_speed(b)/b.length * std::numbers::pi; 
}

//Mass loaded bar with mass at x=L and free at x=0

struct MassLoadedBar {
	Bar bar;
	double mass;
};

inline double bar_mass(const Bar& b) {
	return b.length * b.crossSectionArea * b.density; 
}

// Resonance: tan kL = -(m / m_b) kL, m_b = rho_b * L. No formula exists, so it is solved numerically,
// written as sin(kL) + (m/m_b) kL cos(kL) = 0 to avoid tan's poles. The n-th root lies between
// kL = (n - 1/2) pi (an infinitely heavy mass, a fixed end) and kL = n pi (no mass, a free end).
// f_res = kL * c / (2pi L)

inline double mass_loaded_kL(const Bar& b, double m, int n) {
	const double pi = std::numbers::pi, mu = m/bar_mass(b);
	auto g = [mu](double x) {return std::sin(x) + mu * x * std::cos(x);};
	return Brent(g, (n-0.5) * pi, n *pi);
}

// Nodes of mode n: cos kx = 0 -> x_q = (q-1/2) pi/k q= 1,...,n
inline std::vector<double> mass_loaded_nodes(const Bar& b, double m, int n) {
	const double pi = std::numbers::pi;
	const double k = mass_loaded_kL(b,m,n)/b.length;
	std::vector<double> nodes;
	for (int i = 0; i < n; ++i) {
		nodes.push_back((i - 0.5) * pi/k); 
	}
	return nodes;
}	

//Freely vibrating bar with general B.C

//mech imp  of support at x=0 is zm0 the force acting on this support due to the bar is 
// f0 = -Zm0 u(0,t)

struct GeneralBar {
	std::complex<double> impLeft;
	std::complex<double> impRight;
	Bar bar;  
};

inline std::complex<double> forceLeft(const WaveAmplitude& amps, const GeneralBar& b, double w, double t) {
	return -1.0 * b.impLeft * velocity(amps,b.bar,w,0,t);
}

inline std::complex<double> forceRight(const WaveAmplitude& amps, const GeneralBar& b, double w, double t) {
	return b.impRight * velocity(amps,b.bar,w,b.bar.length,t);
}

// B = (1+(zm0/plc))/(1-(zm0/plc)) A
// tankL = j ((zm0/plc) + (zm/plc))/(1+(zm0/plc)(zm/plc))

inline std::complex<double> mode_shape(const GeneralBar& b, double A) {
	std::complex<double> zm0 = b.impLeft;
	std::complex<double> plc = b.bar.density * phase_speed(b.bar);

	return (1.0+zm0/plc) / (1.0-zm0/plc) * A;
}

//kL = arctan(j ((zm0/plc) + (zm/plc))/(1+(zm0/plc)(zm/plc))) * npi
inline std::complex<double> general_freely_vibrating_kl(const GeneralBar& b, std::size_t n) {
	const std::complex<double> j = {0.0,1.0};
	const std::complex<double> zm0 = b.impLeft;
	const std::complex<double> zm = b.impRight;
	const std::complex<double> plc = b.bar.density * phase_speed(b.bar);

	const std::complex<double> w = j * ((zm0 / plc) + (zm / plc)) / (1.0+(zm0/plc)*(zm/plc));
	return std::atan(w) + n * std::numbers::pi;
}

};