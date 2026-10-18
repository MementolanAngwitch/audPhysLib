#pragma once

#include <stdexcept>
#include <numbers>
#include <complex>
#include <vector>
#include <cmath>
#include <algorithm>
#include <array>
#include "tools/brent_method.hpp"
#include "bars.hpp"

namespace audphys{
struct CrossSection {
	double area; 
	double kappa; //radius of gyration
};

// rectangular cross section kappa = t/\sqrt(12)
// circular cross section kappa = a/2

inline CrossSection circular_section(double a) {
	return { std::numbers::pi * a * a, a / 2.0 }; 
}

inline CrossSection rectangular_section(double w, double t) {
	return { w*t, t/std::sqrt(12.0)};
}
// general equation 1/S int(r^2 dS) = 1/S int(0,x) int (0,y) y^2 dx dy = 1/S int y^2 w(y)dy 

template <class W> 
inline CrossSection general_cross_section(W width, double y_min, double y_max, int n=1000) {
	double area = simpson([&](double y) { return width(y); }, y_min, y_max, n);
	double ybar = simpson([&](double y) { return y*width(y); }, y_min, y_max, n) / area;  //neutral axis
	double I = simpson([&](double y) { return y*width(y); }, y_min, y_max, n);

	return { area, std::sqrt(I / area) };
}

struct Bar {
	double length;
	CrossSection crossSection;
	double youngMod;
	double density;
	double shearMod;

	Bar(double l_, CrossSection cs_,double y_,double d_) : length(l_), crossSection(cs_),youngMod(y_), density(d_){
		if (!(std::isfinite(l_) && l_ > 0)) throw std::invalid_argument("Invalid bar length");
		if (!(std::isfinite(y_) && y_ > 0)) throw std::invalid_argument("Invalid young's modulus");
		if (!(std::isfinite(d_) && d_ > 0)) throw std::invalid_argument("Invalid density");
		if (!(std::isfinite(cs_.area) && cs_.area > 0))   throw std::invalid_argument("Invalid cross-sectional area");
		if (!(std::isfinite(cs_.kappa) && cs_.kappa > 0)) throw std::invalid_argument("Invalid radius of gyration");
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

inline double longitudinal_phase_speed(const Bar& b) {
	return std::sqrt(b.youngMod / b.density);
}

inline double wave_number(const Bar& b, double w) {
	return w / longitudinal_phase_speed(b);
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
	return (n * std::numbers::pi * longitudinal_phase_speed(b)) / b.length;
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
	return ((2.0*n-1.0)/2.0) * longitudinal_phase_speed(b)/b.length * std::numbers::pi; 
}

//Mass loaded bar with mass at x=L and free at x=0

struct MassLoadedBar {
	Bar bar;
	double mass;
};

inline double bar_mass(const Bar& b) {
	return b.length * b.crossSection.area * b.density; 
}

// Resonance: tan kL = -(m / m_b) kL, m_b = rho_b * L. No formula exists, so it is solved numerically,
// written as sin(kL) + (m/m_b) kL cos(kL) = 0 to avoid tan's poles. The n-th root lies between
// kL = (n - 1/2) pi (an infinitely heavy mass, a fixed end) and kL = n pi (no mass, a free end).
// f_res = kL * c / (2pi L)

inline double mass_loaded_kL(const Bar& b, double m, int n) {
	const double pi = std::numbers::pi, mu = m/bar_mass(b);
	auto g = [mu](double x) {return std::sin(x) + mu * x * std::cos(x);};
	return brent(g, (n-0.5) * pi, n *pi);
}

// Nodes of mode n: cos kx = 0 -> x_q = (q-1/2) pi/k q= 1,...,n
inline std::vector<double> mass_loaded_nodes(const Bar& b, double m, int n) {
	const double pi = std::numbers::pi;
	const double k = mass_loaded_kL(b,m,n)/b.length;
	std::vector<double> nodes;
	for (int i = 0; i < n; ++i) {
		nodes.push_back((i + 0.5) * pi/k); 
	}
	return nodes;
}

struct Load {
	double R = 0.0; // Resistance
	double m = 0.0; // Mass
	double s = 0.0; // Stiffness
};

struct LoadedBar{
	Bar bar;
	Load left;
	Load right;
};

inline std::complex<double> impedance(const Load& L, std::complex<double> w) {
	const std::complex<double> j{ 0.0,1.0 };
	std::complex<double> Z = L.R + j * w * L.m;
	if (L.s >0) Z += L.s / (j*w);
	return Z;
}
//Freely vibrating bar with general B.C

//mech imp  of support at x=0 is zm0 the force acting on this support due to the bar is 
// f0 = -Zm0 u(0,t) using convention of ending at the right, starting at the left

inline std::complex<double> forceLeft(const WaveAmplitude& amps, const LoadedBar& b, double w, double t) {
	return -1.0 * impedance(b.left,w) * velocity(amps,b.bar,w,0.0,t);
}

inline std::complex<double> forceRight(const WaveAmplitude& amps, const LoadedBar& b, double w, double t) {
	return impedance(b.right,w) * velocity(amps,b.bar,w,b.bar.length,t);
}

// B = (1+(zm0/plc))/(1-(zm0/plc)) A
// tankL = j ((zm0/plc) + (zm/plc))/(1+(zm0/plc)(zm/plc))

inline std::complex<double> reflection_ratio(const LoadedBar& b, std::complex<double> w) {
	std::complex<double> zm0 = impedance(b.left,w);
	std::complex<double> plc = b.bar.density * b.bar.crossSection.area * longitudinal_phase_speed(b.bar);

	return (1.0+zm0/plc) / (1.0-zm0/plc);
}

struct Termination {
	bool fixed = false;
	Load load{};
};

inline Termination fixed_end() {return {true, Load{}}; }
inline Termination free_end() {return {false, Load{}}; }
inline Termination loaded_end(const Load& L) {return {false, L};}

struct TerminatedBar {
	Bar bar;
	Termination left;
	Termination right;
};

// Reflection coefficient of an end: R = (1 - z) / (1 + z),  z = Z / (rho c S).
// Finite for every end: free R = +1, fixed R = -1, matched resistance (z = 1) R = 0.
inline std::complex<double> reflection(const Termination& e, double z0, std::complex<double> w) {
	if (e.fixed) return -1.0;
	const std::complex<double> z= impedance(e.load,w) / z0;
	return (1.0-z) / (1.0+z);
}

// tan kL = j (z0 + zL) / (1 + z0 zL),  z = Z / (rho c S), both impedances evaluated at w = c kL / L.
// kL is complex when a load has resistance: Re = frequency, Im = decay. Because the impedances depend
// on kL itself, solve by fixed-point iteration: kL <- arctan(rhs(kL)) + n pi, until it stops changing.
// Mode n (n >= 1) has Re(kL) between (n - 1/2) pi and (n + 1/2) pi.
inline std::complex<double> mode_kl(const TerminatedBar& b, std::size_t n) {
	const std::complex<double> j{ 0.0,1.0 };
	const double pi = std::numbers::pi;
	const double c = longitudinal_phase_speed(b.bar);
	const double L = b.bar.length;
	const double plc = b.bar.density * b.bar.crossSection.area * c;

	auto step = [&](std::complex<double> kl) {
		const std::complex<double> w = c * kl / L;
		const std::complex<double> RR = reflection(b.left, plc, w) * reflection(b.right, plc, w);
		if (std::abs(RR) == 0.0) throw std::domain_error("mode_kl an end absorbs every wave, no modes exist");
		return -0.5 * j * std::log(RR) + n * pi;
	};
	std::complex<double> kl = std::max(n*pi, 0.25*pi);
	for (int it = 0; it < 500; ++it) {
		const std::complex<double> next = step(kl);
		if (std::abs(next - kl) <= 1e-14 * std::abs(next)) return next;
		kl = next;
	} 
	throw std::runtime_error("general_freely_vibrating_kl: iteration did not converge");
}

// The lowest N modes with postive frequency, sorted. 
inline std::vector<std::complex<double>> lowest_modes(const TerminatedBar& tb, int N) {
	std::vector<std::complex<double>> out;
	for (int n = 0; static_cast<int>(out.size()) < N; ++n) {
		std::complex<double> kl;
		try { kl = mode_kl(tb,n); }
		catch (const std::runtime_error&) {
			if (n ==0) continue;
			throw;
		}
		if (kl.real() > 1e-9) out.push_back(kl); // skip kl = 0
	}
	std::sort(out.begin(), out.end(), [](auto a, auto b) {return a.real() < b.real(); });
	out.resize(std::min<std::size_t>(out.size(),N));
	return out;
}


// Frequency and decay rate of a mode, from its complex kL:  w = c kL / L
inline double mode_frequency(const Bar& b, std::complex<double> kL) {
	return longitudinal_phase_speed(b) * kL.real() / (2.0 * std::numbers::pi * b.length);
}

inline double mode_decay_rate(const Bar& b, std::complex<double> kL) {
	return longitudinal_phase_speed(b) * kL.imag() / b.length;
}

// Forced bar, 

struct ForcedBar{
	Bar bar;
	double force;
	Load left;
	Load right;
};

inline std::complex<double> complex_force(const ForcedBar& b, double w, double t) {
	return std::polar(b.force, w* t);
}

inline std::complex<double> forced_input_imp(const ForcedBar& b,double w) {
	double plc = b.bar.crossSection.area * b.bar.density * longitudinal_phase_speed(b.bar);
	double kl = wave_number(b.bar,w) * b.bar.length;
	std::complex<double> zm = impedance(b.right, w);
	std::complex<double> j = { 0.0,1.0 };

	return plc * ((zm / plc) + j * std::tan(kl)) / (1.0 + (zm/plc)*j*std::tan(kl));
}
// Using scaled mech imp Zm/plc = R/plc + jX/plc = r+jx
// xtan^2 (kL) + (r^2 + x^2 -1) tan(kL) - x = 0




}