#pragma once

#include <stdexcept>
#include <numbers>
#include <complex>
#include <vector>
#include <cmath>
#include <algorithm>

#include "tools/brent_method.hpp"
#include "1d_wave.hpp"
#include "bars/bars.hpp"
namespace audphys::bars{


//Mass loaded bar with mass at x=L and free at x=0

struct MassLoadedBar {
	Bar bar;
	double mass;
};

struct LoadedBar{
	Bar bar;
	Load left;
	Load right;
};

struct TerminatedBar {
	Bar bar;
	Termination left;
	Termination right;
};

// Forced bar, 

struct ForcedBar{
	Bar bar;
	double force;
	Load left;
	Load right;
};

//Freely vibrating bar with general B.C

//mech imp  of support at x=0 is zm0 the force acting on this support due to the bar is 
// f0 = -Zm0 u(0,t) using convention of ending at the right, starting at the left

inline std::complex<double> force_left(const WaveAmplitude& amps, const LoadedBar& b, double w, double t) {
	return -1.0 * impedance(b.left,w) * velocity(amps,b.bar,w,0.0,t);
}

inline std::complex<double> force_right(const WaveAmplitude& amps, const LoadedBar& b, double w, double t) {
	return impedance(b.right,w) * velocity(amps,b.bar,w,b.bar.length,t);
}

// B = (1+(zm0/plc))/(1-(zm0/plc)) A
// tankL = j ((zm0/plc) + (zm/plc))/(1+(zm0/plc)(zm/plc))

inline std::complex<double> reflection_ratio(const LoadedBar& b, std::complex<double> w) {
	std::complex<double> zm0 = impedance(b.left,w);
	std::complex<double> plc = b.bar.density * b.bar.crossSection.area * longitudinal_phase_speed(b.bar);

	return (1.0+zm0/plc) / (1.0-zm0/plc);
}


// Fixed-fixed (and free-free) modes

// Rigidly fixed at both ends

inline double mode_wave_number(const Bar& b, std::size_t n) {
	return (n * std::numbers::pi / b.length);
}

inline double mode_ang_freq(const Bar& b, std::size_t n) {
	return (n * std::numbers::pi * longitudinal_phase_speed(b)) / b.length;
}

// Free-fixed modes

// free-fixed bar used free-free displacment eqs.

inline double free_fixed_mode_wave_number(const Bar& b, std::size_t n) {
	return ((2.0*n-1.0) * std::numbers::pi)/(2.0 * b.length);
}

inline double free_fixed_mode_ang_freq(const Bar& b, std::size_t n) {
	return ((2.0*n-1.0)/2.0) * longitudinal_phase_speed(b)/b.length * std::numbers::pi; 
}

//  Mass-loaded modes 

// Resonance: tan kL = -(m / m_b) kL, m_b = rho_b * L. No formula exists, so it is solved numerically,
// written as sin(kL) + (m/m_b) kL cos(kL) = 0 to avoid tan's poles. The n-th root lies between
// kL = (n - 1/2) pi (an infinitely heavy mass, a fixed end) and kL = n pi (no mass, a free end).
// f_res = kL * c / (2pi L)

inline double mass_loaded_kL(const Bar& b, double m, int n) {
	const double pi = std::numbers::pi, mu = m/bar_mass(b);
	auto g = [mu](double x) {return std::sin(x) + mu * x * std::cos(x);};
	return brent(g, (n-0.5) * pi, n *pi);
}

// General terminations (complex kL) 

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
	return longitudinal_phase_speed(b) * kL.real() /   b.length;
}

inline double mode_decay_rate(const Bar& b, std::complex<double> kL) {
	return longitudinal_phase_speed(b) * kL.imag() / b.length;
}


//  Fixed-fixed bar displacement 

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

// Free-free bar displacement 

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

// Mass-loaded bar nodes 

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

// Forced bar 

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