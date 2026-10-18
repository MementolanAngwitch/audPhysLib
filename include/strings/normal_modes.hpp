#pragma once

#include <numbers>
#include <vector>
#include <cmath>
#include "1d_wave.hpp"

namespace audphys::string{

// ==== Normal modes of the fixed, fixed string ==================================

// k_n = n pi / L  
inline double mode_wave_number(const String& s,std::size_t n) {
	return n * std::numbers::pi / s.length;
}
// w_n = n pi c / L 
inline double mode_ang_freq(const String& s, std::size_t n) {
	return n * std::numbers::pi * wave_speed(s) / s.length;
}

// Build the modes from initial displacement y0(x) and initial velocity v0(x).
// The integrals use the midpoint rule with n_points sub-intervals.
// An = 2/L   int(0,L) y(x,0) sin(knx) dx
// Bn = 2/wnL int (0,L) u(x,0) sin(knx) dx
template <class F, class G>
NormalModes modes_from_initial(const String& s, F y0, G v0, int n_modes, double beta = 0.0, int n_points = 4000){
	const double L = s.length;
	const double dx = L/n_points;
	NormalModes nm{s,{}};
	for (int n = 1; n <= n_modes; ++n) {
		const double k = mode_wave_number(s, n);
		double sy = 0.0;
		double sv = 0.0;
		for (int i = 0; i < n_points; ++i) {
			const double x = (i + 0.5) * dx;
			const double sn = std::sin(k*x);
			sy += y0(x) * sn;
			sv += v0(x) * sn;
		}
		const double An = 2.0 / L * sy * dx;
		const double Vn = 2.0 / L * sv * dx;
		nm.modes.push_back({ n, An, (Vn + beta * An) / mode_ang_freq(s, n), beta });
	}
	return nm;
}


// Plucked: released at rest from y0
// u(x,0) is zero, check that all Bn banish and An = [8h/*npi^2] sin(npi/2) where h is the initial displacement at center
template <class F>
NormalModes pluck_modes(const String& s, F y0, int n_modes, double beta, int n_points = 4000) {
	return modes_from_initial(s, y0, [](double) {return 0.0;}, n_modes, beta, n_points );
}

template <class G>
NormalModes strike_modes(const String& s, G v0, int n_modes, int n_points = 4000) {
	return modes_from_initial(s,[](double){return 0.0;}, v0, n_modes, 0.0, n_points);
}

// y(x,t) = sum_{n=1, ...} (An cos(wnt) + Bn sin(wnt)) * sin(knx)
inline double position(const NormalModes& nm, double x, double t) {
	double y = 0.0;
	for (const Mode& m : nm.modes) {
		const double w = mode_ang_freq(nm.s,m.n);
		const double k = mode_wave_number(nm.s, m.n);
		y += std::exp(-m.beta * t) *
			 (m.A * std::cos(w*t) +
			   m.B * std::sin(w*t) ) 
			 * std::sin(k * x);
	}
	return y; 
}

// u(x,t) = dy/dt  
inline double velocity(const NormalModes& nm, double x, double t) {
	double u = 0.0;
	for (const Mode& m : nm.modes) {
		const double w = mode_ang_freq(nm.s, m.n);
		const double c = std::cos(w * t);
		const double sn = std::sin(w * t);
		u += std::exp(-m.beta * t) * ((-m.beta * m.A + w * m.B) * c + (-m.beta * m.B - w * m.A) * sn)
		   * std::sin(mode_wave_number(nm.s, m.n) * x);
	}
	return u;
}

// dy/dx (x,t): the string's slope. -T dy/dx is the transverse force the string exerts.   
inline double slope(const NormalModes& nm, double x, double t) {
	double d = 0.0;
	for (const Mode& m : nm.modes) {
		const double w = mode_ang_freq(nm.s, m.n);
		const double k = mode_wave_number(nm.s, m.n);
		d += std::exp(-m.beta * t) * (m.A * std::cos(w * t) + m.B * std::sin(w * t)) * k * std::cos(k * x);
	}
	return d;
}

}