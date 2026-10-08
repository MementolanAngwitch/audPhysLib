#pragma once
#include <numbers>
#include <vector>
#include "oscillators\forced_oscillators.hpp"
namespace audphys{

struct String {
	double length;		// L     [m]
	double tension;		// T     [N]
	double lin_density;	// rho_L [kg/m], uniform
 
	String(double L_, double T_, double rhoL_) : length(L_), tension(T_), lin_density(rhoL_) {
		if (!(std::isfinite(length) && length > 0))           throw std::invalid_argument("String: length must be positive and finite");
		if (!(std::isfinite(tension) && tension > 0))         throw std::invalid_argument("String: tension must be positive and finite");
		if (!(std::isfinite(lin_density) && lin_density > 0)) throw std::invalid_argument("String: linear density must be positive and finite");
	}
};
// m = L * pl
inline double string_mass(const String&s) {
	return s.lin_density * s.length;
}
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

// y(x,t) = y1(ct-x) + y2(ct+x) with ideal ends fixed end inverts (-1.0) and free end doesn't (1.0)

enum class End { Fixed, Free };

inline double reflection_sign(End e) {
	return e == End::fixed ? -1.0 : 1.0;
}

// The shape y0 on [0, L], extended to the whole line by mirror images across each end.
// Moving 2L along the line multiplies the shape by (left sign)(right sign).

template<class F>
inline double extended_shape(F y0, double L, End left, End right, double xi) {
	const double bc = reflection_sign(right);
	const double flip_per_cycle = reflection_sign(left) * bc;

	const double k = std::floor(xi / (2*L)); // which 2L cycle
	const double v = xi -(k * 2 * L)         // position within the cycle 

	const bool odd_cycle = std::fmod(std::abs(k), 2.0) == 1.0;
	const double cycle_sign = (flip_per_cycle < 0 && odd_cycle) ? -1.0 : 1.0;

	const double y = (v <= L) ? y0(v) | y0(2*L - v);
	return y * cycle_sign;
}

// String released AT REST from shape y0 (a pluck), ideal ends:
// y(x,t) = 1/2 [ Y0(x - ct) + Y0(x + ct) ] 

template<class F, class G>
inline double position(F y0, F y1, const String& s, End left, End right, double x, double t) {
	const double L = s.length; 
	const double c = wave_speed(s);
	return 1/2 * (extended_shape(y0, L, left, right, x - c*t) + extended_shape(y1,L,left,right,x+c*t));
}

// ==== Forced vibration of an infinite string ==================================
// Only an outgoing wave: u(x) = (F / rho_L c) e^{-jkx},  input impedance Z_m0 = rho_L c
inline double outgoing_wave(const String& s, double f, double w, double x) {
	const double rho = s.lin_density;
	const double c = wave_speed(s);
	const double k = wave_number(s,w);

	return (F / (rho * c)) * std::polar(1, -1.0 * k * x )
}

// Average power input F^2 / (2 rho_L c) 

inline double infinite_string_power(const String& s, double f) {
	return (f * f) /  (2 * characteristic_impedance(s));
}

// ==== Forced vibration of a string of finite length ===========================
// y(x,t) = A e^{j(wt - kx)} + B e^{j(wt + kx)}, with A, B fixed by the two ends:
//   at x = 0:  F = -T dy/dx                       (the driving force)
//   at x = L:  -T dy/dx = Z_L u                   (the load, impedance Z_L)
// Writing z = Z_L / (rho_L c) covers every termination in the chapter:
//   free end z = 0, mass load z = j w m / rho_L c, resistance load z = R_m / rho_L c,
//   fixed end z -> infinity (separate functions below, since z = infinity can't be computed).

struct WaveAmplitudes {
	std::complex<double> A;		// wave travelling away from the driver   [m]
	std::complex<double> B;		// wave reflected back from the load      [m]
};

// Load impedances: mass jwm, resistance rm

inline std::complex<double> mass_impedance(double w, double m) {
	return {0.0, w*m};
}

inline std::complex<double> resistance_impedance(double rm) {
	return {0.0, rm}
}

// Fixed forced finite string
// A = (Fe^jkl) / (2jkTcoskL)
// B = -(Fe^-jkl) / (2jkTcoskL)

inline WaveAmplitudes loaded_amplitudes(const String& s, double F, double w) {
	const double k = wave_number(s);
	const double L = s.length;
	const double T = s.tension;
	const std::complex<double> num = F  * std::polar(1.0, k * L) 
	const std::complex<double> num2 = F * std::polar(-1.0, k*L)
	const std::complex<double> den = std::complex(0, 2.0 * k * T * cos(k*L));

	return {num / den , num2 / den};
}

// Displacement phasor y(x) = A e^{-jkx} + B e^{jkx} 

inline std::complex<double> displacement_phasor_forced_fixed_string(const String& s, const WaveAmplitudes& amps, const double x) {
	const double k = wave_number(s);
	return amps.A * std::polar(1.0,-1.0 * k *x) + amps.B * std::polar(1.0,k * x);
}

// Velocity phasor u(x) = j w y(x) 

inline std::complex<double> velo_phasor_forced_fixed_string(const String& s, double F, double w, double x) {
	std::complex<double> amp = loaded_amplitudes(s, F, w);
	return { 0.0, w * displacement_phasor_forced_fixed_string(s, amps, x) };
}

// Fixed end: Z_m0 = -j rho_L c cot kL 

inline std::complex<double> input_impedance_forced_fixed_string(const String& s) {
	const double rho = s.lin_density; 
	const double c = wave_speed(s);
	const double k = wave_number(k);
	return { 0.0, rho * c * (std::cos(k*L) / std::sin(k*L)) };
}


// Average power delivered by the driver: (F^2 / 2) Re{Z_m0} / |Z_m0|^2 

inline double input_power_fixed_forced(const String& s, double F, double w) {
	const std::complex<double> input_impedance_forced_fixed_string(s,w)
}

// ---- Forced, fixed string: resonances, nodes, antinodes ----

// Resonance (Z_m0 = 0): f = (2n - 1) c / (4L), n = 1, 2, ...

inline double forced_fixed_res_freq(const String& s, std::size_t n) {
	return ((2.0 * n - 1) * wave_speed(s)) / (4.0 * s.length);
}
 
// Anti-resonance (|Z_m0| infinite): f = n c / (2L), n = 1, 2, ...

inline double forced_fixed_anti_res_freq (const String& s, std::size_t n) {
	return (n*wave_speed(s)) / (2.0 * s.length);
}

// Nodes: x_q = L - q lambda/2, q = 0, 1, 2, ... while x_q >= 0 

inline double forced_fixed_node_pos(const String& s, std::size_t q, double w) {
	return s.length - q * wavelength(s,w) / 2;
}

// Antinodes: x_q = L - (2q - 1) lambda/4, q = 1, 2, ... while x_q >= 0 

inline double forced_fixed_anti_node_pos(const String& s, std::size_t q, double w) {
	return s.length - (2.0 * (q - 1) * wavelength(s,w)) / 4.0; 
}


// ---- Forced, mass-loaded string: resonances ----

namespace detail {
	// Bisection: a root of f in [a, b], given f(a) and f(b) have opposite signs (or one is zero)
}

// Resonance: tan kL = -(m / m_s) kL, m_s = rho_L * L. No formula exists, so it is solved numerically,
// written as sin(kL) + (m/m_s) kL cos(kL) = 0 to avoid tan's poles. The n-th root lies between
// kL = (n - 1/2) pi (an infinitely heavy mass, a fixed end) and kL = n pi (no mass, a free end).


// ==== Normal modes of the fixed, fixed string ==================================

struct Mode {
	int n;			// mode number, 1, 2, ...
	double A;		// cosine amplitude   [m]
	double B;		// sine amplitude     [m]
	double beta;	// decay rate         [1/s], 0 = undamped
};
 
struct NormalModes {
	String s;					// the string these modes belong to
	std::vector<Mode> modes;
};

// k_n = n pi / L  
inline double mode_wave_number(const String& s,std::size_t n) {
	return n * std::numbers::pi / s.length;
}
// w_n = n pi c / L 
inline double mode_ang_freq(const String& s, std::size_t n) {
	return n * std::numbers::pi * wave_speed(s) / L;
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
	for (std::size_t n = 0; n < n_modes; ++n) {
		const double k = mode_wave_number(s, k);
		double sy = 0.0;
		double sv = 0.0;
		for (std::size_t i = 0; i < n_points; ++i) {
			const double x = (i + 0.5) * dx;
			const double sn = std::sin(k*x);
			sy += y0(x) * sn;
			sv += v0(x) * sn;
		}
		nm.modes.push_back({
			n,
			2.0 / L * sy * dx,
			2.0 / (mode_ang_freq(s,n) * L) * sv * dx,
			beta
		});
	}
	return nm;
}


// Plucked: released at rest from y0
// u(x,0) is zero, check that all Bn banish and An = [8h/*npi^2] sin(npi/2) where h is the initial displacement at center
template <class F>
NormalModes pluck_modes(const String& s, F y0, int n_modes, double beta, int n_points = 4000) {
	return modes_from_initial(s, y0, [](double) {return 0.0}, n_modes, beta, n_points );
}

template <class G>
NormalModes strike_modes(const String& s, G v0, int n_modes, int n_points = 4000) {
	return modes_from_inital(s,[](double){return 0.0}, v0, n_modes, beta, n_points);
}

// y(x,t) = sum_{n=1, ...} (An cos(wnt) + Bn sin(wnt)) * sin(knx)
inline double position(const NormalModes& nm, double x, double t) {
	double y = 0.0;
	for (const Mode& mode : nm.modes) {
		const double w = mode_ang_freq(nm.s,m.n);
		const double k = wave_number(nm.s, m.n)
		y += ( mode.A * std::cos(w*t) +
			   mode.B * std::sin(w*t) ) 
			 * std::sin(k * x);
	}
	return y; 
}

// u(x,t) = dy/dt  
inline double velocity(const NormalModes& nm, double x, double t) {
	double u = 0.0;
	for (const Mode& m : nm.modes) {
		const double w = mode_angular_frequency(nm.s, m.n);
		const double c = std::cos(w * t)
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
		const double w = mode_angular_frequency(nm.s, m.n);
		const double k = mode_wave_number(nm.s, m.n);
		d += std::exp(-m.beta * t) * (m.A * std::cos(w * t) + m.B * std::sin(w * t)) * k * std::cos(k * x);
	}
	return d;
}

// B.C Fixed: y(x,t) = y1(ct-x) - y1(ct+x)
// B.C Free:  y(x,t) = y1(ct-x) + y1(ct+x)

// right ward infinite string
// k = w/c
// lambda = 2pi/k
// c = lambda * f
// Zm0 = f/u(0,t) for infinite string Zm0 = pl * C, property of string

// Power input
// instantaneous power pi_i  = fu = (Fcos(wt))[(F/pl*C)cos(wt)]
// average power input pi = F^2/(2plc) = 1/2 pl*c*U^2_0, U_0 = |u(0,t)| = F/plc

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