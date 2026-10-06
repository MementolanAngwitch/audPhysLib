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
double extended_shape(F y0, double L, End left, End right, double xi) {
	const double bc = reflection_sign(right);
	const double flip_per_cycle = reflection_sign(left) * bc;

	const double k = std::floor(xi / (2*L)); // which 2L cycle
	const double v = xi -(k * 2 * L)         // position within the cycle 

	const bool odd_cycle = std::fmod(std::abs(k), 2.0) == 1.0;
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