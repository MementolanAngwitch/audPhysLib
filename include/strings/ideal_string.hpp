#pragma once
#include <numbers>
#include <vector>
#include "oscillators\forced_oscillators.hpp"
namespace audphys{

struct String{
	double length;
	double tension;
 	double lin_density; //assuming uniform density pl
 	String(double l_, double t_, double ld_) : length(l_), tension(t_), lin_density(ld_){
 		if (!(std::isfinite(length) && length>0)) throw std::invalid_argument("length must be postive and finite");
 		if (!(std::isfinite(tension) && tension>0)) throw std::invalid_argument("tension must be postive and finite");
 		if (!(std::isfinite(lin_density) && lin_density>0)) throw std::invalid_argument("linear density must be postive and finite");
 	}
};

struct Mode{
	int n;				//mode number
	double An;
	double Bn;	
	double omega;		
	double decay;		
};

struct StringModes{
	String string;
	std::vector<Mode> modes;
};

// c^2 = T/pl
inline double prop_speed(const String& s) {
	return std::sqrt(s.tension / s.lin_density);
}

// wave eq. y''/t'' = 1/c^2 y''/t''
// gen soln. y(x,t) = y1(ct-x) + y2(ct+x) y1 and y2 are ic

template<class F, class G> // plucked string F and G identical
inline double position(const String& s, F y0, G y1, double x, double t) {
	const double c = prop_speed(s);
	const double L = s.length;
	auto Y0  = [&](double u) {
		double v = std::fmod(u, 2*L);
		if (v<0) v += 2*L; //sign correction
		if (v<=L) return y0(v); 
		return -y0(2*L-v);
	};
	auto Y1  = [&](double u) {
		double v = std::fmod(u, 2*L);
		if (v<0) v += 2*L; //sign correction
		if (v<=L) return y1(v); 
		return -y1(2*L-v);
	};
	return 0.5 * (Y0(x - c * t) + Y1(x + c * t));
}


//Rigidily supported, freely vibrating string
inline double position(double x, double t, const StringModes& sm) {
	// sum_{n=1,inf} (An cos(wnt) + Bn sinwnt) sinknx
	if (x < 0 || x > sm.string.length) throw::std::invalid_argument("position outside of string");
	double y = 0.0;
	for (const Mode& m: sm.modes) {
		y += (m.An * std::cos(m.omega * t) + m.Bn * std::sin(m.omega * t)) * std::sin(m.n * m.omega * t);
	}
	return y;
}



//one-dimensional wave equation d^2y/dx^2 = 1/c^2 d^2y/dt^2 where c^2 = T/Pl wave speed
//general solution              y(x,t) = y1(ct-x) + y2(ct+x)
// 


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