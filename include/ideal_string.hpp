#pragma once
#include <numbers>
#include "forced_oscillators.hpp"
namespace audphys{

struct String{
	double length;
	double tension;
 	double lin_density; //assuming uniform density pl
 	String(double l_, t_, ld_) : length(l_), tension(t_), lin_density(ld_){
 		if (!(std::isfinite(length) && length>0)) throw std::invalid_argument("length must be postive and finite");
 		if (!(std::isfinite(tension) && tension>0)) throw std::invalid_argument("tension must be postive and finite");
 		if (!(std::isfinite(lin_density) && lin_density>0)) throw std::invalid_argument("linear density must be postive and finite");
};

//2.3 one-dimensional wave equation d^2y/dx^2 = 1/c^2 d^2y/dt^2 where c^2 = T/Pl
//2.4 general solution              y(x,t) = y1(ct-x) + y2(ct+x)


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