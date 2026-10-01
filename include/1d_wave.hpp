#pragma once
#include <numbers>
namespace audphys{

struct String{
	double length;
	double tension;
 	double lin_density}; //assuming uniform density pl
 	String(double l_, t_, l_) : length(l), tension(t_), lin_density(l_){
 		if (!(std::isfinite(length) && length>0)) throw std::invalid_argument("length must be postive and finite");
 		if (!(std::isfinite(tension) && tension>0)) throw std::invalid_argument("tension must be postive and finite");
 		if (!(std::isfinite(lin_density) && lin_density>0)) throw std::invalid_argument("linear density must be postive and finite");
};

// 1D wave equation: d^2y / dx^2 = 1/c^2 d^2y/dt^2 where c^2 is the constant T/pl
// gen solution: y(x,t) = y1(ct-x) + y2(ct+x) where y1 and y2 are functions determined by IC

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
}