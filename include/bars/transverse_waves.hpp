#pragma once 
#include "long_waves.hpp"

namespace audphys{

enum class BendingEnd { Clamped, Free, SimplySupported };

inline std::array<int,2> vanishing_orders(BendingEnd e) {
	switch(e) {
		case(BendingEnd::Clamped): 		   return { 0,1 };
		case(BendingEnd::Free):  		   return { 2,3 };
		case(BendingEnd::SimplySupported): return { 0,2 };
	}
	return { 0,0 };
}
// Bending moment -YSkappa^2 / R
// kappa^2 = 1/S int(r^2 dS)

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

// v = std::sqrt(wc*kappa) where c is longitudinal phase speed
// g = w/v
// y(x,t) = (Ae^gx + Be^-gx)e^jwt + Ce^j(wt+gx) + De^j(wt-gx)
// Re(y(x,t)) = [A cosh gx + B sinh gx + C cos gx + D sin gx] * cos(wt+phi)

inline std::complex<double> complex_transverse_displacement( const Bar& b, std::complex<double> A, std::complex<double> B, 
	std::complex<double> C, std::complex<double> D, double kappa, double w, double t, double x) {
	const double v = std::sqrt(w * longitudinal_phase_speed(b) * kappa);
	const double g = w/v;

	return (A * std::exp(g*x) + B * std::polar(-g*x)) * std::polar(w*t) + C * std::polar(w*t + g*x) + D * std::polar(w*t-g*x);
}

inline std::complex<double> transverse_displacement( const Bar& b, std::complex<double> A, std::complex<double> B, 
	std::complex<double> C, std::complex<double> D, double kappa, double w, double t, double x, double phi) {
	const double v = std::sqrt(w * longitudinal_phase_speed(b) * kappa);
	const double g = w/v;

	return (A * std::cosh(g*x) + B * std::sinh(g*x) + C * std::cos(g*x) + D * std::sin(g*x)) * std::cos(w*t + phi);
}


inline double bending_moment(const Bar& b, double k2, double R) { //TODO: Replace R with [1+(dy/dx)^2]^3/2 / d^2y/dx^2
	return -1.0 * b.youngMod * b.crossSection.area * k2 / R;
}


}//audphys