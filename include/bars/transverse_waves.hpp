#pragma once 

#include <stdexcept>
#include <numbers>
#include <vector>
#include <array>
#include <complex>
#include <cmath>

#include "1d_wave.hpp"
#include "tools/bracket_finder.hpp"
#include "tools/brent_method.hpp"

namespace audphys::bars{

// Bending moment -YSkappa^2 / R
// kappa^2 = 1/S int(r^2 dS)

inline double bending_moment(const Bar& b, double R) {
    const double k = b.crossSection.kappa;
    return -b.youngMod * b.crossSection.area * k * k / R;
}

// General solution 

// v = std::sqrt(wc*kappa) where c is longitudinal phase speed
// g = w/v,  w = gv, f = 2pi * gv
// y(x,t) = (Ae^gx + Be^-gx)e^jwt + Ce^j(wt+gx) + De^j(wt-gx)
// Re(y(x,t)) = [A cosh gx + B sinh gx + C cos gx + D sin gx] * cos(wt+phi)

inline std::complex<double> complex_transverse_displacement( const Bar& b, std::complex<double> A, std::complex<double> B, 
	std::complex<double> C, std::complex<double> D, double w, double t, double x) {
	const double kappa = b.crossSection.kappa;
	const double v = std::sqrt(w * longitudinal_phase_speed(b) * kappa);
	const double g = w/v;

	const std::complex<double> j(0.0, 1.0);
	return (A * std::exp(g*x) + B * std::exp(-g*x)) * std::exp(j*w*t) 
	       + C * std::exp(j*(w*t + g*x))+ D * std::exp(j*(w*t - g*x));
}

inline double transverse_displacement( const Bar& b, double A, double B, 
	double C, double D, double w, double t, double x, double phi) {
	const double kappa = b.crossSection.kappa;
	const double v = std::sqrt(w * longitudinal_phase_speed(b) * kappa);
	const double g = w/v;

	return (A * std::cosh(g*x) + B * std::sinh(g*x) + C * std::cos(g*x) + D * std::sin(g*x)) * std::cos(w*t + phi);
}


enum class BendingEnd { Clamped, Free, SimplySupported };

inline std::array<int,2> vanishing_orders(BendingEnd e) {
	switch(e) {
	case(BendingEnd::Clamped): 		   return { 0,1 };
	case(BendingEnd::Free):  		   return { 2,3 };
	case(BendingEnd::SimplySupported): return { 0,2 };
	}
	return { 0,0 };
}

struct BentBar {
	BendingEnd left;
	BendingEnd right;
	Bar bar;
};

enum class Ends { ClampedClamped, FreeFree, ClampedFree }; //TODO: SimplySupported

// Natural frequencies 

// g = w/v,  w = gv, f = 1/2pi * gv * n
// Transverse (bending) vibration of a thin bar.
// Dispersion: omega = kappa * c * g^2.
// sign = +1: clamped-clamped or free-free,  cos(gL)cosh(gL) =  1
// sign = -1: clamped-free,                  cos(gL)cosh(gL) = -1

inline std::vector<double> bending_natural_frequency_helper(const Bar& b, int n, int sign) {
	auto f = [sign](double x){ return std::cos(x) - sign/std::cosh(x);};

	const double kappa = b.crossSection.kappa;
	const double L = b.length;
	const double c = longitudinal_phase_speed(b); 
	const double pi = std::numbers::pi;

	std::vector<double> out;
	out.reserve(n);
	for (int i = 1; i <= n; ++i) {
		// Root i lies in [(k-1)π, kπ]; the sign of f flips across it.
		const int k = (sign > 0) ? i+1 : i;
		const double x = brent(f, (k-1)*pi, k*pi);
		const double g = x/L;
		out.push_back(kappa * c * g * g);
	}
	return out;
} 

inline std::vector<double> symmetric_bending_natural_frequency(const Bar& b, int n) { //TODO: Simply Supported Bars
	return bending_natural_frequency_helper(b,n,1);
}

inline std::vector<double> asymmetric_bending_natural_frequency(const Bar& b, int n) {
	return bending_natural_frequency_helper(b,n,-1);
}

// Mode shapes 

// Mode shape of a bending bar.
inline double mode_shape(Ends e, double g, double L, double x) {
    const double a = g * L, u = g * x;
    switch (e) {
        case Ends::ClampedFree: {
            double s = (std::cosh(a) + std::cos(a)) / (std::sinh(a) + std::sin(a));
            return (std::cosh(u) - std::cos(u)) - s * (std::sinh(u) - std::sin(u));
        }
        case Ends::ClampedClamped: {
            double s = (std::cosh(a) - std::cos(a)) / (std::sinh(a) - std::sin(a));
            return (std::cosh(u) - std::cos(u)) - s * (std::sinh(u) - std::sin(u));
        }
        case Ends::FreeFree: {
            double s = (std::cosh(a) - std::cos(a)) / (std::sinh(a) - std::sin(a));
            return (std::cosh(u) + std::cos(u)) - s * (std::sinh(u) + std::sin(u));
        }
    }
    return 0.0;
}


//interior nodes only, free ends are never nodes, clamped ends always are, take freq and find zeros
inline std::vector<double> node_positions(Ends e, double g, double L) { //TODO: FreelySupported
	
	auto y = [&](double x){ return mode_shape(e,g,L,x); };
	const double eps = 1e-6 * L;
	std::vector<double> nodes;

	for (auto[a,b] : bracket_finder(y,eps, L - eps, 1000)){
		nodes.push_back(brent(y,a,b));
	}
	return nodes;
}

}//audphys