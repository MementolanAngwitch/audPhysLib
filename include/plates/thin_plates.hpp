#pragma once
#include <stdexcept>
#include <cmath>
#include <numbers>
#include "tools/bessel.hpp" 
#include "tools/brent_method.hpp"

namespace audphys::plates{
// Symmetric vibrations of a uniform circular diaphragm
struct CircularThinPlate{
	double p; // volume density
	double sigma; //Poisson's Ratio
	double upsilon; //Young's modulus
	double d; // For this plate the radius of gyration is kappa = d/sqrt(12)
	double a;
};


inline double ga_solver(std::size_t n) {
    auto f = [](double x) {
        return first_kind_bessel(0, x) * modified_first_kind_bessel(1, x)
                   / modified_first_kind_bessel(0, x)
             + first_kind_bessel(1, x);
    };
    const double centre = static_cast<double>(n) * std::numbers::pi;
    return brent(f, centre - std::numbers::pi / 2, centre + std::numbers::pi / 2);
}

inline double circular_thin_phase_speed(const CircularThinPlate& p) {
	return std::sqrt( p.upsilon/(p.p * (1 - p.sigma * p.sigma)));
}

// f_n = g_n^2 / (2pi*a^2) * kappa * c
inline double natural_frequency(const CircularThinPlate& p, std::size_t n) {
	if (n < 1) return 0.0;
	const double g = ga_solver(n) / p.a;
	return g*g / (2.0 * std::numbers::pi) * p.d / (std::sqrt(12)) * circular_thin_phase_speed(p);
}

inline double natural_angular_frequency(const CircularThinPlate& p, std::size_t n) {
	return natural_frequency(p,n) * 2.0 * std::numbers::pi;
}

struct PlateMode{
	double g;
	double a;
	double ratio; // B/A = -J0(ga) / I0(ga)
	double w;
	double A = 1.0;
	double phi = 0.0; // startin from rest, phi = -pi/2 for struck plate
};

inline PlateMode make_mode(const CircularThinPlate& p, std::size_t n) {
	const double ga = ga_solver(n);
	return {
		ga / p.a,
		p.a,
		-first_kind_bessel(0,ga) / modified_first_kind_bessel(0,ga),
		natural_angular_frequency(p,n)
	};
}

//y(r,t) = A[J_0 (g * r) + (B/A) * I_0(gr)] * cos(wt + phi)

inline double displacement(const PlateMode& p, std::size_t n, double r, double t) {
	return p.A * (first_kind_bessel(0, p.g *r) + p.ratio * modified_first_kind_bessel(0,p.g * r)) * std::cos(p.w * t + p.phi);
} 

inline double displacement(const CircularThinPlate& p, std::size_t n, double r, double t) {
	return displacement(make_mode(p,n), n, r, t);
}

// Pitch-to-parameter helpers for mode n at angular frequency w [rad/s].
// Plate frequencies scale as d / a^2 (w = kappa c (ga_n / a)^2, kappa = d / sqrt(12)).

inline double radius_for_ang_freq(const CircularThinPlate& p, double w, std::size_t n = 1) {
	if (!(std::isfinite(w) && w > 0.0)) throw std::invalid_argument("radius_for_ang_freq: angular frequency must be positive and finite");
	if (n < 1) throw std::invalid_argument("radius_for_ang_freq: mode number starts at 1");
	return p.a * std::sqrt(natural_angular_frequency(p, n) / w);
}

inline double thickness_for_ang_freq(const CircularThinPlate& p, double w, std::size_t n = 1) {
	if (!(std::isfinite(w) && w > 0.0)) throw std::invalid_argument("thickness_for_ang_freq: angular frequency must be positive and finite");
	if (n < 1) throw std::invalid_argument("thickness_for_ang_freq: mode number starts at 1");
	return p.d * w / natural_angular_frequency(p, n);
}

}