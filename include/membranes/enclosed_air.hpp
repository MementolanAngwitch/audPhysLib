#pragma once

#include <stdexcept>
#include <cmath>
#include <numbers>

#include "membranes/2d_wave.hpp"
#include "membranes/membranes.hpp"
#include "membranes/membrane_modes.hpp"
#include "membranes/circular_fixed_rim.hpp"
#include "tools/bessel.hpp"
#include "tools/brent_method.hpp"

namespace audphys::membranes{

struct AirVolume{
	double eq_pressure;   // P0
	double volume;        // V0
	double heat_capacity; // gamma
	CircularMode mode;
};

inline double volume(const AirVolume& av, double P) {
	return av.volume * std::pow(av.eq_pressure / P, 1.0 / av.heat_capacity);
}

inline double pressure(const AirVolume& av, double V) {
	return av.eq_pressure * std::pow(av.volume / V, av.heat_capacity);
}

inline double symmetric_normal_mode(const CircularMode& mode, double r) {
	const double k = mode.k;
	return mode.A * ( first_kind_bessel(0, k*r) - first_kind_bessel(0, k*mode.a) );
}


inline double excess_pressure(const AirVolume& av){
	const double a = av.mode.a;
	const double averageDisplacement = av.mode.A * first_kind_bessel(2, av.mode.k * a);
	return -av.heat_capacity * (av.eq_pressure / av.volume) * std::numbers::pi * a * a * averageDisplacement;
}


inline double ka_solver(const AirVolume& av){
	const double a = av.mode.a;
	const double B = std::numbers::pi * std::pow(a, 4) * av.heat_capacity * av.eq_pressure
	               / (av.mode.membrane.tension * av.volume);
	auto f = [B](double x){ return first_kind_bessel(0, x) + B * first_kind_bessel(2, x) / (x*x); };
	const double j21 = bessel_zero(2, 1)[0];
	return brent(f, 2.0, j21);
}

inline CircularMode air_loaded_mode(const AirVolume& av){
	CircularMode md = av.mode;
	md.k = ka_solver(av) / md.a;
	md.w = stretched_phase_speed(md.membrane) * md.k;
	return md;
}

// Tension that puts the air-loaded (0,1) mode at angular frequency w [rad/s].
// The air's stiffness depends on tension, so there is no formula: solve w(T) = w numerically.
// Uses the radius, surface density and air from av; the tension stored in av is ignored.

inline double tension_for_ang_freq(const AirVolume& av, double w) {
	const std::size_t n = 1;
	if (!(std::isfinite(w) && w > 0.0)) throw std::invalid_argument("tension_for_ang_freq: angular frequency must be positive and finite");
	if (n < 1) throw std::invalid_argument("tension_for_ang_freq: mode number starts at 1");
	auto err = [&](double T) {
		AirVolume trial = av;
		trial.mode.membrane.tension = T;
		return air_loaded_mode(trial).w - w;
	};
	double lo = 1e-6, hi = 1.0;          // w rises with T, and tends to 0 as T does
	while (err(hi) < 0.0) hi *= 2.0;
	return brent(err, lo, hi);
}

inline void strike_air_loaded(CircularMode& md, double I, double r0) {
	if (r0 < 0.0 || r0 > md.a) throw std::invalid_argument("strike_air_loaded: point must lie on the membrane (0 <= r0 <= a)");
	const double x = md.k * md.a, J0 = first_kind_bessel(0, x), J1 = first_kind_bessel(1, x);
	const double norm = std::numbers::pi * md.a * md.a * (2.0 * J0 * J0 + J1 * J1 - 4.0 * J0 * J1 / x);
	const double psi0 = first_kind_bessel(0, md.k * r0) - J0;
	md.A = I / md.membrane.ps * psi0 / (damped_angular_frequency(md) * norm);
	md.phi = -std::numbers::pi / 2.0;
}



}