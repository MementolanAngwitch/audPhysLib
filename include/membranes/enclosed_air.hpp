#pragma once

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

}