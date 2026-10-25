#pragma once

#include "2dwave.hpp"
#include "membranes.hpp"
#include "membrane_modes.hpp"
#include "circular_fixed_rim.hpp"

namespace audphys::membrane{

struct AirVolume{
	double eq_pressure;
	double volume;
	double heat_capacity;
	Membrane membrane;
	CircularMode mode;
};

inline double volume(const AirVolume& av, double P) {
	return av.volume * std::pow(av.eq_pressure / P, 1.0 / av.heat_capacity);
}

inline double pressure(const AirVolume& av, double V) {
	return av.eq_pressure * std::pow(av.volume,av.heat_capacity) / std::pow(V, av.heat_capacity);
}

inline double excess_pressure(const AirVolume& av){
	const double a = av.mode.a;
	const double k = av.mode.k;
	double averageDisplacement = av.mode.A * (2.0 * first_kind_bessel(1, k*a)) / (k*a);  
	return -av.heat_capacity * (av.eq_pressure / av.volume) * std::numbers::pi * av.mode.a * av.mode.a * averageDisplacement;
}

inline double symmetric_normal_mode(const CircularMode& mode, double r) {
	const double k = mode.k;
	return mode.A * ( first_kind_bessel(0,k*r) - first_kind_bessel(0,k*a) );
}

inline double ka_solver(const AirVolume& av, double lo, double hi){
	double B = std::numbers::pi * std::pow(av.mode.a, 4) * av.heat_capacity * av.eq_pressure / (av.membrane.tension * av.volume);
	auto f = [B](x){ return first_kind_bessel(0,x) + B * first_kind_bessel(2,x) / (x*x); };
	for (std::size_t i =0; i < 100; ++i) {
		double mid = 0.5*(lo + hi);
		if (f(lo) * f(mid) <= 0) hi = mid;
		else lo = mid;
	}
	return 0.5 * lo * hi;
}

}