#pragma once
#include <cmath>
#include <complex>
#include "tools\bessel.hpp" //TODO: modified bessel functions of the first kind 
#include "tools\brent_method.hpp"
namespace audphys::plates{
// Symmetric vibrations of a uniform circular diaphragm
struct CircularThinPlate{
	double p; // volume density
	double sigma; //Poisson's Ratio
	double upsilon; //Young's modulus
	double d; // For this plate the radius of gyration is kappa = d/sqrt(12)
	double a;
};

inline double ga_solver(const CircularThinPlate& p, std::size_t n) {
	if (n ==0) return 0.0;
	auto f = [](double ga){return first_kind_bessel(0,ga) / first_kind_bessel(1,ga) + 
					modified_first_kind_bessel(0,ga)/modified_first_kind_bessel(1,ga); };
	if (n == 1) return 3.20;
	return brent(f, 3.19 + (n-1)*std::numbers::pi, 3.20 + n*std::numbers::pi);
}

inline double fundamental_natural_frequency(const CircularThinPlate& p) {
	return 0.47 * (p.d / (p.a*p.a)) * std::sqrt( p.upsilon / (p.p*(1-p.sigma*p.sigma))); 
}

inline double natural_frequency(const CircularThinPlate& p, std::size_t n) {
	if (n==0) return 0.0;
	if (n==1) return fundamental_natural_frequency(p);
	const double f1 = fundamental_natural_frequency(p);
	const double g1 = ga_solver(p,1) / p.a; 
	const double g = ga_solver(p,n) / p.a;
	return f1 * (g/g1 * g/g1);
}
//TODO: Displacement

}