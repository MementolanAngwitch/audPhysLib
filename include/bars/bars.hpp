#pragma once
#include <stdexcept>
#include <numbers>
#include <complex>
#include <vector>
#include <cmath>
#include <algorithm>
#include <array>
#include "tools/brent_method.hpp"

namespace audphys {

struct CrossSection {
	double area; 
	double kappa; //radius of gyration
};

struct Bar {
	double length;
	CrossSection crossSection;
	double youngMod;
	double density;

	Bar(double l_, CrossSection cs_,double y_,double d_) : length(l_), crossSection(cs_),youngMod(y_), density(d_){
		if (!(std::isfinite(l_) && l_ > 0)) throw std::invalid_argument("Invalid bar length");
		if (!(std::isfinite(y_) && y_ > 0)) throw std::invalid_argument("Invalid young's modulus");
		if (!(std::isfinite(d_) && d_ > 0)) throw std::invalid_argument("Invalid density");
		if (!(std::isfinite(cs_.area) && cs_.area > 0))   throw std::invalid_argument("Invalid cross-sectional area");
		if (!(std::isfinite(cs_.kappa) && cs_.kappa > 0)) throw std::invalid_argument("Invalid radius of gyration");
	};
};

struct WaveAmplitude {
	std::complex<double> A;
	std::complex<double> B;
};

struct Load {
	double R = 0.0; // Resistance
	double m = 0.0; // Mass
	double s = 0.0; // Stiffness
};

struct Termination {
	bool fixed = false;
	Load load{};
};

inline Termination fixed_end() {return {true, Load{}}; }
inline Termination free_end() {return {false, Load{}}; }
inline Termination loaded_end(const Load& L) {return {false, L};}

}

