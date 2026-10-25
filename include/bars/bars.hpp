#pragma once

#include <stdexcept>
#include <cmath>
#include <complex>
#include <vector>

namespace audphys::bars{

struct CrossSection{
	double area;
	double kappa;
};

struct Bar {
	double length;
	CrossSection crossSection;
	double youngMod;
	double density;
	double shearMod;
};

inline double bar_mass(const Bar& b) {
	return b.length * b.crossSection.area * b.density;
}

struct Load {
	double R = 0.0; //Resistance
	double m = 0.0; //Mass
	double s = 0.0; //Stiffness
};

struct Termination {
	bool fixed = false;
	Load load{};
};

struct WaveAmplitude {
	std::complex<double> A;
	std::complex<double> B;
};

struct WaveModeAmplitudes {
	std::vector<WaveAmplitude> amps;
};


}