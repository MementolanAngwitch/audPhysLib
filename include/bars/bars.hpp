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

// One end of a bar: rigidly fixed, or a resistance + mass + stiffness load (all zero = free)
struct End {
	bool fixed = false; // rigid: infinite impedance
	double R = 0.0;     // resistance [kg/s]
	double m = 0.0;     // mass       [kg]
	double s = 0.0;     // stiffness  [N/m]
};

struct WaveAmplitude {
	std::complex<double> A;
	std::complex<double> B;
};

struct WaveModeAmplitudes {
	std::vector<WaveAmplitude> amps;
};


}