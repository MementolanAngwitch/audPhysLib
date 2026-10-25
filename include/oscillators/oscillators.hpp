#pragma once
#include <cmath>
#include <stdexcept>

namespace audphys::oscillators {

struct SimpleOscillator { // should probably be a class but whatever
    double m;   // kg
    double s;   // N/m
    SimpleOscillator(double m_, double s_) : m(m_), s(s_) {
        if (!(std::isfinite(m) && m > 0)) throw std::invalid_argument("mass must be positive and finite");
        if (!(std::isfinite(s) && s > 0)) throw std::invalid_argument("stiffness must be positive and finite");
    }
};
struct InitialConditions {  double x0; double u0; }; 
struct AmplitudePhase { double A; double phi; };

struct DampedOscillator {
	SimpleOscillator base;	
	double Rm;				// kg/s
	double beta;            //B = Rm / 2m

	DampedOscillator(const SimpleOscillator& b, double Rm_) : base(b), Rm(Rm_), beta(Rm_ / (2 * b.m)) {
		if (!(std::isfinite(Rm) && Rm >= 0)) throw std::invalid_argument("DampedOscillator: Rm must be non-negative and finite");
	}
};

// Forced oscillations: driving force f = F e^{jwt}   (w = omega, j = imaginary unit)
struct Drive { double F; double w; };	// force amplitude [N], angular frequency [rad/s]

}
