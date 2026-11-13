#pragma once
#include "input.hpp"

namespace audphys::membranes{
class RectangularMembrane{
public:
	float get_lx() const {return Lx;};
	float get_ly() const {return Ly;};
	float get_tension() const {return tension;};
	float get_sigma() const {return sigma;};

	float mode_frequencies(std::size_t m, std::size_t n); // w_mn = sqrt(T/sigma) * sqrt( (m^2/Lx^2)+(n^2/Ly^2) ), for Lx = Ly, account for degenerate


private:
	float Lx;		// x length [m]
	float Ly;		// y length [m]
	float tension;  // tension [N]
	float sigma;    // area density [kg/m^2]
}

class CircularMembrane{
public:
	float get_radius() const {return a;};
	float get_tension() const {return T;};
	float get_sigma() const {return sigma;};

	float mode_frequencies(std::size_t m); //Solutions are Bessel functions of order m
private:
	float a;		// radius [m]
	float tension;  // tension [N]
	float sigma; 	// area density [kg/m^2]
}
}