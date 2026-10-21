#pragma once 

namespace audphys::membrane {

struct Membrane{
	double ps; // Surface density
	double tension;  // Tension per unit length
	double thickness; 
};

struct MembraneMode {
	int m, n;
	double k;					// k_mn: sets both the shape and omega_mn = c k_mn
	std::complex<double> A;		// amplitude and phase, from initial conditions
	double beta;				// decay rate 
};

}