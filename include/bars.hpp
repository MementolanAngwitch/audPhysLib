#pragma once
#include "input.hpp"

namespace audphys::bars{

class Bar{
public:
	Bar(float L_, float S_, float rho_) : L(L_), S(S_), rho(rho_) {};
	float get_rho() const {return rho;};
	float get_length() const {return L;};
	float get_area() const {return S;};
	float get_young_mod() const { return E;};
private:
	float rho; 			   // Volumetric density [kg/m^2]
	float L; 			   // Length [m] 
	float S; 			   // Cross Sectional Area [m^2]
	std::complex<float> E; // Young's modulus [N/m^2]
};

class Air {
public:
	Air(float rho_, float eta_) : rho(rho_), eta(eta_) {};
	float get_rho() const { return rho; };
	float get_eta() const { return eta; };
private:
	float rho = 1.2;    // air density [kg/m^3]
	float eta = 1.5e-5; // kinematic viscosity [m^2/s]
};

class PlayBar{
public:
	PlayBar()
	// Euler Bernoulli beam (long, thin bar or rod)
	// longitudinal normal modes of vibrations depends upon the end conditions 
	// fixed-fixed or free-free the mode frequenices are gtiven by fn = n cL/2L n = 1,2,3,...
	// fixed-free 												   fm = m cL/4L m = 1,2,3,...
	// for non-perfect support cotkL = m/M kL as in string case

	// for bending waves: (constants should be precomputed in constructor)
	// free-free:	  tan(wL/2v) = +- tanh(wL/2v) => wn = (pi^2 K)/(4L^2) * std::sqrt(E/rho) [3.011^2,5^2,7^2,..,(2n+1)^2]
	// free-clamped:  cot(wL/2v) = +- tanh(wL/2v) => wn = (pi^2 K)/(4L^2) * std::sqrt(E/rho) [1.194^2, 2.988^2, 5^2, ... (2n-1)^2] 
	// simple-simple: wn = (2piK)/(2L^2) sqrt(E/rho) m^2
	// solve for rest by free:		  d^2y/dx^2 = 0, d^3y/dx^3 = 0
	//					 supported:   y=0, d^2y,dx^2 = 0
	// 					 clamped end: y = 0, dy/dx = 0

	// Timoshenko beam, considers rotary inertia and shear stress, preferred for thick bars
	// rotary intertia rho i and shear deformation kappa G A
	// cutof frequency w_c = sqrt(kGA/rhoI)
	// eignsolve bullshit
private:
	Bar b;
	Air a; 
}


}