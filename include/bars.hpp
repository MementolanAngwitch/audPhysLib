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
private:
	Bar b;
	Air a; 
}


}