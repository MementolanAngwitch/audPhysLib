/*
--------- Timoshenko Beam Normal Modes ----------------------------



*/

#pragma once
#include "input.hpp"

// Although this is named bars, can also be used to modelling stiff strings. 
namespace audphys::bars{

class Bar{
public:
	Bar(float L_, float S_, float rho_) : L(L_), S(S_), rho(rho_) {};
	float get_rho() const {return rho;};
	float get_length() const {return L;};
	float get_area() const {return S;};
	float get_young_mod() const { return E;};
private:
	float rho; 			   // Volumetric density [kg/m^3]
	float L; 			   // Length [m] 
	float S; 			   // Cross Sectional Area [m^2]
	std::complex<float> E; // Young's modulus [N/m^2]
	float G; 			   // Shear Modulus [GPa]
	float I; 		       // Moment of inertia [m^4]
	float kappa; 		   // shear correction factor
};

class Air {
public:
	Air(float rho_, float eta_) : rho(rho_), eta(eta_) {};
	float get_rho() const { return rho; };
	float get_eta() const { return eta; };
private:
	float rho = 1.2f;    // air density [kg/m^3]
	float eta = 1.5e-5f; // kinematic viscosity [m^2/s]
};
// We will use Timoshenko beam, which considers rotary inertia and shear stress
// Let k' = kappa
// rotary intertia rho i and shear deformation k'GA
// critical frequency  w_c = sqrt(kGA/rhoI)
// critical wavenumber a_c = 1/k sqrt(k'GI + 1), b_c = ~b_c = 0
//                         = 1/k sqrt(1/upsilon^2 + 1)
// We need only consider 4 b.c: fixed-fixed, clamped-clamped, hinged-hinged, and clamped-fixed as normalized wavenumbers are related by:
// a_n^free-sliding = 1/2 a_{2n-1}^free-free
// a_n^free-hinged = 1/2 a_{2n}^free-free
// a_n^clamped-sliding = 1/2 a_{2n-1}^clamped-clamped
// a_n^clamped-hinged = 1/2 a_{2n}^clamped-clamped
// a_n^hinged-sliding = 1/2 a_{2n-1}^hinged-hinged
// a_n^hinged-hinged = 1/2 a_{2n}^hinged-hinged
// a_n^sliding-sliding = 1/2 a_{2n}^sliding-sliding
// a_n^sliding-hinged = 1/2 a_{2n-1}^sliding-sliding

// Timoshenko beam nodes (Han, Benaroya, Wei)
// upsilon^2 = 2(1+v)/(k')
// (upsilon^2 b^2 + a^2)(a^2 upsilon^2 + b^2) / (a^2-b^2)(1+upsilon^2) = s^2
// b = i ~b
// Frequency equations for Timoshenko when a < a_c
// fixed-fixed
// (a^2 - b^2)(a^2+b^2+upsilon^2ab-ab)(a^2+b^2-upsilon^2ab+ab) / (2ab)(b^2+upsilon^2a^2)(a^2+upsilon^2b^2) * sin(a)sinh(b) -cos(a)cosh(b)+1 = 0
// clamped-clamped
// (a^2 + ~b^2)[(a^2 - ~b^2)^2+(a(~b)upsilon^2 - a(~b))^2] / 2a(~b)(-(~b)^2+upsilon^2 a^2)(a^2-upsilon^2(~b)^2) * sin(a)sinh(b)-cos(a)cosh(b) + 1 = 0
// hinged-hinged 
// sin(a)sinh(b) = 0
// clamped-fixed
// (a^2-b^2)sin(a)sinh(b)-ab(a^4+a^4 upsilon^4+4upsilon^4+b^4)/(b^2 +upsilon^2 a^2)(a^2+upsilon^2 b^2) cos(a)cosh(b)

// Frequency equations for Timoshenko when a > a_c
// fixed-fixed
// (a^2 + ~b^2)[(a^2-(~b)^2)^2 + (a(~b)upsilon^2 - a(~b))^2] / (2a(~b))(-(~b)^2+upsilon^2 a^2)(a^2-upsilon^2(~b)^2) * sin(a)sin(~b)-cos(a)cos(~b)+1 = 0
// clamped-clamped
// (a^2+(~b)^2)[(a^2-(~b)^2)^2 + (a(~b)upsilon^2 - a(~b))^2] / (2a(~b))(-(~b)^2+upsilon^2 a^2)(a^2-upsilon^2(~b)^2) * sin(a)sin(~b)-cos(a)cos(~b)+1 = 0
// hinged-hinged
// sin(a)sin(~b) = 0
// clamped-fixed
// (a^2+(~b)^2)sin(a)sin(~b) - a(~b) (a^4+a^4 upsilon^4 - 4upsilon^2 a^2(~b)^2 + (~b)^4 upsilon^4+(~b)^4)/[(-(~b)^2 + (upsilon)^2 a^2)(a^2-upsilon^2(~b)^2)] 
// * cos(a)cos(~b)-2a(~b) = 0
// Solve for a for a < ac by IVP using euler-bernoulli as IV
// Solve for a for a > ac by IVP using a = ac as IV.
class PlayBar{
public:
	PlayBar();
private:
	Bar b;
	Air a; 
}


}