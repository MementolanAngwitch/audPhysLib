/*
--------- Timoshenko Beam Normal Modes ----------------------------
  rotary intertia rho i and shear deformation k'GA
  upsilon^2 = 2(1+v)/(k')
  critical frequency  w_c = sqrt(k'GA/rhoI)
  critical wavenumber a_c = 1/k sqrt(k'GI + 1), b_c = ~b_c = 0
                          = 1/k sqrt(1/upsilon^2 + 1)
  w-bar = w*L*sqrt(rho/E)

 Dispertion equations:
  a^2 = 1/2(w-bar^2 (1+upsilon^2)+D)
  b^2 = 1/2(-(w-bar)^2 (1+upsilon^2)+D)
  D = sqrt( w-bar^4 (1 - upsilon^2)^2 + 4 w-bar^2 s^2 )


  We need only consider 4 b.c: free-free, clamped-clamped, hinged-hinged, and clamped-free as normalized wavenumbers are related by:
  a_n^free-sliding = 1/2 a_{2n-1}^free-free
  a_n^free-hinged = 1/2 a_{2n}^free-free
  a_n^clamped-sliding = 1/2 a_{2n-1}^clamped-clamped
  a_n^clamped-hinged = 1/2 a_{2n}^clamped-clamped
  a_n^hinged-sliding = 1/2 a_{2n-1}^hinged-hinged
  a_n^hinged-hinged = 1/2 a_{2n}^hinged-hinged
  a_n^sliding-sliding = 1/2 a_{2n}^sliding-sliding
  a_n^sliding-hinged = 1/2 a_{2n-1}^sliding-sliding

  Timoshenko beam modes (Han, Benaroya, Wei)
  upsilon^2 = 2(1+v)/(k')
  (upsilon^2 b^2 + a^2)(a^2 upsilon^2 + b^2) / (a^2-b^2)(1+upsilon^2) = s^2
  For b = i(~b)
  (-upsilon^2 (~b)^2+a^2)(upsilon^2 a^2 - (~b)^2)/ [(a^2+(~b)^2)(1+upsilon^2)] = s^2
  b = i ~b above the cutofff
  Frequency equations for Timoshenko when a < a_c
  free-free
  (a^2 - b^2)(a^2+b^2+upsilon^2ab-ab)(a^2+b^2-upsilon^2ab+ab) / (2ab)(b^2+upsilon^2a^2)(a^2+upsilon^2b^2) * sin(a)sinh(b) -cos(a)cosh(b)+1 = 0
  clamped-clamped
  (a^2 + ~b^2)[(a^2 - ~b^2)^2+(a(~b)upsilon^2 - a(~b))^2] / 2a(~b)(-(~b)^2+upsilon^2 a^2)(a^2-upsilon^2(~b)^2) * sin(a)sinh(b)-cos(a)cosh(b) + 1 = 0
  hinged-hinged 
  sin(a)sinh(b) = 0
  clamped-free
  (a^2+(~b)^2)sin(a)sin(~b) - a(~b)(a^4 + a^4 upsilon^4 - 4 upsilon^2 a^2 (~b)^2 + (~b)^4 upsilon^4 + (~b)^4)
  / [(-(~b)^2 + upsilon^2 a^2)(a^2 - upsilon^2 (~b)^2)] * cos(a)cos(~b) - 2a(~b) = 0

  Frequency equations for Timoshenko when a > a_c
  free-free
  (a^2 + ~b^2)[(a^2-(~b)^2)^2 + (a(~b)upsilon^2 - a(~b))^2] / (2a(~b))(-(~b)^2+upsilon^2 a^2)(a^2-upsilon^2(~b)^2) * sin(a)sin(~b)-cos(a)cos(~b)+1 = 0
  clamped-clamped
  (a^2 - b^2)(upsilon^2 a^2 + upsilon^2 b^2 + upsilon^2 ab - ab)(upsilon^2 a^2 + upsilon^2 b^2 - upsilon^2 ab + ab)
  / (2ab)(b^2 + upsilon^2 a^2)(a^2 + upsilon^2 b^2) * sin(a)sinh(b) - cos(a)cosh(b) + 1 = 0  hinged-hinged
 sin(a)sin(~b) = 0
clamped-free
 (a^2+(~b)^2)sin(a)sin(~b) - a(~b)(a^4 + a^4 upsilon^4 - 4 upsilon^2 a^2 (~b)^2 + (~b)^4 upsilon^4 + (~b)^4)
  / [(-(~b)^2 + upsilon^2 a^2)(a^2 - upsilon^2 (~b)^2)] * cos(a)cos(~b) - 2a(~b) = 0
 Solve for a for a < ac by IVP using euler-bernoulli as IV
 Solve for a for a > ac by IVP using a = ac as IV.

Frequency recovery
 wL sqrt(rho/E) = sqrt((a^2-b^2)/(1+upsilon^2)) below the cutoff, sub (a^2 +(~b)^2) above the cutoff
 Sort by w before passing

 Using frequency scan at fixed s and upsilon
 1. Step the chart frequency (-w) from << to the to top of the band
 2. At each step compute a,b (or ~b) from the dispersion formulas 
 3. Evaluate the frequency equation, refine by bisection
 4. Reject a sign change if |F| doesn't actually shirnk toward zero, indicates a pole
 5. Recover w from w-bar, sort modes by frequency and get each mode's coefficeints from the null vector of [F]

 Phase 1: Modal analysis (design time, double precision)
 For each mode n solve for: frequency w_n, damping ratio, and gain g_n = W_n (x_strike) * W_n(x_pickup) from normalized mode shapes
Phase 2: Modal synthesis (runtime, float)
yn[k] = 2r cos(wnT)yn[k-1]-r^2 yn[k-2] + gn(x)[k]

URTHER CONSIDERATIONS:
Xylophone bars are not perfectly uniform, requires the use of finite element analysis
 Damping?


*/

#pragma once
#include "input.hpp"
#include "tools.hpp"

// Although this is named bars, can also be used to modelling stiff strings. 
namespace audphys::bars{

struct Dispersion{
	double a;
	std::complex<double> b;
	bool above;
}

class Bar{
public:
	Bar(double L_, double A_, double rho_, double v_) : L(L_), A(A_), rho(rho_), v(v_){};
	double get_rho() const {return rho;};
	double get_length() const {return L;};
	double get_area() const {return A;};
	double get_young_mod() const { return E;};
	double get_poisson_ratio() const {return v;};
	double get_kappa() const {return kappa;};
	double get_shear() const {return shear;};
private:
 	double rho; 			   // Volumetric density [kg/m^3]
	double L; 			   // Length [m] 
	double A; 			   // Cross Sectional Area [m^2]
	double E; // Young's modulus [N/m^2]
	double G; 			   // Shear Modulus [GPa]
	double I; 		       // Moment of inertia [m^4]
	double kappa; 		   // Dimensionaless radius of gyration
	double v; 			   // Poisson Ratio
	double shear; 		   // Shear Factor
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

enum class BoundaryConditions {
	Free, Clamped, Hinged;
};

class BarModes{ 
public:
	BarModes(Bar, double);
	Dispersion dispersion(double,double,double);
private:
	Bar bar;
	std::vector<double> a;    // switch to std::array for embedded
	std::vector<double> b;
	std::vector<double> w_bar; // w * L * sqrt(rho/E)
	std::vector<double> w;
	double f_max;
};


}