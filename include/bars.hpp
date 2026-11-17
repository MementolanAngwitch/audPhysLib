/*
--------- Timoshenko Beam Normal Modes ----------------------------
 Reference: S. M. Han, H. Benaroya, T. Wei, "Dynamics of transversely vibrating beams
 using four engineering theories", J. Sound Vib. 225(5), 935-988 (1999).
 Numbers in [brackets] are that paper's equations and tables.

 Notation (dimensionless unless starred). u^2 below means upsilon^2.
  s        slenderness L* / k*,  s^2 = L*^2 A* / I*
  k        dimensionless radius of gyration = 1/s   (not the shear factor)
  k'       shear factor [Table 3]
  u^2      upsilon^2 = E* / (k'G*)  (= 2(1+v)/k' for isotropic, v = Poisson ratio) [81]
  w-bar    w* L* sqrt(rho* / E*)
  a, b     wave numbers; above the cutoff b = i ~b [75]

 Cutoff [76-79]
  w_c*     = sqrt(k'G*A* / (rho* I*)),   w-bar_c = s/u
  a_c      = (1/k) sqrt(k'GI + 1) = s sqrt(1 + 1/u^2)   (G, I dimensionless here),  b_c = ~b_c = 0

 Dispersion
  D        = sqrt( w-bar^4 (1-u^2)^2 + 4 w-bar^2 s^2 )
  a^2      = 1/2 ( w-bar^2 (1+u^2) + D)
  b^2      = 1/2 (-w-bar^2 (1+u^2) + D)          b^2 < 0 above the cutoff: ~b^2 = -b^2
  inverse  w-bar^2 = 2a^4 / (X + sqrt(X^2 - 4 u^2 a^4)),  X = a^2 (1+u^2) + s^2

 Wave-number relations [110, 111]
  (u^2 b^2 + a^2)(u^2 a^2 + b^2) / [(a^2 - b^2)(1+u^2)] = s^2
  (a^2 - u^2 ~b^2)(u^2 a^2 - ~b^2) / [(a^2 + ~b^2)(1+u^2)] = s^2

 Frequency recovery [115]
  w-bar^2 = (a^2 - b^2)/(1+u^2) below the cutoff,  (a^2 + ~b^2)/(1+u^2) above

 Frequency equations, a < a_c [Table 9]   (b, sinh/cosh)
  free-free
   (a^2-b^2)(a^2+b^2+u^2ab-ab)(a^2+b^2-u^2ab+ab) / [2ab(b^2+u^2a^2)(a^2+u^2b^2)] sin a sinh b - cos a cosh b + 1 = 0
  clamped-clamped
   (a^2-b^2)(u^2a^2+u^2b^2+u^2ab-ab)(u^2a^2+u^2b^2-u^2ab+ab) / [2ab(b^2+u^2a^2)(a^2+u^2b^2)] sin a sinh b - cos a cosh b + 1 = 0
  hinged-hinged
   sin a sinh b = 0
  clamped-free
   (a^2-b^2) sin a sinh b - ab(a^4 + a^4u^4 + 4u^2a^2b^2 + b^4u^4 + b^4) / [(b^2+u^2a^2)(a^2+u^2b^2)] cos a cosh b - 2ab = 0

 Frequency equations, a > a_c [Table 10]  (~b, sin/cos)
  free-free
   (a^2+~b^2)[(a^2-~b^2)^2 + (a~bu^2 - a~b)^2] / [2a~b(u^2a^2-~b^2)(a^2-u^2~b^2)] sin a sin ~b - cos a cos ~b + 1 = 0
  clamped-clamped
   (a^2+~b^2)[(u^2a^2-u^2~b^2)^2 + (u^2a~b - a~b)^2] / [2a~b(u^2a^2-~b^2)(a^2-u^2~b^2)] sin a sin ~b - cos a cos ~b + 1 = 0
  hinged-hinged
   sin a sin ~b = 0
  clamped-free
   (a^2+~b^2) sin a sin ~b - a~b(a^4 + a^4u^4 - 4u^2a^2~b^2 + ~b^4u^4 + ~b^4) / [(u^2a^2-~b^2)(a^2-u^2~b^2)] cos a cos ~b - 2a~b = 0

 Symmetry [Section 3.1, Table 4]: only FF, CC, HH, CF need solving.
  free-sliding  / free-hinged     = symmetric / antisymmetric modes of free-free
  clamped-sliding / clamped-hinged = symmetric / antisymmetric modes of clamped-clamped
  hinged-sliding                   = symmetric modes of hinged-hinged
  sliding-sliding                  = same frequencies as hinged-hinged
  The full beam has slenderness 2s and w-bar 2w (a halves). Rigid-body modes (a = 0) excluded.
  Modes are split by parity at the midpoint, not by index, because frequency curves can cross.

 Solution method (frequency scan at fixed s, u)
  1. Scan in a on a fine grid (uniform in a, so roots are evenly spaced), mapping each a to w-bar
     with the inverse dispersion. Below and above the cutoff are scanned separately.
  2. Each frequency equation is multiplied through by its own denominator (and divided by cosh b
     below the cutoff). This keeps the roots, removes the poles and avoids overflow.
  3. Refine each sign change with Brent's method.
  4. Mode shapes come from the null vector of the 4x4 boundary-condition matrix for the requested
     ends, built with the bounded basis sin, cos, e^{-bx}, e^{-b(1-x)} below the cutoff.
  5. Normalize so that  integral_0^1 (W^2 + Psi^2 / s^2) dx = 1   [129]
  6. Sort by frequency. Hinged-hinged has two modes per a_n (second branch above the cutoff).

 Phase 1: modal analysis (design time, double). Per mode: w_n, damping ratio, gain
          g_n = W_n(x_strike) W_n(x_pickup) from the normalized shapes.
 Phase 2: modal synthesis (runtime, float):
          y_n[k] = 2r cos(w_n T) y_n[k-1] - r^2 y_n[k-2] + g_n x[k],   r = exp(-zeta_n w_n T)

 FURTHER CONSIDERATIONS
  Undercut bars (marimba, xylophone, vibraphone) are not uniform: needs a 1-D Timoshenko FEM, maybe 3d modelling.
  Damping is not part of Timoshenko theory: supply zeta_n per mode or a loss factor (zeta = eta/2).
*/

#pragma once
#include "input.hpp"
#include "tools.hpp"

// Although this is named bars, can also be used to modelling stiff strings. 
namespace audphys::bars{

enum class BoundaryConditions { Free, Clamped, Hinged, Sliding };

struct Dispersion {
	double a;      // oscillating wave number
	std::complex<double> b;      // b below the cutoff, ~b above it
	bool above;    // true when a > a_c
};

class Bar {
public:
	
	Bar(double L_, double A_, double I_, double E_, double G_, double rho_, double kprime_)
		: L(L_), A(A_), I(I_), E(E_), G(G_), rho(rho_), kprime(kprime_) {}

	static Bar isotropic(double L, double A, double I, double E, double v, double rho, double kprime) {
		return Bar(L, A, I, E, E / (2.0 * (1.0 + v)), rho, kprime);
	}

	double get_length() const { return L; }
	double get_area() const { return A; }
	double get_inertia() const { return I; }
	double get_young_mod() const { return E; }
	double get_shear_mod() const { return G; }
	double get_rho() const { return rho; }
	double get_kprime() const { return kprime; }

	double get_s2() const { return L * L * A / I; }              // slenderness squared
	double get_upsilon2() const { return E / (kprime * G); }    // upsilon^2
	double get_wave_speed() const { return std::sqrt(E / rho); } // longitudinal wave speed c [m/s]

private:
	double L;      // Length [m]
	double A;      // Cross-sectional area [m^2]
	double I;      // Second moment of area [m^4]
	double E;      // Young's modulus [Pa]
	double G;      // Shear modulus [Pa]
	double rho;    // Volumetric density [kg/m^3]
	double kprime; // Shear factor [-]
};

// Shear factors k'
inline double kprime_circle(double v)          { return 6.0 * (1.0 + v) / (7.0 + 6.0 * v); }
inline double kprime_rectangle(double v)       { return 10.0 * (1.0 + v) / (12.0 + 11.0 * v); }
inline double kprime_thin_round_tube(double v) { return 2.0 * (1.0 + v) / (4.0 + 3.0 * v); }
inline double kprime_thin_square_tube(double v){ return 20.0 * (1.0 + v) / (48.0 + 39.0 * v); }
inline double kprime_hollow_circle(double v, double m) { // m = r_inner / r_outer
   const double q = (1.0 + m * m) * (1.0 + m * m);
   return 6.0 * (1.0 + v) * q / ((7.0 + 6.0 * v) * q + (20.0 + 12.0 * v) * m * m);
}

class Air {
public:
   Air(float rho_, float eta_) : rho(rho_), eta(eta_) {};
   float get_rho() const { return rho; };
   float get_eta() const { return eta; };
private:
   float rho = 1.2f;    // air density [kg/m^3]
   float eta = 1.5e-5f; // kinematic viscosity [m^2/s]
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