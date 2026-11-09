#pragma once
#include <cmath>
#include <cstddef>
#include <numbers>
#include <stdexcept>
#include <vector>
#include <concepts>
#include <complex>
#include "tools.hpp"
#include "input.hpp"

// TODO: Route class Input into strings, final syntax should like something like 
// gen(String, Air, Input) -> std::vector<float> with size t/dt. 
namespace audphys::strings {

enum class BoundaryCondition { Fixed, Free }; // TODO: Replace BoundaryCondition with Support since real supports are not perfectly fixed

class Support {
	double m; 
};

class String {
public:
	String(double L_, double mu_, double T_, double r_, std::complex<double> E_,
	       BoundaryCondition left_ = BoundaryCondition::Fixed,
	       BoundaryCondition right_ = BoundaryCondition::Fixed)
		: L(L_), mu(mu_), T(T_), r(r_), E(E_), left(left_), right(right_) {}
	double get_length() const { return L; }
	double get_mu() const { return mu; }
	double get_rho() const { return mu / std::numbers::pi * r * r; }
	double get_tension() const { return T; }
	double wavespeed() const { return std::sqrt(T / mu); }
	double get_radius() const { return r;}
	std::complex<double> get_young_mod() const { return E;}
	std::complex<double> impedance() const { return std::sqrt(T * mu); }
	BoundaryCondition get_left() const { return left; }
	BoundaryCondition get_right() const { return right; }
	// Folds a position back onto the string and returns it.
	double check_length(double x) const {
		if (std::abs(x) > L) x = L - (x - L);
		if (x < 0) x = -x;
		return x;
	}
private:
	double L;        // length [m]
	double mu;       // linear density [kg/m]
	double T;        // tension [N]
	double r;        // radius [m]
	std::complex<double> E;  		 // Young's modulus [N/m^2]
	std::size_t parts = 0; // # of divisions of the string
	BoundaryCondition left;  // negative direction
	BoundaryCondition right; // positive direction
};

class Air {
public:
	Air(double rho_, double eta_) : rho(rho_), eta(eta_) {};
	double get_rho() const { return rho; };
	double get_eta() const { return eta; };
private:
	double rho = 1.2;    // air density [kg/m^3]
	double eta = 1.5e-5; // kinematic viscosity [m^2/s]
};


// Damping decay times: A = A_0 * e^(-t/tau)
double tau_air(Air air, String s, double w, double v);

double tau_internal(String s, double w, double v);

double tau_support(String s, double w);

double tau_total(Air air, String s, double w, double v);


// q_n'' + (2/tau_n) q_n' + w_n^2 q_n = (2 / (mu L)) * phi_n(x0) * F(t),
// with mode shape phi_n(x) = sin(n pi x / L).
// On a fixed-fixed string
class PlayString {
public: 
	PlayString(String s_, Air air_, std::size_t n) : s(s_), air(air_), q(n, 0.0), qdot(n, 0.0), w(n), tau(n) {
		const double L = s.get_length();
		for (std::size_t i = 0; i < n; ++i) {
			w[i] = (i+1) * std::numbers::pi * s.wavespeed() / L; // f_n = n v/2L -> wn = f_n * 2pi = (n*v*pi)/L
			tau[i] = tau_total(air_, s_, w[i], 0.0);
		}
	}
	double shape(std::size_t i, double x) const;
	void pluck(double x0, double h);
	void strike(double x0, double J);
	void step(double dt, double x0, double F);
	double displacement(double x) const;

private:
	String s;
	Air air; 
	std::vector<double> tau; 			  // decay time [s]
	std::vector<double> q;    // modal displacement [m]
	std::vector<double> qdot; // modal velocity [m/s]
	std::vector<double> w;	  // angular frequency [rad/s]
};
// A hammer with a compressible felt tip. Force exists only while the felt is squeezed:
// F = K * compression^p  if compression > 0, else 0.
class Hammer {
public:
	Hammer(double M_, double K_, double p_, double x0_, double v0_) 
	 : M(M_), K(K_), x0(x0_), v(v0_) {};
	void step(PlayString& s, double dt);
private:
	double M; 	    // hammer mass [kg]
	double K;       // felt stiffness [n/m^p]
	double p; 	    // felt stiffness exponent 
	double x0; 	    // strike point [m]
	double v;       // hammer velocity [m/s], postive toward string
	double y = 0.0; // hammer position [m], starts touching the string at rest
};

template<class F, class G>
requires std::invocable<F, double> && std::invocable<G, double>
class TransverseWaves { // Takes infinite string
public:
	TransverseWaves(String s_, F f_, G g_) : s(s_), f(f_), g(g_) {}
	double displacement(double x, double t) const {
		const double c = s.wavespeed();
		return f(x - c * t) + g(x + c * t);
	}
private:
	String s;
	F f; // Right-bound wave
	G g; // Left-bound wave
};

// Antiquated, keep for reference and solution checking
template<class Shape>
requires std::invocable<Shape, double>
class PluckedString { // Initial displacement, zero initial velocity
public:
	PluckedString(Air air, String s_, std::size_t n_, Shape shape_) : s(s_), n(n_), shape(shape_) {
		if (!(s.get_right() == BoundaryCondition::Fixed && s.get_left() == BoundaryCondition::Fixed))
			throw std::invalid_argument("PluckedString requires fixed ends");
		make_modes();
	}
	void make_modes() {
		const double L = s.get_length();
		for (std::size_t i = 0; i < n; ++i) {
			const double k = static_cast<double>(i + 1); // mode number 1..n
			B.push_back((2.0 / L) * simpson([&](double x) {
				return shape(x) * std::sin(k * std::numbers::pi * x / L);
			}, 0.0, L));
			A.push_back(0.0); // A_n = 0 for zero initial velocity
		}
	}
	double mode_displacement(std::size_t i, double x, double t) const {
		x = s.check_length(x);
		const double w_i = (i + 1) * std::numbers::pi * s.wavespeed() / s.get_length();
		return (A[i] * std::sin(w_i * t) + B[i] * std::cos(w_i * t) ) * std::sin( (w_i*x)/s.wavespeed());
	}
	double displacement(double x, double t) const {
		x = s.check_length(x);
		double sum = 0.0;
		for (std::size_t i = 0; i < n; ++i) sum += mode_displacement(i, x, t);
		return sum;
	}
private:
	Air air;
	String s;
	std::size_t n; // # of modes. TODO: currently an input
	Shape shape;   // Displacement function y(x,0)
	std::vector<double> A; // Mode contribution from y'
	std::vector<double> B; // Mode contribution from y
};

// Antiquated, keep for reference and solution checking
template<class Shape>
requires std::invocable<Shape, double>
class StruckString { // Zero initial displacement, initial velocity
public:
	StruckString(String s_, std::size_t n_, Shape shape_) : s(s_), n(n_), shape(shape_) {
		if (!(s.get_right() == BoundaryCondition::Fixed && s.get_left() == BoundaryCondition::Fixed))
			throw std::invalid_argument("StruckString requires fixed ends");
		make_modes();
	}
	void make_modes() {
		const double L = s.get_length();
		for (std::size_t i = 0; i < n; ++i) {
			const double k = static_cast<double>(i + 1); // mode number 1..n
			const double w_i = k * std::numbers::pi * s.wavespeed() / L;
			A.push_back((2.0 / (L * w_i)) * simpson([&](double x) {
				return shape(x) * std::sin(k * std::numbers::pi * x / L);
			}, 0.0, L));
			B.push_back(0.0); // B_n = 0 for zero initial displacement
		}
	}
	double mode_displacement(std::size_t i, double x, double t) const {
		x = s.check_length(x);
		const double w_i = (i + 1) * std::numbers::pi * s.wavespeed() / s.get_length();
		return (A[i] * std::sin(w_i * t) + B[i] * std::cos(w_i * t)) * std::sin(w_i * x / s.wavespeed());
	}
	double displacement(double x, double t) const {
		x = s.check_length(x);
		double sum = 0.0;
		for (std::size_t i = 0; i < n; ++i) sum += mode_displacement(i, x, t);
		return sum;
	}
private:
	String s;
	std::size_t n; // # of modes. TODO: currently an input
	Shape shape;   // Velocity function y'(x,0)
	std::vector<double> A; // Mode contribution from y'
	std::vector<double> B; // Mode contribution from y
};
/// Antiquated, keep for reference and solution checking
class HammerString {
public:
	HammerString(String s_, double v0_, double M_, double x0_) : s(s_), v0(v0_), M(M_), x0(x0_) {}
	double contact_displacement(double t) const;
	double displacement(double x , double t) const; // TODO: only true for long strings
private:
	String s;
	double v0; // Initial velocity [m/s]
	double M;  // Hammer mass [kg]
	double x0; // Strike point [m]
};

// class BowedString {
//     TODO
// };

} 