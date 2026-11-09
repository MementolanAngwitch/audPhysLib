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
	float m; 
};

class String {
public:
	String(float L_, float mu_, float T_, float r_, std::complex<float> E_,
	       BoundaryCondition left_ = BoundaryCondition::Fixed,
	       BoundaryCondition right_ = BoundaryCondition::Fixed)
		: L(L_), mu(mu_), T(T_), r(r_), E(E_), left(left_), right(right_) {}
	float get_length() const { return L; }
	float get_mu() const { return mu; }
	float get_rho() const { return mu / std::numbers::pi * r * r; }
	float get_tension() const { return T; }
	float wavespeed() const { return std::sqrt(T / mu); }
	float get_radius() const { return r;}
	std::complex<float> get_young_mod() const { return E;}
	std::complex<float> impedance() const { return std::sqrt(T * mu); }
	BoundaryCondition get_left() const { return left; }
	BoundaryCondition get_right() const { return right; }
	// Folds a position back onto the string and returns it.
	float check_length(float x) const {
		if (std::abs(x) > L) x = L - (x - L);
		if (x < 0) x = -x;
		return x;
	}
private:
	float L;        // length [m]
	float mu;       // linear density [kg/m]
	float T;        // tension [N]
	float r;        // radius [m]
	std::complex<float> E;   // Young's modulus [N/m^2]
	std::size_t parts = 0;   // # of divisions of the string
	BoundaryCondition left;  // negative direction
	BoundaryCondition right; // positive direction
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


// Damping decay times: A = A_0 * e^(-t/tau)
float tau_air(Air air, String s, float w, float v);

float tau_internal(String s, float w, float v);

float tau_support(String s, float w);

float tau_total(const Air& air, const String& s, float w, float v);

// sin(k+1 theta) = 2cos(theta) * sin(k theta) - sin((k-1)theta) where theta = pi x /L
void fill_shapes(double x, double L, set::vector<double>& out);
// q_n'' + (2/tau_n) q_n' + w_n^2 q_n = (2 / (mu L)) * phi_n(x0) * F(t),
// with mode shape phi_n(x) = sin(n pi x / L).
// On a fixed-fixed string
class PlayString {
public: 
	PlayString(const String& s_, const Air& air_, std::size_t n, double x_listen) 
			   : s(s_), air(air_), q(n, 0.0), qdot(n, 0.0), w(n), tau(n) {
		const float L = s.get_length();
		for (std::size_t i = 0; i < n; ++i) {
			w[i] = (i+1) * std::numbers::pi * s.wavespeed() / L; // f_n = n v/2L -> wn = f_n * 2pi = (n*v*pi)/L
			tau[i] = tau_total(air_, s_, w[i], 0.0);
		}
		fill_shapes(x_listen, L, phi_listen);
	}

	void set_contact(double x) {x_contact = x; fill_shapes(x, s.get_length(), phi_contact);};
	void pluck(float x0, float h);
	void strike(float x0, float J);
	void step(float dt, float x0, float F);
	float output() const { return sum(phi_listen); }
	double contact_displacement const { return sum(phi_contact); }

private:
	String s;
	Air air; 
	std::vector<float> tau;  // decay time [s]
	std::vector<float> q;    // modal displacement [m]
	std::vector<float> qdot; // modal velocity [m/s]
	std::vector<float> w;	 // angular frequency [rad/s]
	std::vector<float> phi_contact(n), phi_listen(n); // modes at contact and at listen
	double x_contact = 0.0;
};
// A hammer with a compressible felt tip. Force exists only while the felt is squeezed:
// F = K * compression^p  if compression > 0, else 0.
class Hammer {
public:
	Hammer(float M_, float K_, float p_, float x0_, float v0_) 
	 : M(M_), K(K_), x0(x0_), v(v0_) {};
	void step(PlayString& s, float dt);
private:
	float M; 	    // hammer mass [kg]
	float K;       // felt stiffness [n/m^p]
	float p; 	    // felt stiffness exponent 
	float x0; 	    // strike point [m]
	float v;       // hammer velocity [m/s], postive toward string
	float y = 0.0; // hammer position [m], starts touching the string at rest
};

template<class F, class G>
requires std::invocable<F, float> && std::invocable<G, float>
class TransverseWaves { // Takes infinite string
public:
	TransverseWaves(String s_, F f_, G g_) : s(s_), f(f_), g(g_) {}
	float displacement(float x, float t) const {
		const float c = s.wavespeed();
		return f(x - c * t) + g(x + c * t);
	}
private:
	String s;
	F f; // Right-bound wave
	G g; // Left-bound wave
};

// Antiquated, keep for reference and solution checking
template<class Shape>
requires std::invocable<Shape, float>
class PluckedString { // Initial displacement, zero initial velocity
public:
	PluckedString(Air air, String s_, std::size_t n_, Shape shape_) : s(s_), n(n_), shape(shape_) {
		if (!(s.get_right() == BoundaryCondition::Fixed && s.get_left() == BoundaryCondition::Fixed))
			throw std::invalid_argument("PluckedString requires fixed ends");
		make_modes();
	}
	void make_modes() {
		const float L = s.get_length();
		for (std::size_t i = 0; i < n; ++i) {
			const float k = static_cast<float>(i + 1); // mode number 1..n
			B.push_back((2.0 / L) * simpson([&](float x) {
				return shape(x) * std::sin(k * std::numbers::pi * x / L);
			}, 0.0, L));
			A.push_back(0.0); // A_n = 0 for zero initial velocity
		}
	}
	float mode_displacement(std::size_t i, float x, float t) const {
		x = s.check_length(x);
		const float w_i = (i + 1) * std::numbers::pi * s.wavespeed() / s.get_length();
		return (A[i] * std::sin(w_i * t) + B[i] * std::cos(w_i * t) ) * std::sin( (w_i*x)/s.wavespeed());
	}
	float displacement(float x, float t) const {
		x = s.check_length(x);
		float sum = 0.0;
		for (std::size_t i = 0; i < n; ++i) sum += mode_displacement(i, x, t);
		return sum;
	}
private:
	Air air;
	String s;
	std::size_t n; // # of modes. TODO: currently an input
	Shape shape;   // Displacement function y(x,0)
	std::vector<float> A; // Mode contribution from y'
	std::vector<float> B; // Mode contribution from y
};

// Antiquated, keep for reference and solution checking
template<class Shape>
requires std::invocable<Shape, float>
class StruckString { // Zero initial displacement, initial velocity
public:
	StruckString(String s_, std::size_t n_, Shape shape_) : s(s_), n(n_), shape(shape_) {
		if (!(s.get_right() == BoundaryCondition::Fixed && s.get_left() == BoundaryCondition::Fixed))
			throw std::invalid_argument("StruckString requires fixed ends");
		make_modes();
	}
	void make_modes() {
		const float L = s.get_length();
		for (std::size_t i = 0; i < n; ++i) {
			const float k = static_cast<float>(i + 1); // mode number 1..n
			const float w_i = k * std::numbers::pi * s.wavespeed() / L;
			A.push_back((2.0 / (L * w_i)) * simpson([&](float x) {
				return shape(x) * std::sin(k * std::numbers::pi * x / L);
			}, 0.0, L));
			B.push_back(0.0); // B_n = 0 for zero initial displacement
		}
	}
	float mode_displacement(std::size_t i, float x, float t) const {
		x = s.check_length(x);
		const float w_i = (i + 1) * std::numbers::pi * s.wavespeed() / s.get_length();
		return (A[i] * std::sin(w_i * t) + B[i] * std::cos(w_i * t)) * std::sin(w_i * x / s.wavespeed());
	}
	float displacement(float x, float t) const {
		x = s.check_length(x);
		float sum = 0.0;
		for (std::size_t i = 0; i < n; ++i) sum += mode_displacement(i, x, t);
		return sum;
	}
private:
	String s;
	std::size_t n; // # of modes. TODO: currently an input
	Shape shape;   // Velocity function y'(x,0)
	std::vector<float> A; // Mode contribution from y'
	std::vector<float> B; // Mode contribution from y
};
/// Antiquated, keep for reference and solution checking
class HammerString {
public:
	HammerString(String s_, float v0_, float M_, float x0_) : s(s_), v0(v0_), M(M_), x0(x0_) {}
	float contact_displacement(float t) const;
	float displacement(float x , float t) const; // TODO: only true for long strings
private:
	String s;
	float v0; // Initial velocity [m/s]
	float M;  // Hammer mass [kg]
	float x0; // Strike point [m]
};

// class BowedString {
//     TODO
// };

} 