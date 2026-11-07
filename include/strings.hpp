#pragma once
#include <cmath>
#include <cstddef>
#include <numbers>
#include <stdexcept>
#include <vector>
#include <concepts>
#include "tools.hpp"
namespace audphys::strings {

enum class BoundaryCondition { Fixed, Free };

class String {
public:
	String(double L_, double mu_, double T_,
	       BoundaryCondition left_ = BoundaryCondition::Fixed,
	       BoundaryCondition right_ = BoundaryCondition::Fixed)
		: L(L_), mu(mu_), T(T_), left(left_), right(right_) {}
	double get_length() const { return L; }
	double get_mu() const { return mu; }
	double get_tension() const { return T; }
	double wavespeed() const { return std::sqrt(T / mu); }
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
	double E = 0.0;  // Young's modulus [N/m^2]
	std::size_t parts = 0; // # of divisions of the string
	BoundaryCondition left;  // negative direction
	BoundaryCondition right; // positive direction
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

template<class Shape>
requires std::invocable<Shape, double>
class PluckedString { // Initial displacement, zero initial velocity
public:
	PluckedString(String s_, std::size_t n_, Shape shape_) : s(s_), n(n_), shape(shape_) {
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
	Shape shape;   // Displacement function y(x,0)
	std::vector<double> A; // Mode contribution from y'
	std::vector<double> B; // Mode contribution from y
};

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

class HammerString {
public:
	HammerString(String s_, double v0_, double M_, double x0_) : s(s_), v0(v0_), M(M_), x0(x0_) {}
	double contact_displacement(double t) const;
	double displacement(double t, double x) const; // TODO: only true for long strings
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