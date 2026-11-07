#pragma once
#include <cmath>
#include <vector>
#include <concepts>
#include "include/tools.hpp"
namespace audphys::strings{

enum class BoundaryCondition { Fixed, Free };

class String {
public:
 String(double L_, double mu_, double T_, BoundaryCondition left_ = BoundaryCondition::Fixed, BoundaryCondition right_ = BoundaryCondition::Fixed) :
 		L(L_), mu(mu_), T(T_),  left(left_), right(right_) {};
 double get_length() const {return L;};
 double get_mu() const {return mu;};
 double get_tension() const {return T;};
 double wavespeed() const {return std::sqrt(T/mu);};
 BoundaryCondition get_left() const {return left;};
 BoundaryCondition get_right() const {return right;};
private: 
 double L; // length [m]
 double mu; // liner density [kg/m]
 double T; // Tension [N]
 std::size_t parts; // # of divisions of the string
 BoundaryCondition left;
 BoundaryCondition right; // let left be negative direction, right be positive direction
};

template<class F, class G>
requires std::invocable<F, double> && std::invocable<G,double>
class TransverseWaves { // Takes infintie string
public: 
 TransverseWaves(String s_, F f_, G g_) : s(s_), f(f_), g(g_) {};
 double displacement(double x, double t) const {
 	const double c = s.wavespeed();
 	return f(x-c*t) + g(x+c*t);	
 }
private:
 String s;
 F f; //Right Bound Wave
 G g; //Left Bound Wave
};

class PluckingString {
 public:
 	PluckedString(String s_, double n_, double x0_, double v0_) : s(s_), n(n_), x0(x0_), v0(v0_) {make_modes(); };
 	{ if (!(s.get_right() == BoundaryCondition::Fixed && s.get_left() == BoundaryCondition::Fixed)) throw; };
 	void make_modes();
 	double mode_displacement(std::size_t n, double x, double t); 
 	double displacement(double x, double t);
 private:
 	std::vector<double> A;
 	std::vector<double> B;
 	double n;
 	String s;
 	x0;
 	v0;
};

}