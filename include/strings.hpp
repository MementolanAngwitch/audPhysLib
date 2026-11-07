#pragma once
#include <cmath>
#include <vector>
#include <concepts>
#include "tools.hpp"
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

template<class Shape> 
requires std::invocable<Shape, double> 
class PluckedString {
 public:
 	PluckedString(String s_, double n_, Shape shape_) : s(s_), n(n_), shape(shape_){ 
 	 	if (!(s.get_right() == BoundaryCondition::Fixed && s.get_left() == BoundaryCondition::Fixed)) throw; 
    	make_modes();
    }
 	void make_modes(double x) {
		for (std::size_t i = 0; i < n; ++i) {
			const double w_i = i * std::numbers::pi * s.wavespeed() / s.get_length(); 
			A.push_back( (2.0 / (w_i * s.L)) * simpson([&]() { Shape(x) * std::sin((i*std::numbers::pi*x)/s.L);}, 0, s.L) );
			B.push_back(0.0); //B_n = 0 for string starting at rest
		}
	}
 	double mode_displacement(std::size_t i, double x, double t) {
		const double w_i = i * std::numbers::pi * s.wavespeed() / s.get_length(); 
		return ( A[i] * std::sin(w_i * t) + B[i] * std::cos(w_i * t) ) * std::sin((w_i*x)/s.wavespeed());
	}

 	double displacement(double x, double t) {
		double sum = 0.0;
		for (std::size_t i = 0; i < n; ++i) {
			sum += mode_displacement(i, x, t);
		}
		return sum;
	}

 private:
 	std::vector<double> A;
 	std::vector<double> B;
 	double n;
 	String s;
 	Shape shape;
};

}