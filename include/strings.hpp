#pragma once
#include <cmath>
#include <vector>

namespace audphyslib::strings{

enum class BoundaryCondition { Fixed, Free };

class String() {
 String(double L_, double mu_, ,double T_, BoundaryCondition left_, BoundaryCondition right_) :
 		L(L_), mu(mu_), T(T_),  left(left_), right(right_) {};
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
class TransverseWaves() {
public: 
 double displacement(F f, G g, double x, double t) {
 	c = std::sqrt(c.T / c.mu);
 	return f(c*t - x) + g (c*t + x);
 }
private:
 String s;
}

}