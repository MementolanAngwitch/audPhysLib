#pragma once
#include <cmath>
#include <numbers>
#include <cassert>
#include <stdexcept>
#include <complex>
#include <vector>

namespace audphys::oscillators {
class Oscillator{
public:
 Oscillator(double m_, double K_, double R_, double F_) : m(m_), K(K_), R(R_), F(F_) {};
 const double get_mass(){return m;};
 const double get_K(){return K;};
private: 
 double m; // mass [kg]
 double K; // spring constant [k/m]
 double R = 0.0; // Mechnical Resistance
 double F = 0.0; // Force [N]
};
	
class WaveMotion {
public:
// ---- Constructors
 WaveMotion(Oscillator o,double x0_, double v0_, double a0_, double w0_);

// ---- Eq of motion
 double period(); // 1/f [s]
 double natural_frequency(); // w0 = \sqrt(K/m) [rad/s]
 double displacement(); // x = Acos(w0 t + phi) [m]
 double velocity(); // v = -w0 A sin(w0t + phi) [m/s]
 double acceleration(); // a = -w0^2 A cos(w0t + phi) [m/s^2]
 double amplitude(); // A = sqrt(x0^2 + (v0/w0)^2) [m]
 double phi(); // tan-1(-v0/w0x0) [rad]

// ---- Superposition
 WaveMotion superposition(const std::vector<WaveMotion>&); // For each, check for same/different frequency,
 														   // Use A = sqrt( (sum(A_n*cos(phi_n)))^2 + ((sum(A_n*sin(phi_n)))^2)
 														   // Use tan(phi) = (sum(A_n*sin(phi_n)))/(sum(A_n*cos(phi_n)))
 														   // For different frequncies write w_i+1 = w_i + delta(w)_i
 WaveMotion operator+(const WaveMotion&);
 WaveMotion operator*(const double);
 WaveMotion operator-(const WaveMotion&);

// ---- Energy 
 double potential_energy(); //1/2 Kx^2 [J]
 double kinetic_energy(); // 1/2 mU^2 [J]
 double total_energy(); // 

// ---- Visualization
private:
 struct InitialConditions { double x0; double v0; double a0; double w0; };
 Oscillator osc; 
 double A;
 double Phi;
};

inline double toAngular(double f) { return 2.0 * std::numbers::pi * f; }; //[hz to rad]
inline double fromAngular(double w) { return w / (2.0 * std::numbers::pi); }; // [rad to hz]

} //namespace audphys