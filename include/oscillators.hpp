#pragma once
#include <cmath>
#include <numbers>
#include <cassert>
#include <stdexcept>
#include <complex>
#include <vector>


namespace audphys::oscillators {

class EnclosedAir{
// gamma = 1.4 for air
// K = gamma * Pa * S / L
public:
 EnclosedAir(double Pa_, double m_, double S_, double L_, double gamma_) : Pa(Pa_), m(m_), S(S_), L(L_), gamma(gamma_) {};
private:
 double Pa; //atmospheric pressure [N/m^2]
 double m; // mass [kg]
 double S; // cross section area [m]
 double L; // length [m^2]
 double gamma = 1.4;
};

class HelmholtzResonator {
// m = rho sL
// K = rho S^2 c^2 / V
public:
 HelmholtzResonator(double rho_, double S_, double L_, double c_) : rho(rho_), S(S_), L(L_), c(c_) {}; 
private:
 double rho; // pressure in container [N/m^2]
 double S; //Cross section area of neck [m^2]
 double L;
 double c = 343.3; // Speed of Sound
};

class Impulse { // we only consider impulse force such as a pluck or strike, wind instruments will be modelled by air speed
 Impulse(double F_, double t_, double w_) : F(F_), t(t_), w(w_) {};
private: 
 double F=0.0; // Force at t=0 [N]
 double t=0.0; // Time of force application [s]
 double w=0.0; // Sinusoidal wave [rad]
};

class Oscillator{
public:
 Oscillator(double m_, double K_, double R_) : m(m_), K(K_), R(R_) {};
 Oscillator(EnclosedAir a);
 Oscillator(HelmholtzResonator r);
 double get_mass() const {return m;};
 double get_K() const {return K;};
private: 
 double m; // mass [kg]
 double K; // spring constant [N/m]
 double R = 0.0; // Mechnical Resistance
};
	
class WaveMotion {
public:
// ---- Constructors
 WaveMotion(Oscillator o,Impulse F_, double x0_, double v0_, double a0_) : osc(o), F(F_), x0(x0_), v0(v0_), a0(a0_) {};

// ---- Eq of motion
 double period(double) const; // 1/f [s]
 //double natural_frequency() const; // w0 = \sqrt(K/m) [rad] if R > 0, called damped natural frequency
 double displacement(double) const; // x = Acos(w0 t + phi) [m] function takes only oscillator and calls best helper function
 double velocity(double) const; // v = -w0 A sin(w0t + phi) [m/s] function takes only oscillator and calls best helper function
 double acceleration(double); // a = -w0^2 A cos(w0t + phi) [m/s^2] function takes only oscillator and calls best helper function

 double get_A() const; //For time dependant amplitudes, return A0
 double get_phi() const; 


// ---- Superposition
 WaveMotion superposition(const std::vector<WaveMotion>&); // For each, check for same/different frequency,
 														   // Use A = sqrt( (sum(A_n*cos(phi_n)))^2 + ((sum(A_n*sin(phi_n)))^2)
 														   // Use tan(phi) = (sum(A_n*sin(phi_n)))/(sum(A_n*cos(phi_n)))
 														   // For different frequncies write w_i+1 = w_i + delta(w)_i
 WaveMotion operator+(const WaveMotion&);
 WaveMotion operator*(const double);
 WaveMotion operator-(const WaveMotion&);

// ---- Energy 
 double potential_energy() const; //1/2 Kx^2 [J]
 double kinetic_energy() const; // 1/2 mU^2 [J]
 double total_energy() const; // PE + KE constant for free oscillator [J]
 							  // For damped oscillator rate of energy loss is -2* alpha* m * v^2 
 							  // Work of forced oscillator w = F(t)*d

// ---- Visualization
private:
 // Input variables
 double x0; //[m]
 double v0; //[m/s]
 double a0; //[m/s^2]
 Impulse F;
 Oscillator osc; 
 // Oscillator types
 bool damped = false;
 bool forced = false;
 // Derived variables
 double w0 = 0.0;  // w0 = \sqrt(K/m) [rad] if R > 0, called damped natural frequency
 double A = 0.0;   // A = sqrt(x0^2 + (v0/w0)^2) [m]
 double phi = 0.0; // tan-1(-v0/w0x0) [rad]
//---- Damped Helper Functions
 double alpha() const; // R / 2m
 double damped_displacement_helper() const; // A = e^(-alpha*t) cos(wd * t + phi), 
    										// For critically damped, alpha = w0, xc = x0 *(1+alpha*t) e^(-alpha*t)
//---- Forced Helper Functions
 double forced_displacement_helper() const;
 double quality_factor() const; // Q = omega_0 / 2 alpha() where w0 is undamped nat frequency
//---- Impedance and reactance
 std::complex<double> mechnaical_reactance() const;
 std::complex<double> mechanical_impedance() const;
};

inline double toAngular(double f) { return 2.0 * std::numbers::pi * f; }; //[hz to rad]
inline double fromAngular(double w) { return w / (2.0 * std::numbers::pi); }; // [rad to hz]

} 