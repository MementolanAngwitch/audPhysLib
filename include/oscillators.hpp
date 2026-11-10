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
 friend class Oscillator;
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
 friend class Oscillator;
public:
 HelmholtzResonator(double rho_, double S_, double L_, double c_) : rho(rho_), S(S_), L(L_), c(c_) {}; 
private:
 double rho; // air density [kg/m^3]
 double S;   // cross-section area of neck [m^2]
 double L;   // neck length [m]
 double V;   // cavity volume [m^3]
 double c = 343.3; // Speed of Sound
};

class Impulse { // we only consider impulse force such as a pluck or strike, wind instruments will be modelled by air speed
 friend class WaveMotion;
public:
 Impulse(double F_, double t_, double w_) : F(F_), t(t_), w(w_) {};
private: 
 double F=0.0; // Force at t=0 [N]
 double t=0.0; // Time of force application [s]
 double w=0.0; // Sinusoidal wave [rad]
};

class Oscillator{
 friend class WaveMotion;
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
 WaveMotion(Oscillator o, Impulse F_, double x0_, double v0_, double a0_);

// ---- Eq of motion
 double period(double) const; // 1/f [s]
 //double natural_frequency() const; // w0 = \sqrt(K/m) [rad] if R > 0, called damped natural frequency
 double displacement(double) const; // x = Acos(w0 t + phi) [m] function takes only oscillator and calls best helper function
 double velocity(double) const; // v = -w0 A sin(w0t + phi) [m/s] function takes only oscillator and calls best helper function
 double acceleration(double) const; // a = -w0^2 A cos(w0t + phi) [m/s^2] function takes only oscillator and calls best helper function

 double get_A() const; //For time dependant amplitudes, return A0
 double get_phi() const; 

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
 double mechanical_reactance() const;
 std::complex<double> mechanical_impedance() const;
};


}