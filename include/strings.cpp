#include "strings.hpp"

namespace audphys::strings {

double HammerString::contact_displacement(double t) const {
	const double tau = M * s.wavespeed() / (2.0 * s.get_tension());
	return v0 * tau * (1.0 - std::exp(-t / tau));
}

double HammerString::displacement(double x, double t) const {
	x = s.check_length(x);
	return contact_displacement(t - ((x - x0) / s.wavespeed()));
}

// Damping decay times: A = A_0 * e^(-t/tau)
double tau_air(Air air, String s, double f, double v) {
	const double rho = s.get_mu() * std::numbers::pi * s.get_radius() * s.get_radius();
	const double M = (s.get_radius()/2.0) * std::sqrt(2.0*std::numbers::pi*f / air.get_eta());
	return rho/(2.0 * std::numbers::pi * air.get_rho() * f) * (2.0*M*M)/(2.0 * std::sqrt(2.0) * M + 1.0) * 2.0 * std::numbers::pi;
} 

double tau_internal(String s, double f, double v) {
	const double E1 = s.get_young_mod().real();
	const double E2 = s.get_young_mod().imag();
	return (1.0 / (std::numbers::pi * f)) * E1/E2  * 2.0 * std::numbers::pi;
}

double tau_support(String s, double f) {
	const double G = std::pow(s.impedance(), -1).real(); 
	return std::pow(8.0 * s.get_mu() * s.get_length() * f * f * G,-1);
}

double tau_total(Air air, String s, double f, double v) {
	return std::pow( 1.0/tau_air(air,s,f,v)  + 1.0/tau_internal(s,f,v) + 1.0/tau_support(s,f),-1) * 2.0 * std::numbers::pi;
}

}