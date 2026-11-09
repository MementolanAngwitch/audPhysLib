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
double tau_air(Air air, String s, double w, double v) {
	const double rho = s.get_mu() * std::numbers::pi * s.get_radius() * s.get_radius();
	const double M = (s.get_radius()/2.0) * std::sqrt(2.0*std::numbers::pi*w / air.get_eta());
	return rho/(2.0 * std::numbers::pi * air.get_rho() * w) * (2.0*M*M)/(2.0 * std::sqrt(2.0) * M + 1.0);
} 

double tau_internal(String s, double w, double v) {
	const double E1 = s.get_young_mod().real();
	const double E2 = s.get_young_mod().imag();
	return (1.0 / (std::numbers::pi * w)) * E1/E2;
}

double tau_support(String s, double w) {
	const double G = std::pow(s.impedance(), -1).real(); 
	return std::pow(8.0 * s.get_mu() * s.get_length() * w * w * G,-1);
}

double tau_total(Air air, String s, double w, double v) {
	return std::pow( 1.0/tau_air(air,s,w,v)  + 1.0/tau_internal(s,w,v) + 1.0/tau_support(s,w),-1);
}

double PlayString::shape(std::size_t i, double x) const {
	return std::sin((i+1) * std::numbers::pi * x / s.get_length());
}

void PlayString::strike(double x0, double J){
	const double g = 2.0 / (s.get_mu() * s.get_length());
	for (std::size_t i = 0; i < qdot.size(); ++i) {
		qdot[i] += g * J * shape(i, x0);
	};
}

double PlayString::displacement(double x) const {
	double y = 0.0;
	for (std::size_t i = 0; i < q.size(); ++i) {
		y += q[i] * shape(i,x);
	}
	return y;
}

void PlayString::step(double dt, double x0 = 0.0, double F = 0.0) {
	const double g = 2.0 / (s.get_mu() * s.get_length());
	for (std::size_t i = 0; i < q.size(); ++i) {
		const double a = (g * F * shape(i,x0)) - (2.0 * qdot[i]/ tau[i]) - (w[i] * w[i] * q[i]); // force * damping * restoring force
		qdot[i] += a * dt;
		q[i] += qdot[i] * dt;
	}
}

void Hammer::step(PlayString& s, double dt) {
	const double compression = y - s.displacement(x0);
	const double F = compression > 0.0 ? K * std::pow(compression, p) : 0.0;
	v -= F / M * dt; //string pushes back on the hammer
	y += v * dt; 
	s.step(dt,x0,F); // Hammer pushes on the string
}

}