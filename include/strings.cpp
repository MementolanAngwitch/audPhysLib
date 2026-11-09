#include "strings.hpp"

namespace audphys::strings {
void fill_shapes(double x, double L, set::vector<double>& out) {
	const double theta = std::numbers::pi * x/L;
	const double c = 2.0 * std::cos(theta);
	double prev = 0.0;
	double curr = std::sin(theta);
	for (std::size_t i = 0; i < out.size(); ++i) {
		out[i] = curr;
		const double next = c * curr - prev;
		prev = curr;
		curr = next;
	}
}

// Damping decay times: A = A_0 * e^(-t/tau)
float tau_air(Air air, String s, float w, float v) {
	const float rho = s.get_mu() * std::numbers::pi * s.get_radius() * s.get_radius();
	const float M = (s.get_radius()/2.0) * std::sqrt(2.0*std::numbers::pi*w / air.get_eta());
	return rho/(2.0 * std::numbers::pi * air.get_rho() * w) * (2.0*M*M)/(2.0 * std::sqrt(2.0) * M + 1.0);
} 

float tau_internal(String s, float w, float v) {
	const float E1 = s.get_young_mod().real();
	const float E2 = s.get_young_mod().imag();
	return (1.0 / (std::numbers::pi * w)) * E1/E2;
}

float tau_support(String s, float w) {
	const float G = std::pow(s.impedance(), -1).real(); 
	return std::pow(8.0 * s.get_mu() * s.get_length() * w * w * G,-1);
}

float tau_total(Air air, String s, float w, float v) {
	return std::pow( 1.0/tau_air(air,s,w,v)  + 1.0/tau_internal(s,w,v) + 1.0/tau_support(s,w),-1);
}

// Pluck at the contact point, pulled to h
void PlayString::pluck(double x0, double h) {
	const double L = s.get_length();
	for (std::size_t i = 0; i < q.size(); ++i) {
		const double k = (i + 1) * std::numbers::pi;
		q[i] += 2.0 * h * L * L / (k * k * x_contact * (L - x_contact)) * phi_contact[i]; 
	}
}

void PlayString::strike(float x0, float J){
	const float g = 2.0 / (s.get_mu() * s.get_length());
	for (std::size_t i = 0; i < qdot.size(); ++i) {
		qdot[i] += g * J * phi_contact[i];
	}
}

void PlayString::step(double dt, double F = 0.0) {
	const double g = 2.0 / (s.get_mu() * s.get_length());
	for (std::size_t i = 0; i < q.size(); ++i) {
		const double a = g * F * phi_contact[i] - 2.0 * qdot[i] / tau[i] - w[i] * w[i] * q[i];
		qdot[i] += a * dt;
		q[i] += qdot[i] * dt;
	}
}

void Hammer::step(PlayString& s, float dt) {
	const float compression = y - s.contact_displacement();
	const float F = compression > 0.0 ? K * std::pow(compression, p) : 0.0;
	v -= F / M * dt; //string pushes back on the hammer
	y += v * dt; 
	s.step(dt,x0,F); // Hammer pushes on the string
}


float HammerString::contact_displacement(float t) const {
	const float tau = M * s.wavespeed() / (2.0 * s.get_tension());
	return v0 * tau * (1.0 - std::exp(-t / tau));
}

float HammerString::displacement(float x, float t) const {
	x = s.check_length(x);
	return contact_displacement(t - ((x - x0) / s.wavespeed()));
}

}