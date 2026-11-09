#include "strings.hpp"

namespace audphys::strings {
void fill_shapes(float x, float L, std::vector<float>& out) {
	const float theta = std::numbers::pi * x/L;
	const float c = 2.0 * std::cos(theta);
	float prev = 0.0;
	float curr = std::sin(theta);
	for (std::size_t i = 0; i < out.size(); ++i) {
		out[i] = curr;
		const float next = c * curr - prev;
		prev = curr;
		curr = next;
	}
}

// Damping decay times: A = A_0 * e^(-t/tau)
float tau_air(const Air& air, const String& s, float w) {
	const float rho = s.get_rho();
	const float M = (s.get_radius()/2.0) * std::sqrt(w/air.get_eta());
	return rho/(air.get_rho() * w) * (2.0*M*M)/(2.0*std::sqrt(2.0)*M+1.0);
}

float tau_internal(const String& s, float w) {
	const float E1 = s.get_young_mod().real();
	const float E2 = s.get_young_mod().imag();
	return (2.0 / w) * E1/E2;
}

float tau_support(const String& s, float G) {
	const float f1 = s.wavespeed() / (2.0 * s.get_length());
	return 1.0 / (8.0 * s.get_mu() * s.get_length() * f1 * f1 * G);
}

float tau_total(const Air& air, const String& s, float w, float G) {
	return 1.0 / (1.0/tau_air(air,s,w) + 1.0/tau_internal(s,w) + 1.0/tau_support(s,G));
}

// Pluck at the contact point, pulled to h
void PlayString::pluck(float x0, float h) {
	set_contact(x0);
	const float L = s.get_length();
	for (std::size_t i = 0; i < q.size(); ++i) {
		const float k = (i + 1) * std::numbers::pi;
		q[i] += 2.0 * h * L * L / (k * k * x_contact * (L - x_contact)) * phi_contact[i]; 
	}
}

void PlayString::strike(float x0, float J){
	set_contact(x0); // FIXED: x0 was ignored
	const float g = 2.0 / (s.get_mu() * s.get_length());
	for (std::size_t i = 0; i < qdot.size(); ++i) {
		qdot[i] += g * J * phi_contact[i];
	}
}

void PlayString::step(float dt, float F) {
	const float g = 2.0 / (s.get_mu() * s.get_length());
	for (std::size_t i = 0; i < q.size(); ++i) {
		const float a = g * F * phi_contact[i] - 2.0 * qdot[i] / tau[i] - w[i] * w[i] * q[i];
		qdot[i] += a * dt;
		q[i] += qdot[i] * dt;
	}
}

float Hammer::step(PlayString& s, float dt) {
	const float compression = y - s.contact_displacement();
	const float F = compression > 0.0 ? K * std::pow(compression, p) : 0.0;
	v -= F / M * dt; //string pushes back on the hammer
	y += v * dt; 
	s.step(dt,F); // Hammer pushes on the string
	return F;
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