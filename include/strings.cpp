#include "strings.hpp"

namespace audphys::strings {
void fill_shapes(float x, float L, std::vector<float>& out) { //loses accuracy at high mode counts
	const float theta = std::numbers::pi * x/L;
	const float c = 2.0 * std::cos(theta);
	float prev = 0.0;
	float curr = std::sin(theta);
	for (std::size_t i = 0; i < out.size(); ++i) {
		out[i] = curr;
		const double next = c * curr - prev;
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
	return s.get_length() / (2.0 * s.get_tension() * G);
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
		const float B = 2.0 * h * L * L / (k * k * x_contact * (L - x_contact)) * phi_contact[i];
		q[i]      += B;
		q_prev[i] += B * back[i];
	}
}

void PlayString::strike(float x0, float J){
	set_contact(x0); // FIXED: x0 was ignored
	const float g = 2.0 / (s.get_mu() * s.get_length());
	for (std::size_t i = 0; i < q.size(); ++i) {
		q_prev[i] -= kick[i] * g * J * phi_contact[i];  // q unchanged: string still flat
	}
}

void PlayString::step(float, float F) {
	for (std::size_t i = 0; i < q.size(); ++i) {
		const double next = c1[i] * q[i] - c2[i] * q_prev[i] + b[i] * phi_contact[i] * F;
		q_prev[i] = q[i];
		q[i] = next;
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

}