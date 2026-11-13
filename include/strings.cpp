#include "strings.hpp"

namespace audphys::strings {

// ---- Mode frequencies and decay ---------------------------------------------

std::complex<float> support_shift(const String& s, std::size_t n) {
	const float w0 = n * s.fundamental();  // rigid-end frequency; first order is enough
	const std::complex<float> j(0.0f, 1.0f);
	const std::complex<float> d = j * s.impedance()
	                            * (s.get_left().admittance(w0) + s.get_right().admittance(w0));
	// Written as !(x < limit) so a NaN (support at its own resonance) also throws.
	if (!(std::abs(d) < 0.5f))
		throw std::invalid_argument("support_shift: support too soft for the perturbation model");
	return d;
}

float mode_frequency(const String& s, std::size_t n) {
	const float shift = support_shift(s, n).real();
	return s.wavespeed() * (n * std::numbers::pi_v<float> + shift) / s.get_length();
}

float tau_air(const Air& air, const String& s, float w) {
	const float M = (s.get_radius() / 2.0f) * std::sqrt(w / air.get_eta());
	return s.get_rho() / (air.get_rho() * w) * (2.0f * M * M) / (2.0f * std::numbers::sqrt2_v<float> * M + 1.0f);
}

float tau_internal(const String& s, float w) {
	const std::complex<float> E = s.get_young_mod();
	return (2.0f / w) * E.real() / E.imag();
}

float tau_support(const String& s, std::size_t n) {
	return s.get_length() / (s.wavespeed() * support_shift(s, n).imag());
}

float tau_total(const Air& air, const String& s, std::size_t n) {
	const float w = mode_frequency(s, n);
	// Rigid or loss-free pieces give 1/inf = 0 and drop out.
	return 1.0f / (1.0f / tau_air(air, s, w) + 1.0f / tau_internal(s, w) + 1.0f / tau_support(s, n));
}

void fill_shapes(float x, float L, std::vector<float>& out) {  // loses accuracy at very high mode counts
	const float theta = std::numbers::pi_v<float> * x / L;
	const float c = 2.0f * std::cos(theta);
	float prev = 0.0f;
	float curr = std::sin(theta);
	for (std::size_t i = 0; i < out.size(); ++i) {
		out[i] = curr;
		const float next = c * curr - prev;
		prev = curr;
		curr = next;
	}
}

// ---- Modal simulation ------------------------------------------------------------

PlayString::PlayString(const String& s_, const Air& air, float w_max, float dt_, float x_listen)
	: s(s_), dt(dt_) {
	const float L = s.get_length();
	const float nyquist = std::numbers::pi_v<float> / dt; //highest frequency that can be reproduced without error
	const std::size_t n = static_cast<std::size_t>(std::min(w_max, nyquist) / s.fundamental());

	w.assign(n, 0.0f); tau.assign(n, 0.0f);
	q.assign(n, 0.0); q_prev.assign(n, 0.0);
	c1.assign(n, 0.0); c2.assign(n, 0.0); b.assign(n, 0.0);
	back.assign(n, 0.0); kick.assign(n, 0.0);
	phi_contact.assign(n, 0.0f); phi_listen.assign(n, 0.0f);

	const double g = 2.0 / (s.get_mu() * L);  // force -> modal acceleration
	for (std::size_t i = 0; i < n; ++i) {
		w[i]   = mode_frequency(s, i + 1);
		tau[i] = tau_total(air, s, i + 1); // i + 1 by convention so w[1] is the fundamental frequency
		const double a  = 1.0 / tau[i];                          // decay rate
		const double wd = std::sqrt(double(w[i]) * w[i] - a * a); // damped frequency
		const double R  = std::exp(-a * dt);
		const double th = wd * dt;

		//q[k+1] = 2R cos(w_d dt) q[k] - R^2 q[k-1] + b F,   R = e^(-dt/tau)
		c1[i]   = 2.0 * R * std::cos(th);
		c2[i]   = R * R;
		b[i]    = g * dt * R * std::sin(th) / wd;                                  // times phi_contact * F
		back[i] = std::exp(a * dt) * (std::cos(th) - a * std::sin(th) / wd);       // q[-1] / q[0], at rest
		kick[i] = std::exp(a * dt) * std::sin(th) / wd;                            // q[-1] per unit velocity
	}
	fill_shapes(x_listen, L, phi_listen);
}

void PlayString::set_contact(float x) {
	if (x < 0.0f || x > s.get_length()) throw std::invalid_argument("set_contact: x outside [0, L]");
	x_contact = x;
	fill_shapes(x, s.get_length(), phi_contact);
}

void PlayString::pluck(float x0, float h) {
	const float L = s.get_length();
	if (!(x0 > 0.0f && x0 < L)) throw std::invalid_argument("pluck: x0 must be strictly inside (0, L)");
	set_contact(x0);
	for (std::size_t i = 0; i < q.size(); ++i) {
		const double k = (i + 1) * std::numbers::pi;
		const double B = 2.0 * h * L * L / (k * k * x0 * (L - x0)) * phi_contact[i];
		q[i]      += B;
		q_prev[i] += B * back[i];
	}
}

void PlayString::strike(float x0, float J) {
	set_contact(x0);
	const double g = 2.0 / (s.get_mu() * s.get_length());
	for (std::size_t i = 0; i < q.size(); ++i) {
		q_prev[i] -= kick[i] * g * J * phi_contact[i];  // q unchanged: the string is still flat
	}
}

void PlayString::step(float F) {
	for (std::size_t i = 0; i < q.size(); ++i) {
		const double next = c1[i] * q[i] - c2[i] * q_prev[i] + b[i] * phi_contact[i] * F;
		q_prev[i] = q[i];
		q[i] = next;
	}
}

// ---- Exciters ---------------------------------------------------------------------

float Hammer::step(PlayString& s) {
	if (s.contact() != x0) s.set_contact(x0);
	const float dt = s.timestep();
	const float compression = y - static_cast<float>(s.contact_displacement());
	const float F = compression > 0.0f ? K * std::pow(compression, p) : 0.0f;
	v -= F / M * dt;  // string pushes back on the hammer
	y += v * dt;
	s.step(F);        // hammer pushes on the string
	return F;
}

} // namespace audphys::strings