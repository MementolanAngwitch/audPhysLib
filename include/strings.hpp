#pragma once
#include <cmath>
#include <complex>
#include <cstddef>
#include <numbers>
#include <stdexcept>
#include <vector>
#include "tools.hpp"

// TODO: Route class Input into strings, final syntax should look something like
// gen(String, Air, Input) -> std::vector<float> with size t/dt.
namespace audphys::strings {

// ---- Shared types -----------------------------------------------------------

// One end of a string, modelled as a mass-spring-damper the string pulls on.
// Mechanical impedance  Z(w) = R + j(w m - K/w)   [N s/m]
// Admittance            Y(w) = 1/Z                [m/(N s)]
// Support::rigid() = perfectly fixed end (Y = 0).
class Support {
public:
	static Support rigid() { return Support(); }
	Support(float m_, float K_, float R_) : m(m_), K(K_), R(R_), fixed(false) {}
	bool is_rigid() const { return fixed; }
	std::complex<float> admittance(float w) const {
		if (fixed) return {0.0f, 0.0f};
		return 1.0f / std::complex<float>(R, w * m - K / w);
	}
private:
	float m = 0.0f;  // moving mass of the support [kg]
	float K = 0.0f;  // stiffness [N/m]
	float R = 0.0f;  // mechanical resistance [N s/m]
	bool fixed = true;
	Support() = default;  // private: use Support::rigid()
};

class String {
public:
	String(float L_, float mu_, float T_, float r_, std::complex<float> E_,
	       Support left_ = Support::rigid(), Support right_ = Support::rigid())
		: L(L_), mu(mu_), T(T_), r(r_), E(E_), left(left_), right(right_) {}
	float get_length() const { return L; }
	float get_mu() const { return mu; }
	float get_rho() const { return mu / (std::numbers::pi_v<float> * r * r); }
	float get_tension() const { return T; }
	float get_radius() const { return r; }
	float wavespeed() const { return std::sqrt(T / mu); }
	float impedance() const { return std::sqrt(T * mu); }                          // Z0 [N s/m]
	float fundamental() const { return std::numbers::pi_v<float> * wavespeed() / L; } // rigid-end w1 [rad/s]
	std::complex<float> get_young_mod() const { return E; }
	const Support& get_left() const { return left; }
	const Support& get_right() const { return right; }
	// Mirrors any position back onto [0, L] (reflection at both ends).
	float check_length(float x) const {
		x = std::fmod(std::abs(x), 2.0f * L);
		return x > L ? 2.0f * L - x : x;
	}
private:
	float L;                // length [m]
	float mu;               // linear density [kg/m]
	float T;                // tension [N]
	float r;                // radius [m]
	std::complex<float> E;  // Young's modulus, imaginary part = internal loss [N/m^2]
	Support left;           // end at x = 0
	Support right;          // end at x = L
};

class Air {
public:
	Air() = default;
	Air(float rho_, float eta_) : rho(rho_), eta(eta_) {}
	float get_rho() const { return rho; }
	float get_eta() const { return eta; }
private:
	float rho = 1.2f;     // air density [kg/m^3]
	float eta = 1.5e-5f;  // kinematic viscosity [m^2/s]
};

// ---- Mode frequencies and decay ---------------------------------------------

// Each mode n has wavenumber k_n with k_n L = n pi + delta_n. For supports much stiffer
// or heavier than the string (|delta_n| << 1), each end contributes
//     delta = j Z0 Y(w_n)
// The real part shifts the pitch; the imaginary part is energy lost into the support.
// Mode shapes are kept as sin(n pi x / L), which is correct to the same order.

// Total delta_n from both ends
std::complex<float> support_shift(const String& s, std::size_t n);

// w_n = c (n pi + Re delta_n) / L  [rad/s]
float mode_frequency(const String& s, std::size_t n);

// Decay times for amplitude A = A0 e^(-t/tau)
float tau_air(const Air& air, const String& s, float w);
float tau_internal(const String& s, float w);
float tau_support(const String& s, std::size_t n);  // L / (c Im delta_n)
float tau_total(const Air& air, const String& s, std::size_t n);

// sin((k+1) theta) = 2 cos(theta) sin(k theta) - sin((k-1) theta), theta = pi x / L
void fill_shapes(float x, float L, std::vector<float>& out);

// ---- Modal simulation ------------------------------------------------------------

// q_n'' + (2/tau_n) q_n' + w_n^2 q_n = (2 / (mu L)) phi_n(x0) F(t),  phi_n(x) = sin(n pi x / L)
// Advanced by recurrence relationship
//     q[k+1] = 2R cos(w_d dt) q[k] - R^2 q[k-1] + b F,   R = e^(-dt/tau)
// The time step is fixed at construction because the coefficients depend on it.
class PlayString {
public:
	PlayString(const String& s_, const Air& air, float w_max, float dt_, float x_listen);

	std::size_t modes() const { return q.size(); }
	float frequency(std::size_t i) const { return w[i]; }   // [rad/s]
	float decay_time(std::size_t i) const { return tau[i]; } // [s]
	float timestep() const { return dt; }
	float get_length() const { return s.get_length(); }

	void set_contact(float x);
	float contact() const { return x_contact; }

	// Excitations add to whatever the string is already doing (superposition).
	void pluck(float x0, float h);   // released from rest in a triangle of height h [m]
	void strike(float x0, float J);  // impulse J [N s]
	void step(float F = 0.0f);       // advance one sample with force F [N] at the contact point

	float output() const { return static_cast<float>(sum(phi_listen)); }  // displacement at x_listen [m]
	double contact_displacement() const { return sum(phi_contact); }

private:
	String s;
	float dt;
	float x_contact = 0.0f;           // 0 = no contact point yet (all shapes vanish)
	std::vector<float> w;             // mode frequency [rad/s]
	std::vector<float> tau;           // decay time [s]
	std::vector<double> q;            // modal displacement now [m]
	std::vector<double> q_prev;       // modal displacement one sample ago [m]
	std::vector<double> c1, c2, b;    // recurrence coefficients (double: c1 ~ 2, c2 ~ 1)
	std::vector<double> back, kick;   // turn pluck / strike into a value for q_prev
	std::vector<float> phi_contact;   // mode shapes at the contact point
	std::vector<float> phi_listen;    // mode shapes at the listening point

	double sum(const std::vector<float>& phi) const {
		double y = 0.0;
		for (std::size_t i = 0; i < q.size(); ++i) y += q[i] * phi[i];
		return y;
	}
};

// ---- Exciters ---------------------------------------------------------------------

// A hammer with a compressible felt tip. Force exists only while the felt is squeezed:
// F = K * compression^p if compression > 0, else 0.
class Hammer {
public:
	Hammer(float M_, float K_, float p_, float v0_, float x0_)
		: M(M_), K(K_), p(p_), x0(x0_), v(v0_) {}
	float step(PlayString& s);  // returns the contact force [N]
private:
	float M;        // hammer mass [kg]
	float K;        // felt stiffness [N/m^p]
	float p;        // felt stiffness exponent
	float x0;       // strike point [m]
	float v;        // hammer velocity [m/s], positive toward the string
	float y = 0.0f; // hammer position [m], starts touching the string at rest
};

// Free vibration for N samples.
inline std::vector<float> render(PlayString& str, std::size_t N) {
	std::vector<float> out(N);
	for (std::size_t k = 0; k < N; ++k) { str.step(); out[k] = str.output(); }
	return out;
}

} 