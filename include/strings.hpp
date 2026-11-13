#pragma once
#include <cmath>
#include <cstddef>
#include <numbers>
#include <stdexcept>
#include <vector>
#include <concepts>
#include <complex>
#include "tools.hpp"
#include "input.hpp"

// TODO: Route class Input into strings, final syntax should like something like 
// gen(String, Air, Input) -> std::vector<float> with size t/dt. 
namespace audphys::strings {

enum class BoundaryCondition { Fixed, Free }; // TODO: Replace BoundaryCondition with Support since real supports are not perfectly fixed

class Support { //cotkL = m/M kL where M = \mu * L is the total mass of the string, solve for k and feed
				//kn = npi/L 
	float m; 
};

class String {
public:
	String(float L_, float mu_, float T_, float r_, std::complex<float> E_,
	       BoundaryCondition left_ = BoundaryCondition::Fixed,
	       BoundaryCondition right_ = BoundaryCondition::Fixed)
		: L(L_), mu(mu_), T(T_), r(r_), E(E_), left(left_), right(right_) {}
	float get_length() const { return L; };
	float get_mu() const { return mu; };
	float get_rho() const { return mu / (std::numbers::pi * r * r); };
	float get_tension() const { return T; };
	float wavespeed() const { return std::sqrt(T / mu); };
	float get_radius() const { return r;};
	float fundamental() const { return std::numbers::pi_v<float> * wavespeed() / L; } //[rad/s]
	std::complex<float> get_young_mod() const { return E;};
	std::complex<float> impedance() const { return std::sqrt(T * mu); };
	BoundaryCondition get_left() const { return left; };
	BoundaryCondition get_right() const { return right; };
	// Folds a position back onto the string and returns it.
	float check_length(float x) const {
		if (std::abs(x) > L) x = L - (x - L);
		if (x < 0) x = -x;
		return x;
	}
private:
	float L;        // length [m]
	float mu;       // linear density [kg/m]
	float T;        // tension [N]
	float r;        // radius [m]
	std::complex<float> E;   // Young's modulus [N/m^2]
	std::size_t parts = 0;   // # of divisions of the string
	BoundaryCondition left;  // negative direction
	BoundaryCondition right; // positive direction
};

class Air {
public:
	Air() {};
	Air(float rho_, float eta_) : rho(rho_), eta(eta_) {};
	float get_rho() const { return rho; };
	float get_eta() const { return eta; };
private:
	float rho = 1.2;    // air density [kg/m^3]
	float eta = 1.5e-5; // kinematic viscosity [m^2/s]
};


// Damping decay times: A = A_0 * e^(-t/tau)
float tau_air(const Air& air, const String& s, float w);

float tau_internal(const String& s, float w);

// Support (bridge) losses. With a frequency-independent conductance G at both ends,
// every mode loses the same fraction of energy per round trip, so tau is the same for all modes:
// tau = 1 / (8 mu L f1^2 G) = L / (2 T G). Frequency dependence enters only through G(f).
float tau_support(const String& s, float G);

float tau_total(const Air& air, const String& s, float w, float G);

// sin(k+1 theta) = 2cos(theta) * sin(k theta) - sin((k-1)theta) where theta = pi x /L
void fill_shapes(float x, float L, std::vector<float>& out);

// q_n'' + (2/tau_n) q_n' + w_n^2 q_n = (2 / (mu L)) * phi_n(x0) * F(t),
// with mode shape phi_n(x) = sin(n pi x / L).
// On a fixed-fixed string
class PlayString {
public: 
	PlayString(const String& s_, const Air& air_, float G, float f_max, float dt, float x_listen)
			   : s(s_), air(air_) {
		const float L = s.get_length();
		const std::size_t n = static_cast<std::size_t>(f_max / (s.fundamental() / (2.0f * std::numbers::pi_v<float>)));
		const float g = 2.0f / (s.get_mu() * L);
		tau.assign(n, 0.0f); w.assign(n, 0.0f); q.assign(n, 0.0f); q_prev.assign(n, 0.0f);
		c1.assign(n, 0.0f); c2.assign(n, 0.0f); b.assign(n, 0.0f); back.assign(n, 0.0f); kick.assign(n, 0.0f);
		phi_contact.assign(n, 0.0f); phi_listen.assign(n, 0.0f);
		for (std::size_t i = 0; i < n; ++i) {
			w[i]   = (i+1) * std::numbers::pi * s.wavespeed() / L;
			tau[i] = tau_total(air_, s_, w[i], G);
			const double a  = 1.0 / tau[i];                    // decay rate alpha = 1/tau
			const double wd = std::sqrt(w[i]*w[i] - a*a);      // damped frequency
			const double R  = std::exp(-a * dt);
			const double th = wd * dt;
			c1[i]   = 2.0 * R * std::cos(th);
			c2[i]   = R * R;
			b[i]    = g * dt * R * std::sin(th) / wd;          // times phi_contact[i] * F
			back[i] = std::exp(a * dt) * (std::cos(th) - a * std::sin(th) / wd); // q[-1]/q[0], released at rest
			kick[i] = std::exp(a * dt) * std::sin(th) / wd;    // q[-1] per unit velocity kick
		}
		fill_shapes(x_listen, L, phi_listen);
	}
	std::size_t modes() const { return q.size(); } 

	void set_contact(double x) {x_contact = x; fill_shapes(x, s.get_length(), phi_contact);};
	void pluck(float x0, float h);
	void strike(float x0, float J);
	void step(float dt, float F =0.0);
	float output() const { return sum(phi_listen); };

	float get_length() const { return s.get_length(); };
	double contact_displacement() const { return sum(phi_contact); };

private:
	String s;
	Air air;
	std::vector<float> tau;    // decay time [s]
	std::vector<double> q;    // modal displacement now [m]
	std::vector<double> q_prev; // modal displacement one sample ago [m]
	std::vector<double> c1, c2, b;  // recurrence coefficients, fixed per mode
	std::vector<float> back, kick; // for turning pluck/strike into two samples
	std::vector<float> w;	 // angular frequency [rad/s]
	std::vector<float> phi_contact;
	std::vector<float> phi_listen; // modes at contact and at listen
	double x_contact = 0.0;

	double sum(const std::vector<float>& phi) const {
		double y = 0.0;
		for (std::size_t i = 0; i < q.size(); ++i) y += q[i] * phi[i];
		return y;
	}
};
// A hammer with a compressible felt tip. Force exists only while the felt is squeezed:
// F = K * compression^p  if compression > 0, else 0.
class Hammer {
public:
	Hammer(float M_, float K_, float p_, float v0_) 
	 : M(M_), K(K_), p(p_), v(v0_) {};
	float step(PlayString& s, float dt);
private:
	float M; 	    // hammer mass [kg]
	float K;       // felt stiffness [n/m^p]
	float p; 	    // felt stiffness exponent 
	float x0; 	    // strike point [m]
	float v;       // hammer velocity [m/s], postive toward string
	float y = 0.0; // hammer position [m], starts touching the string at rest
};

inline std::vector<float> render(PlayString& str, std::size_t N, float dt){
	std::vector<float> out(N);
	for (std::size_t k = 0; k < N; ++k) { str.step(dt); out[k] = str.output(); }
	return out;
}
}