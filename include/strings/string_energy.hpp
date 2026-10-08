#pragma once 

namespace audphys{
	// Energy per unit length: 1/2 rho_L u^2 (kinetic) + 1/2 T (dy/dx)^2 (potential)
	// Energy of one mode: E_n = 1/4 m_s w_n^2 (A_n^2 + B_n^2), modes add (they are orthogonal)

	// Energy per unit length at (x, t) 

	inline double energy_density(const NormalModes& nm, double x, double t) {
		const double u = velocity(nm, x, t), d = slope(nm, x, t);
		return 0.5 * nm.s.lin_density * u * u + 0.5 * nm.s.tension * d * d;
	}

	// E_n at time t  
	inline double mode_energy(const String& s, const Mode& m, double t = 0.0) {
		const double w = mode_angular_frequency(s, m.n);
		return 0.25 * string_mass(s) * w * w * (m.A * m.A + m.B * m.B) * std::exp(-2 * m.beta * t);
	}
	// Total energy at time t

	inline double total_energy(const NormalModes& nm, double t = 0.0) {
		double E = 0.0;
		for (const Mode& m : nm.modes) E += mode_energy(nm.s, m, t);
		return E;
	}
}