#pragma once

#include <numbers>
#include <vector>
#include <algorithm>
#include <cmath>
#include "2dwave.hpp"
#include "membrane_modes.hpp"

namespace audphys::membrane {

struct RectangularMode {
	Membrane membrane;
	std::size_t n,m;
	double Lx, Lz; // side lengths
	double kx, kz; // n pi/Lx, m pi/Lz
	double w;
	double A = 1.0, phi = 0.0, beta = 0.0; //picking phi =0.0 so y(x,z,0) = 0
};

inline RectangularMode rectangular_mode(const Membrane& mem, double Lx, double Lz, std::size_t n, std::size_t m) {
	const double pi = std::numbers::pi; // w = kc = sqrt(kx^2 + kz^2) c
	const double kx =  n*pi/Lx;
	const double kz = m*pi/Lz;
	return { mem,n,m,Lx,Lz,kx,kz, stretched_phase_speed(mem) * std::sqrt(kx*kx + kz*kz)};
}

inline std::vector<RectangularMode> rectangular_modes(const Membrane& mem, double Lx, double Lz, std::size_t n_max, std::size_t m_max) {
	std::vector<RectangularMode> out;
	out.reserve(n_max * m_max);
	for (std::size_t n = 1; n <= n_max; ++n) {
		for (std::size_t m = 1; m <= m_max; ++m) {
			out.push_back(rectangular_mode(mem, Lx, Lz, n, m));
		}
	}
	std::sort(out.begin(), out.end(), [](const RectangularMode& p, const RectangularMode& q) { return p.w < q.w; });
	return out;
}

inline double rectangular_mode_shape(const RectangularMode& md, double x, double z) {
	return std::sin(md.kx * x) * std::sin(md.kz * z);
}

inline std::vector<double> nodal_lines_x(const RectangularMode& md) {
	std::vector<double> out;
	for (int q = 1; q < md.n; ++q) out.push_back(q * md.Lx / md.n);
	return out;
}
 
inline std::vector<double> nodal_lines_z(const RectangularMode& md) {
	std::vector<double> out;
	for (int p = 1; p < md.m; ++p) out.push_back(p * md.Lz / md.m);
	return out;
}

// For rectangular fixed rim, the real standing wave associated with the n,m mode can be written as 
// ynm(x,z,t) = Anm sin(k_xn x) sin(K_zm z) sin(w_nm t + phi_nm)
// Anm = 4V/(LxLz) (sin(k_xn x_0) sin(k_zm z_0))/w_nm, V = I/ps

inline void amplitude(std::vector<RectangularMode>& mds, double I, double x0, double z0) {
	//Check that x0,z0 is on the membrane
	for (RectangularMode& md : mds) {
		const double V = I / md.membrane.ps;
		md.A = 4.0 * V / (md.Lx * md.Lz) * (std::sin(md.kx * x0) * std::sin(md.kz * z0)) / md.w;
	}
}

inline double mode_displacement(const RectangularMode& md, double x, double z, double t) {
	return md.A * std::sin(md.kx * x)*std::sin(md.kz * z)*sin(md.w * t + md.phi);
}

inline double total_mode_displacement(const std::vector<RectangularMode>& mds, double x, double z, double t) {
	double sum = 0.0;
	for (const RectangularMode& m : mds) {
		sum += mode_displacement(m, x,z,t);
	}
	return sum;
}

}