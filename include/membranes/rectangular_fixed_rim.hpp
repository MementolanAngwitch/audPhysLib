#pragma once

#include "membranes.hpp"
#include "2dwave.hpp"

#include <numbers>
#include <complex>

namespace audphys::membrane{
struct RectangularMembrane {
	Membrane membrane;
	Shape shape = Shape::Rectangle;
	double Lx;
	double Lz;
};

// standing wave eq y(x,z,t) = A sin(k_x x)sin(k_z z) e^jwt
inline std::complex<double> rectangular_displacement(const RectangularMembrane& membrane, double w, double x, double z, double t, std::size_t n, std::size_t m) {
	const double pi = std::numbers::pi;
	const double kx = n * pi / membrane.Lx;
	const double kz = m * pi / membrane.Lz;
	const double A = 1.0; // TODO: Max amplitude is determined by initial excitaiton

	return A * std::sin(kx * x) * sin(kz * z) * std::polar(1.0, w*t);
}
}