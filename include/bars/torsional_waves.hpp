#pragma once 
#include <numbers>
#include "long_waves.hpp"

namespace audphys{
// For circular cross section
// stress = shear_mod * radius * (stress_angle / plate_length)
// df = shear_mod (plate_width * plate_thickness) * radius * (stress_angle / plate_length)
// dt = shear_mod (2pi * r^3) (stress_angle / plate_length) * tube_thickness
// t = shear_mod 1/2 pi a^4 (stress_angle / plate_length)
// c^2 = shear_mod/density
// torsional bc: (1) fixed end, phi = 0.
// (2) free end, t = 0 dphi/dx = 0 
// (3) mass-loaded, dphi/dx = (dphi^2/dx^2)I where I is inertia of the load about the axis of the bar.
// solutions for longitudinal waves and waves on strings apply here

inline double stress(const Bar& b, double stress_angle, double plate_length) {
	const double r = std::sqrt(b.crossSection.area / std::numbers::pi);
	return b.shearMod * r * stress_angle / plate_length;
}

inline double force(const Bar& b, double plate_width, double plate_length, double plate_thickness, double stress_angle) {
	const double r = std::sqrt(b.crossSection.area / std::numbers::pi);
	return b.shearMod * (plate_width * plate_thickness) * r * (stress_angle / plate_length);
}

inline double torsion(const Bar& b, double stress_angle, double plate_length, double tube_thickness) {
	const double r = std::sqrt(b.crossSection.area / std::numbers::pi);
	return b.shearMod * (2.0 * std::numbers::pi * r * r *r) * (stress_angle / plate_length) * tube_thickness;
}

inline double total_torsion(const Bar& b, double stress_angle, double plate_length, double a) {
	const double r = std::sqrt(b.crossSection.area / std::numbers::pi);
	if (a > r) a = r;
	return b.shearMod * 1.0 / 2.0 * std::numbers::pi * a * a * a * a * (stress_angle/plate_length);
}

inline double torsional_phase_speed(const Bar& b) {
	return std::sqrt(b.shearMod / b.density);
}

}