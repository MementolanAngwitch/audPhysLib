#include "strings.hpp"

namespace audphys::strings {

double HammerString::contact_displacement(double t) const {
	const double tau = M * s.wavespeed() / (2.0 * s.get_tension());
	return v0 * tau * (1.0 - std::exp(-t / tau));
}

double HammerString::displacement(double t, double x) const {
	x = s.check_length(x);
	return contact_displacement(t - ((x - x0) / s.wavespeed()));
}

} 