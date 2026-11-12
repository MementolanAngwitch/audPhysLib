#include "analytic.hpp"


float HammerString::contact_displacement(float t) const {
	const float tau = M * s.wavespeed() / (2.0 * s.get_tension());
	return v0 * tau * (1.0 - std::exp(-t / tau));
}

float HammerString::displacement(float x, float t) const {
	if(t < 0) return 0.0f;
	x = s.check_length(x);
	return contact_displacement(t - (std::abs(x - x0) / s.wavespeed()));
}