#include "analytic.hpp"

namespace audphys::strings {

float HammerString::contact_displacement(float t) const {
	if (t < 0.0f) return 0.0f;  // the hammer hasn't hit yet
	const float tau = M * s.wavespeed() / (2.0f * s.get_tension());
	return v0 * tau * (1.0f - std::exp(-t / tau));
}

float HammerString::displacement(float x, float t) const {
	x = s.check_length(x);
	// The wave reaches x after |x - x0| / c; before that, contact_displacement returns 0.
	return contact_displacement(t - std::abs(x - x0) / s.wavespeed());
}

} // namespace audphys::strings