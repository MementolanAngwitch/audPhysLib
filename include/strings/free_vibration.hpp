#pragma once

#include <cmath>
#include "1d_wave.hpp"

namespace audphys::strings{

// The shape y0 on [0, L], extended to the whole line by mirror images across each end.
// Moving 2L along the line multiplies the shape by (left sign)(right sign).

template<class F>
inline double extended_shape(F y0, double L, End left, End right, double xi) {
	const double bc = reflection_sign(right);
	const double flip_per_cycle = reflection_sign(left) * bc;

	const double k = std::floor(xi / (2*L)); // which 2L cycle
	const double v = xi -(k * 2 * L);         // position within the cycle 

	const bool odd_cycle = std::fmod(std::abs(k), 2.0) == 1.0;
	const double cycle_sign = (flip_per_cycle < 0 && odd_cycle) ? -1.0 : 1.0;

	const double y = (v <= L) ? y0(v) : bc* y0(2*L - v);
	return y * cycle_sign;
}

// String released AT REST from shape y0 (a pluck), ideal ends:
// y(x,t) = 1/2 [ Y0(x - ct) + Y0(x + ct) ] 

template<class F>
inline double position(F y0, const String& s, End left, End right, double x, double t) {
	const double L = s.length; 
	const double c = wave_speed(s);
	return 0.5 * (extended_shape(y0, L, left, right, x - c*t) + extended_shape(y0,L,left,right,x+c*t));
}

}
 