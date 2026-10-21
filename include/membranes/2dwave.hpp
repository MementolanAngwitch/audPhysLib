#pragma once

#include <cmath>
#include "membranes.hpp"

namespace audphys::membrane{

inline double stretched_phase_speed(const Membrane& m){
	return std::sqrt(m.tension / m.ps);
}

inline double stretched_wave_number(const Membrane& m, double w) {
	return w / stretched_phase_speed(m);
}
}