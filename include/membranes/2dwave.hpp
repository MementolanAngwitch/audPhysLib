#pragma once

#include <cmath>
#include "membranes.hpp"

namespace audphys::membrane{

inline double stretched_phase_speed(const Membrane& m){
	return std::sqrt(m.T / m.ps);
}
}