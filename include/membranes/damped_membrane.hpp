#pragma once 

#include <complex>
#include <cmath>
#include "2dwaves.hpp"
#include "membranes.hpp"
#include "rectangular_fixed_rim.hpp"
#include "circular_fixed_rim.hpp"

namespace audphys::membrane{

inline double mode_shape(const MembraneMode& mode) {

}
//modeShape * e^(-beta*t) * e^(j(wd)t)
inline double mode_displacement(const MembraneMode& mode, double t) {
	
}

}