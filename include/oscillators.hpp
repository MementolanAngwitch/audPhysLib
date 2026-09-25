#pragma once 
#include <cmath>
#include <numbers>
#include <cassert>
namespace audphys {

inline double natural_angular_freq(double s, double m) {
	assert(s > 0 && m > 0); //negative stiffness and negative mass are nonsensical
	return std::sqrt(s / m);
}

inline double natural_frequency(double nat_ang_freq) {
	assert(nat_ang_freq > 0);
	return nat_ang_freq / (2 * std::numbers::pi);
}

inline double period(double f0) {
	assert(f0 > 0);
	return 1.0/f0;
}



}