#include "strings.hpp"

double PluckedString::mode_displacement(std::size_t i, double x, double t) {
	const double w_i = i * std::numbers::pi * s.wavespeed() / s.get_length(); 
	return ( a[i] * std::sin(w_i * t) + b[i] * std::cos(w_i * t) ) * std::sin((w_i*x)/s.wavespeed());
}

double PluckedString::displacement(double x, double t) {
	sum = 0.0;
	for (std::size_t i = 0; i < n; ++i) {
		sum += mode_displacement(i, x, t);
	}
	return sum;
}

void PluckedString::make_modes() {
	for (std::size_t i = 0; i < n; ++i) {
		const double w_i = i * std::numbers::pi * s.wavespeed() / s.get_length(); 
		A.push_back( 2.0 / (w_i * s.L)  * simpson([&](x)->double{ } )
		B.push_back(0.0);
	}
}