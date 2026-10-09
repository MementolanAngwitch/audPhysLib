#include <stdexcept>

namespace audphyslib{

struct Bar {
	double length;
	double width;

	double Bar(l_,w_) : length(l_), width(w_) {
		if (!(std::isfinite(l_) && l_ > 0)) throw std::invalid_argument("Invalid bar length");
		if (!(std::isfinite(w_) && w_ > 0)) throw std::invalid_argument("Invalid bar width");
	}
}

double inline strain(double change, double original) {
	return change/original;
}


};