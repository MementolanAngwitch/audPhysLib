#pragma once

namespace audphys {

template <class F>
inline std::vector<std::pair<double,double>> bracket_finder(F f, double lo, double hi, int steps) {
	std::vector<std::pair<double,double>> out;
	const double h = (hi - lo) / steps;
	double a = lo, fa = f(a);
	for (int i = 1; i <= steps; ++i) {
		double b = lo + i * h, fb = f(b);
		if (fa * fb < 0) {
			out.emplace_back(a,b);
		}
		a = b;
		fa = fb;
	}
	return out;
}

}