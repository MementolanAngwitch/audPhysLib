#pragma once

template<class F>
double brent(F f, double x0, double x1, double TOL = 1e-15, std::size_t N = 1000);
template <class F>
std::vector<std::pair<double,double>> bracket_finder(F f, double lo, double hi, int steps);