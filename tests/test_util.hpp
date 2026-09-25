#pragma once 
#include <cmath>
#include <algorithm>

inline bool close(double a, double b, double tol = 1e-9) { return std::abs(a - b) <= tol * std::max(1.0, std::abs(b)); }