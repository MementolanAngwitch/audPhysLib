#pragma once 

namespace audphys{
// int[a,b] f(x) dx = h/3[ f(x0) + 4*sum{j=1,n} f(x_2j-1) + 2*sum{j=1, n-1} f(x_2j) +f(x2_n)]
// for h = (b-a)/2n and xj = a+jh
template<class F>
inline double simpson(F f, double a, double b, int n) {
	double h = (b-a)/(2.0*n);
	double x0 = a, xn = x0 + 2.0*h;
	double oddSum = 0.0, evenSum = 0.0;
	for (int j = 0; j < n; ++j) {
		oddSum += f(x0 + (2.0*j - 1.0) * h);
		if (j < n-1) evenSum += f(x0 + (2.0*j) * h);
	}
	return h/3.0 * ( f(x0)+ 4.0 * oddSum + 2.0 * evenSum + f(xn) );
}
}