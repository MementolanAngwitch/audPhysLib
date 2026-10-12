#include <cmath>
#include <cstdio>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <iostream>

namespace audphys{
	/*
	Method begins with: stopping tolerance TOL > 0 and points a and b such that f(a)f(b) < 0
	At each iteration the method maintains a,b,c such that b \neq c and
	1. f(b)f(c) < 0 so that a solution lies between b and c if f is continuous
	2. |f(b)| <= |f(c)| so that b can be regarded as the current approximate solution
	3. either a is distinct from b and c, or a = c and is the immediate past value of b

	Each iteration processeds as follow:
	1. if |b-c| <= TOL then the method returns b as the approximate solution
	2. Otherwise, the method determines a trial point b* as follows:
		i. If a = c, then b is determined by linear secant interpolation b* = (af(b)-bf(a))/(f(b)-f(a))
		ii. Otherwise a,b, and c are distinct, and b* is determined using inverse quadratic interpolation:
			- Determine alpha, beta, gamma such that p(y) = alpha y^2 + beta y + gamma satisfies p(f(a)) = a,
			  p(f(b)) = b, p(f(c)) = c
			- Set b* = gamma
	3. If necessary, b* is adjusted or replaced with bisection point
	4. Once b* is finalized, a,b,c and b* are used to determiend the new values of a,b,c

	Rules for step 3:
	Write m = 1/2(c-b) for the step from b to the middle of the bracket. Let delta = 2 epsilon |b| + 1/2 TOL
	be the smallest step worth taking where epsilon \approx 2.2 * 10^-16 is the gap between neighboring double values near 1
	1. Only interpolate if the last step helped. Interpolation is attempted only if |f(a)| > |f(b)| meaning the previous step
	   actually reduced |f|, and the step taken two iterations ago was at least delta. Otherwise b* = b + m, a bisection step
	2. The trial point must stay in the near part of the bracket. b* is accepted only if it lies between b and the point three 
	   quarters of the way from b to c:
	   									|b* - b| < 3/4 |c-b|
	   Interpolation can throw a point outside the bracket, or almost onto its fare end, which wastes the step
	3. The steps must keep shrinking. b* is also accepted only if 
										|b* - b| < 1/2 |e|
	   where e is the step taken two iterations earlier. This forces the step size to at least halve every two iterations, which
	   is as fast bisection. If interpolation is slow, this rule hands over to bisection
	4. Adjustment. If the accepted step is smaller than delta, its lengthened to exactly delta, toward c:
										b* = b + delta * sign(m)
	
	Rules for step 4
	1. The old best becomes the previous point, the trial point becomes the new best: a <- b then b <- b*
	2. If f(b) now has the same sign as f(c), the root is no longer between them. But the old b, now stored in a, had the 
	   opposite sign to c, so its on the other side of the root. Set c <- a, the step history is also reset here, e <- b-a 
	   because the bracket has changed shape.
	3. If |f(c)| < |f(b)| swap b and c and set a to the point that was b before the swap. After this swap, a=c, so the next iteration
	   has only two distinct points and uses the secant. 
	NOTE: 3/4 and 1/2 are choices not deep constants but ensures speed.
	*/

template<class F>
double Brent(F f, double x0, double x1, double TOL = 1e-15, bool log = false) {
	const double eps = std::numeric_limits<double>::epsilon();
	double prev = x0, best = x1, other = 0.0;

	double fprev = f(x0), fbest = f(x1), fother = 0.0;
	if (fprev*fbest > 0) throw std::invalid_argument("f(x0)f(x1) exceeds 0");

	if (std::abs(fprev) < std::abs(fbest)) {
		other = best;
		fother = fbest;
		best = prev;
		fbest = fprev;
	} 
	else{
		other = prev; fother = fprev;
	}
	while (true) {
		if (std::abs(best - other) < TOL) return best;
		if (other == prev) { //b* = (af(b)-bf(a))/(f(b)-f(a))
			double canidate = (prev * fbest - best * fprev) / (fbest - fprev);
		}
		else {
			// Determine alpha, beta, gamma such that p(y) = alpha y^2 + beta y + gamma 
			// satisfies p(f(a)) = a, p(f(b)) = b, p(f(c)) = c
			// Set b* = gamma
			// p(y) = a (y-fb)(y-fc)/(fa-fb)(fa-fc) + b (y-fa)(y-fc)/(fb-fa)(fb-fc) + c (y-fa)(y-fb)/(fc-fa)(fc-fb)
			// set y = 0 and solve for b
			// b* = a fbfc/(fa-fb)(fa-fc) + b fafc/(fb-fa)(fb-fc) + c fafb/(fc-fa)(fc-fb) 

			double candiate = prev * (fbest * fother) / ((fprev - fbest)*(fprev - fother))
							  + best * (fprev * fother) / ((fbest - fprev)*(fbest - fother)) 
							  + other * (fprev * fbest) / ((fother - fprev)*(fother - fbest));

		}
		
	}
}

}; //namespace audphyslib