#pragma once

#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>

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
	   quarters of the way from b to c also b* - b has same sign as c-b,
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
double brent(F f, double x0, double x1, double TOL = 1e-15) {
	const double eps = std::numeric_limits<double>::epsilon();
	
	double fx0 = f(x0), fx1 = f(x1);
	if (fx0 == 0.0) return x0;
	if (fx1 == 0.0) return x1;
	if (fx0*fx1 > 0) throw std::invalid_argument("Brent Method: fx0 and fx1 must have opposite signs");
	double best = x1, fbest = fx1, other = x0, fother = fx0;
	if (std::abs(fx0) < std::abs(fx1)) { best = x0, fbest = fx0, other = x1, fother = fx1;}
	
	double prev = other, fprev = fother;


	//Iteration checking
	double d = other -best; //last step
	double e = d; //step before last
	for (int iter = 0; iter < 200; ++iter) {
		const double step = 0.5 * (other - best);
		const double delta = 2.0 * eps * std::abs(best) + 0.5 * TOL;
		//step 1
		if (fbest == 0.0 || std::abs(best - other) <= TOL + 4.0 * eps * std::abs(best)) return best;
		double candidate;
		if (std::abs(e) >= delta && std::abs(fprev) > std::abs(fbest)) {	// rule 1: interpolate only if the last step helped
			if (other == prev) {		// secant: b* = (af(b)-bf(a))/(f(b)-f(a))
				candidate = (prev * fbest - best * fprev) / (fbest - fprev);
			} else {					// inverse quadratic: b* = gamma
			// Determine alpha, beta, gamma such that p(y) = alpha y^2 + beta y + gamma 
			// satisfies p(f(a)) = a, p(f(b)) = b, p(f(c)) = c
			// Set b* = gamma
			// p(y) = a (y-fb)(y-fc)/(fa-fb)(fa-fc) + b (y-fa)(y-fc)/(fb-fa)(fb-fc) + c (y-fa)(y-fb)/(fc-fa)(fc-fb)
			// set y = 0 and solve for b
			// b* = a fbfc/(fa-fb)(fa-fc) + b fafc/(fb-fa)(fb-fc) + c fafb/(fc-fa)(fc-fb) 
				candidate = prev  * (fbest * fother) / ((fprev  - fbest) * (fprev  - fother))
				          + best  * (fprev * fother) / ((fbest  - fprev) * (fbest  - fother))
				          + other * (fprev * fbest)  / ((fother - fprev) * (fother - fbest));
			}
			const double trial = candidate - best;
			const bool toward_other = (trial > 0) == (step > 0);						// rule 2: direction
			if (toward_other && std::abs(trial) < 0.75 * std::abs(other - best)		// rule 2: near part
			                 && std::abs(trial) < 0.5 * std::abs(e)) {					// rule 3: shrinking
				e = d;  d = trial;			// accepted: shift the history
			} else {
				candidate = best + step;	// rejected: bisect
				d = e = step;
			}
		} else {
			candidate = best + step;		// rule 1 failed: bisect
			d = e = step;
		}
		if (std::abs(candidate - best) < delta) {			// rule 4: at least delta, toward c
			candidate = best + std::copysign(delta, step);
		}

		//step 4
		prev = best; fprev = fbest;
		best = candidate; fbest = f(candidate);

		if (fbest * fother > 0) {
			other = prev;  fother = fprev;
			d = e = best-prev;
			//IMPLEMENT: step history is also reset here, e <- b-a 
		}

		if (std::abs(fother) < std::abs(fbest)) {
			double temp = best, ftemp = fbest;
			best = other; fbest = fother;
			other = temp; fother = ftemp;
			prev = temp; fprev = ftemp;
		}
	}
	return best;
}

} //namespace audphys