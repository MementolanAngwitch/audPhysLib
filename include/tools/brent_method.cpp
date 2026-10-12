#include <cmath>
#include <iostream>

int main() {
	/*
	Method begins with: stopping tolerance delta > 0 and points a and b such that f(a)f(b) < 0
	At each iteration the method maintains a,b,c such that b \neq c and
	1. f(b)f(c) < 0 so that a solution lies between b and c if f is continuous
	2. |f(b)| <= |f(c)| so that b can be regarded as the current approximate solution
	3. either a is distinct from b and c, or a = c and is the immediate past value of b

	Each iteration processeds as follow:
	1. if |b-c| leq sigma
	*/
}