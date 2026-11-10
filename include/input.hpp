// The resonator of the instrument is simulated leaving the excitation and physical paramters to user input, there are four types of 
// excitation for an instrument: A strike, a pluck, a bow, and a breath. 
// These exictation can be summarized (for the most part) by 4 quantites:
// 1. The force
// 2. The location
// 3. The duration
// 4. The damping
#pragma once
#include <concepts>

class Force{
public: 
	Force(double F_, double t_) : F(F_), t(t_) {};
	double get_F() const {return F;};
	double get_t() const {return t;};
private:
	double F; 		 // Force [N]
	double t;	 	 // Time [s]
};

class Location{
public: 
	Location(double x_, double z_) : x(x_), z(z_) {};
	double get_x() const {return x;};
	double get_z() const {return z;};
private:
	double x; // [m]
	double z; // [m]
};

class Input{
public:
	Input(Force f_, Location l_) : f(f_), l(l_) {};
private:
	Force f;
	Location l;
};

template<class E>
concept Exciter = requires(E e, const PlayString& s, float dt) {
	{e.force}
}