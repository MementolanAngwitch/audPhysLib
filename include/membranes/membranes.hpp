#pragma once 

namespace audphys {

enum class Shape{ //TODO: How to allow arbitrary shape definitions? Bit mapping?
	Circle, Rectangle
};

struct Membrane{
	Shape shape
	double ps; // Surface density
	double T;  // Tension per unit length
};



}