#pragma once 

namespace audphys::membrane {

enum class Shape{ //TODO: How to allow arbitrary shape definitions? Bit mapping?
	Circle, Rectangle
};

struct Membrane{
	double ps; // Surface density
	double tension;  // Tension per unit length
	double thickness; 
};



}