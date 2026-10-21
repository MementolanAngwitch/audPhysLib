#pragma once

#include <cmath>

namespace audphys::membrane {

struct Membrane {
	double ps;		// surface density
	double tension; // tension per unit length of edge
}
}