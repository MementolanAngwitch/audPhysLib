#pragma once

#include <stdexcept>
#include <cmath>
#include <complex>
#include <vector>

namespace audphys::string{

struct String {
	double length;		// L     [m]
	double tension;		// T     [N]
	double lin_density;	// rho_L [kg/m], uniform
 
	String(double L_, double T_, double rhoL_) : length(L_), tension(T_), lin_density(rhoL_) {
		if (!(std::isfinite(length) && length > 0))           throw std::invalid_argument("String: length must be positive and finite");
		if (!(std::isfinite(tension) && tension > 0))         throw std::invalid_argument("String: tension must be positive and finite");
		if (!(std::isfinite(lin_density) && lin_density > 0)) throw std::invalid_argument("String: linear density must be positive and finite");
	}
};
// m = L * pl
inline double string_mass(const String&s) {
	return s.lin_density * s.length;
}

enum class End { Fixed, Free };

struct WaveAmplitudes {
	std::complex<double> A;		// wave travelling away from the driver   [m]
	std::complex<double> B;		// wave reflected back from the load      [m]
};

struct Mode {
	int n;			// mode number, 1, 2, ...
	double A;		// cosine amplitude   [m]
	double B;		// sine amplitude     [m]
	double beta;	// decay rate         [1/s], 0 = undamped
};
 
struct NormalModes {
	String s;					// the string these modes belong to
	std::vector<Mode> modes;
};

}