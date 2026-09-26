#pragma once 
#include <cmath>
#include <numbers>
#include <cassert>
namespace audphys {

// omega = 2pi * f
inline double natural_angular_freq(double s, double m) {
	assert(s > 0 && m > 0); //negative stiffness and negative mass are nonsensical
	return std::sqrt(s / m);
}

inline double natural_angular_freq(double nat_freq) {
	assert(nat_freq > 0);
	return nat_freq * 2 * std::numbers::pi;
}

inline double natural_frequency(double s, double m) {
	assert(s > 0 && m > 0);
	return std::sqrt(s/m) / (2 * std::numbers::pi);
}

inline double natural_frequency(double nat_ang_freq) {
	assert(nat_ang_freq > 0);
	return nat_ang_freq / (2 * std::numbers::pi);
}
// T = 1/f
inline double period(double f0) {
	assert(f0 > 0);
	return 1.0/f0;
}

// Wave equation: x = A1 cos (omega * t) + A2 cos (omega * t)

struct InitialConditions { double omega0; double x0; double u0; };
struct AmplitudePhase    { double omega0; double A;  double phi; };
// A = (x^2 + (u/omega)^2)^1/2, phi = tan^-1 (-u/(x * omega))
inline AmplitudePhase to_amp_phase(InitialConditions ic) {
	double A2 = ic.u0 / ic.omega0;
	return {
		ic.omega0,
		std::sqrt(ic.x0 * ic.x0 + A2 * A2),
		atan2(-A2, ic.x0)
	};
}
// x = x0 cos(omega * t) + (u/omega) sin (omega * t)
inline double position(InitialConditions ic, double t) {
	assert (ic.omega0 > 0);
	return ic.x0 * std::cos(ic.omega0 * t) + (ic.u0 / ic.omega0) * std::sin(ic.omega0 * t);
}
// x = A cos (omega * t + phi)
inline double position(AmplitudePhase ap, double t) {
	return ap.A * std::cos(ap.omega0 * t + ap.phi);
}
// U = omega * A
inline double speedAmplitude(AmplitudePhase ap) {
	return ap.omega0 * ap.A;
}
// u = -U sin(omega0 * t + phi)
inline double speed(AmplitudePhase ap, double t) {
	return -1.0 * speedAmplitude(ap) * std::sin(ap.omega0 * t + ap.phi);
}
//a = -omega * U * cos(omega * t + phi)
inline double acceleration(AmplitudePhase ap, double t) {
	return -1.0 * ap.omega0 * speedAmplitude(ap) * cos(ap.omega0 * t + ap.phi);
} 
//Energy

}