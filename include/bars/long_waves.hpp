#include <stdexcept>
#include <numbers>
#include <complex>

namespace audphyslib{

struct Bar {
	double length;
	double width;
	double ygMod;
	double density;

	double Bar(l_,w_, y_) : length(l_), width(w_), ygMod(y_) {
		if (!(std::isfinite(l_) && l_ > 0)) throw std::invalid_argument("Invalid bar length");
		if (!(std::isfinite(w_) && w_ > 0)) throw std::invalid_argument("Invalid bar width");
		if (!(std::isfinite(w_) && y_ > 0)) throw std::invalid_argument("Invalid young's modulus");
		if (!(std::isfinite(d_) && d_ > 0)) throw std::invalid_argument("Invalid density")
	}
}
// using y to represent longitudinal displacement
double inline strain(const Bar& b,double change) {
	return change/b.length;
}
// hooke's law
double inline stress(const Bar& b, double strain) {
	return -1.0 * b.ygMod * strain;
}


// eq of mtion d^2 y/ d^2 x = 1/c^2 d^2 y/dt^2

double inline phase_speed(const Bar& b) {
	return std::sqrt(b.ygMod / b.density);
}

double inline wave_number(const Bar& b, double w) {
	return w / phase_speed(b);
}

// y(x,t) = A exp(j(wt-kx)) + B exp(j(wt+kx))

struct WaveAmplitude {
	std::complex<double> A;
	std::complex<double> B;
};

struct ModeWaveAmplitudes {
	std::vector<WaveAmplitude> amps;
}

double inline position(const WaveAmplitude& amps, const Bar& b, double w, double x, double t) {
	const double k = wave_number(b,w);
	return std::polar(A,w*t-k*x) + std::polar(B,w*t+k*x)
}

// Rigidly fixed at both ends

double inline mode_wave_number(const Bar& b, std::size_t n) {
	return (n * std::numbers::pi / b.length);
}

double inline mode_ang_freq(const Bar& b, std::size_t n) {
	return (n * std::numbers::pi * phase_speed(b)) / b.length;
}

std::complex<double> inline complex_displacement(const Bar& b, const ModeWaveAmplitudes& AmpList, std::size_t n, double w, double x, double t) {
	return -2.0 * std::complex(0, std::polar(AmpList.amps[n].A, mode_ang_freq(b,n)) * std::sin(mode_wave_number(b,n) * x));
}

double inline real_displacement(const Bar& b, const WaveAmplitudes)
	
};