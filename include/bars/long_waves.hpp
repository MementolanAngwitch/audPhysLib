#include <stdexcept>
#include <numbers>
#include <complex>
#include <vector>

namespace audphyslib{

struct Bar {
	double length;
	double crossSectionArea;
	double youngMod;
	double density;

	Bar(double l_,double a_,double y_,double d_) : length(l_), crossSectionArea(w_), younggMod(y_), density(d_){
		if (!(std::isfinite(l_) && l_ > 0)) throw std::invalid_argument("Invalid bar length");
		if (!(std::isfinite(w_) && w_ > 0)) throw std::invalid_argument("Invalid bar cross sectional area");
		if (!(std::isfinite(w_) && y_ > 0)) throw std::invalid_argument("Invalid young's modulus");
		if (!(std::isfinite(d_) && d_ > 0)) throw std::invalid_argument("Invalid density");
	};
};
// using y to represent longitudinal displacement
inline double strain(const Bar& b,double change) {
	return change/b.length;
}
// hooke's law
inline double stress(const Bar& b, double strain) {
	return -1.0 * b.ygMod * strain;
}


// eq of mtion d^2 y/ d^2 x = 1/c^2 d^2 y/dt^2

inline double phase_speed(const Bar& b) {
	return std::sqrt(b.ygMod / b.density);
}

inline double wave_number(const Bar& b, double w) {
	return w / phase_speed(b);
}

// y(x,t) = A exp(j(wt-kx)) + B exp(j(wt+kx))

struct WaveAmplitude {
	std::complex<double> A;
	std::complex<double> B;
};

struct WaveModeAmplitudes {
	std::vector<WaveAmplitude> amps;
};

inline std::complex<double> position(const WaveAmplitude& amps, const Bar& b, double w, double x, double t) {
	const double k = wave_number(b,w);
	return amps.A * std::polar(1.0,w*t-k*x) + amps.B * std::polar(1.0,w*t+k*x);
}

// Rigidly fixed at both ends

inline double mode_wave_number(const Bar& b, std::size_t n) {
	return (n * std::numbers::pi / b.length);
}

inline double mode_ang_freq(const Bar& b, std::size_t n) {
	return (n * std::numbers::pi * phase_speed(b)) / b.length;
}

inline std::complex<double> fixed_complex_displacement(const Bar& b, const WaveModeAmplitudes& AmpList, std::size_t n, double x, double t) {
	const std::complex<double> j{ 0.0, 1.0 };
	return -2.0 * j * AmpList.amps[n].A
	     * std::polar(1.0, mode_ang_freq(b, n) * t)		// e^{j w_n t}
	     * std::sin(mode_wave_number(b, n) * x);
}

inline double fixed_real_displacement(const Bar& b, const WaveModeAmplitudes& AmpList, std::size_t n, double x, double t) {
	const double wn = mode_ang_freq(b, n);
	return std::real(AmpList.amps[n].A * std::cos(wn * t) + AmpList.amps[n].B * std::sin(wn*t)) * std::sin(mode_wave_number(b,n)*x);
}

inline double total_fixed_real_displacement(const Bar& b, const WaveModeAmplitudes& AmpList, double x, double t) {
	std::complex<double> out = 0.0;
	for (std::size_t i = 0; i < AmpList.amps.size(); ++i) {
		out += fixed_real_displacement(b,AmpList,i,x,t);
	}
	return std::real(out);
}

//free, free bar

inline std::complex<double> free_complex_displacement(const Bar& b, const WaveModeAmplitudes& AmpList, std::size_t n, double x, double t) {
	const std::complex<double> j{0.0, 1.0};
	return 2 * AmpList.amps[n].A 
		   * std::polar(1.0, mode_ang_freq(b,n) * t)
		   * std::cos(wave_number(b,n) * x);
}

inline double free_real_displacement(const Bar& b. const WaveModeAmplitudes& AmpList, std::size_t n, double x, double t) {
	const double wn = mode_ang_freq(b,n);
	return std::real(
		AmpList.amps[n].A * std::cos(wn*t) 
		+ AmpList.amps[n].B * std::sin(wn*t))
		* std::cos(wave_number(b,n) * x);
}

inline double total_free_real_displacement(const Bar& b, const WaveModeAmplitudes& Amplist, double x, double t) {
	std::complex<double> out = 0.0;
	for (std::size_t i = 0; i < AmpList.amps.size(); ++i ) {
		out += free_real_displacement(b, Amplist, i, x,t);
	}
	return std::real(out);
}
// free-fixed bar used free-free displacment eqs.

inline double free_fixed_mode_wave_number(const Bar& b, std::size_t n) {
	return ((2*n-1) * std::numbers::pi)/(2 * b.length);
}

inline double free_fixed_mode_ang_freq(const Bar& b, std::size_t n) {
	return ((2*n-1)/2) * phase_speed(b)/b.length * std::numbers::pi; 
}

//Mass loaded bar with mass at x=L and free at x=0

inline double bar_mass 
};