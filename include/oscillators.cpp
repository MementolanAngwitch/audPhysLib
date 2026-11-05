#include "oscillators.hpp"

using namespace::audphys{

Oscillator::Oscillator(EnclosedAir a) {

}

WaveMotion::WaveMotion(Oscillator o,double x0_, double v0_, double a0_) : osc(o), x0(x0_), v0(v0_), a0(a0_){
	if (o.R != 0) damped = true; // allowing negative damping
	if (o.F != 0) forced = true; // allowing negative force
	if (!damped) {
		w0 = std::sqrt(o.K, o.m);
		A = std::sqrt(x0*x0 + (v0/w0) * (v0/w0));
		phi = std::atan2(-v0/(w0*x0));
	}
	if (damped) {
		w0 = std::sqrt(o.K, o.m);
		w0 = std::sqrt(w0*w0 - alpha()*alpha());
		A = x0;
	}
	if (forced) {
		if (F.w == w0) {
			A = quality_factor() * F.F / (o.m * w0 * w0);
		}
		A = F.F / (F.w * mechanical_impedance());
	}
}


WaveMotion::alpha() const {
	return o.R / (2.0 * o.m);
}

std::complex<double> WaveMotion::mechnaical_reactance() const {
	return F.w * o.m - o.K / F.w; 
}
std::complex<double> WaveMotion::mechanical_impedance() const {
	return { o.R, mechanical_reactance();};
}

double WaveMotion::period(double f) const {
	return 1.0/f;
}


}