#include "oscillators.hpp"

using namespace::audphys{

Oscillator::Oscillator(EnclosedAir a) {
	K = a.gamma * a.Pa * a.S / a.L;
	m = a.m;
};	

Oscillator::Oscillator(HelmholtzResonator r) {
	m = r.rho * r.S * r.L;
	K = r.rho * r.S * r.S * r.c * r.c / (r.S * r.L);
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

double WaveMotion::get_A() const { return A; }

double WaveMotion::get_phi() const { return phi; }

double WaveMotion::quality_factor() const { return w0 / (2.0 * alpha()); }

double WaveMotion::displacement(double t) const {
	if (!damped && !forced) {
		return A * std::cos(w0 * t + phi);
	}
	if (damped && !forced) {
		return A * std::exp(-alpha() * t) * std::cos(w0 * t + phi);
	}
	if (forced) {
		return A * std::exp(-alpha() * t) * std::cos(w0 * t + phi) + F.F /(F.w * mechanical_impedance()) * std::sin(F.w * t + phi);
	}
}

double WaveMotion::velocity(double t) const {
	f (!damped && !forced) {
		return -A * w0 * std::sin(w0 * t + phi);
	}
	if (damped && !forced) {
		return A * alpha() * w0 * std::exp(-alpha() * t) * std::sin(w0 * t + phi);
	}
	if (forced) {
		return A * alpha() * w0 * std::exp(-alpha() * t) * std::sin(w0 * t + phi) + F.F /(F.w * mechanical_impedance()) * F.w *std::cos(F.w * t + phi);
	}
}

double WaveMotion::acceleration(double t) const {
	f (!damped && !forced) {
		return -A * w0 * w0 * std::cos(w0 * t + phi);
	}
	if (damped && !forced) {
		return A * alpha() * alpha() * w0 * w0 * std::exp(-alpha() * t) * std::cos(w0 * t + phi);
	}
	if (forced) {
		return A * alpha() * alpha() * w0 * w0* std::exp(-alpha() * t) * std::cos(w0 * t + phi) + F.F /(F.w * mechanical_impedance()) * F.w * F.w *std::sin(F.w * t + phi);
	}
}

}