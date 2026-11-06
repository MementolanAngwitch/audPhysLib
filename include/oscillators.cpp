#include "oscillators.hpp"

namespace audphys::oscillators {

Oscillator::Oscillator(EnclosedAir a) {
	K = a.gamma * a.Pa * a.S / a.L;
	m = a.m;
};	

Oscillator::Oscillator(HelmholtzResonator r) {
	m = r.rho * r.S * r.L;
	K = r.rho * r.S * r.S * r.c * r.c / r.V;
}

WaveMotion::WaveMotion(Oscillator o, Impulse F_, double x0_, double v0_, double a0_) : x0(x0_), v0(v0_), a0(a0_), F(F_), osc(o) {
	if (o.R != 0) damped = true; // allowing negative damping
	if (F.F != 0) forced = true; // allowing negative force
	if (!damped) {
		w0 = std::sqrt(o.K / o.m);
		A = std::sqrt(x0*x0 + (v0/w0) * (v0/w0));
		phi = std::atan2(-v0/w0, x0);
	}
	if (damped) {
		w0 = std::sqrt(o.K / o.m);
		w0 = std::sqrt(w0*w0 - alpha()*alpha());
		A  = std::sqrt(x0*x0 + std::pow((v0 + alpha()*x0) / w0, 2));
		phi = std::atan2(-(v0 + alpha()*x0) / w0, x0);
	}
	if (forced) {
		if (F.w == w0) {
			A = quality_factor() * F.F / (o.m * w0 * w0);
		}
		A = F.F / (F.w * std::abs(mechanical_impedance()));
	}
}


double WaveMotion::alpha() const {
	return osc.R / (2.0 * osc.m);
}

double WaveMotion::mechanical_reactance() const {
	return F.w * osc.m - osc.K / F.w;
}
std::complex<double> WaveMotion::mechanical_impedance() const {
	return { osc.R, mechanical_reactance() };
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
	const double env = A * std::exp(-alpha() * t);
	const double transient = env * std::cos(w0 * t + phi);
	if (!forced) {
		return transient;
	}
	const double X = F.F / (F.w * std::abs(mechanical_impedance()));
	const double theta = std::arg(mechanical_impedance());
	return transient + X * std::sin(F.w * t - theta);
}

double WaveMotion::velocity(double t) const {
	if (!damped && !forced) {
		return -A * w0 * std::sin(w0 * t + phi);
	}
	const double env = A * std::exp(-alpha() * t);
	const double transient = -env * (alpha() * std::cos(w0 * t + phi) + w0 * std::sin(w0 * t + phi));
	if (!forced) {
		return transient;
	}
	const double X = F.F / (F.w * std::abs(mechanical_impedance()));
	const double theta = std::arg(mechanical_impedance());
	return transient + X * F.w * std::cos(F.w * t - theta);
}

double WaveMotion::acceleration(double t) const {
	if (!damped && !forced) {
		return -A * w0 * w0 * std::cos(w0 * t + phi);
	}
	const double env = A * std::exp(-alpha() * t);
	const double transient = env * ((alpha() * alpha() - w0 * w0) * std::cos(w0 * t + phi)
	                                + 2.0 * alpha() * w0 * std::sin(w0 * t + phi));
	if (!forced) {
		return transient;
	}
	const double X = F.F / (F.w * std::abs(mechanical_impedance()));
	const double theta = std::arg(mechanical_impedance());
	return transient - X * F.w * F.w * std::sin(F.w * t - theta);
}
}