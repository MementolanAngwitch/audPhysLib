#include "bars.hpp"

BarModes(Bar bar_, double f_max_) : bar(bar_), f_max(f_max_) {
	double upsilon2 = 2.0 * (1.0 + get_poisson_ratio()) / bar.get_shear();
	double a_c = 1/bar.get_kappa() * std::sqrt(1.0 / (upsilon * upsilon) + 1.0);
	double s2 = L*L*A/I;
	double c = std::sqrt(bar.get_young_mod() / bar.get_rho() );
	double wbar_c = s / std::sqrt(upsilon2);
	double wbar_max = f_max * std::numbers::pi * 2.0 * bar.get_length() / c;


}

Dispersion BarModes::dispersion(double wbar, double s2, double g2) {
	D = std::sqrt(std::pow(wbar,4) * (1.0-g2) * (1.0-g2) + 4.0 * wbar * wbar * s2);
	double a2 = 0.5*(wbar*wbar*(1.0+g2) + D);
	double b2 = 0.5*(-1.0*wbar*wbar*(1.0+g2) + D);
	double a = std::sqrt(a2);
	if (b2> 0) {
		return {a, std::sqrt(b), False}; 
	}
	else{
		return {a, std::complex(0,std::sqrt(b)), True};
	}
}