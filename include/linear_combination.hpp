#pragma once
#include "forced_oscillators.hpp"
#include <vector>

namespace audphys{
//Linear Combination
//Same angular frequencies
// x = A cos(wt+phi)
// A = [(SUM(An*cos(phi_n))^2)^2 + (SUM(A_n*sin(phi_n)))^2]^(1/2)
// tan(phi) = SUM(A_n * sin(phi_n)) / SUM(A_n * cos(phi_n))

inline AmplitudePhase combined_amp_phase(const std::vector<AmplitudePhase>& parts) {
	double A = 0;
	double cosA = 0;
	double sinA = 0;
	double phi = 0;
	for (const auto& p : parts) {
		cosA += p.A * std::cos(p.phi);
		sinA += p.A * std::sin(p.phi);
	}
	A = std::sqrt( cosA * cosA + sinA * sinA);
	phi = atan2(sinA, cosA);
	return {A, phi};
}

// different angular frequenices
// w2 = w1 + delta(w)
// x = Ae^(j(w1*t + phi))
// A = [A1^2 +A2^2 + 2A1A2 * cos(phi1 - phi2 - delta*w*t)]^1/2
// tan(phi) = (A1 * sin(phi_1) + A2 * sin(phi_2 + delta*w*t))/(A1 * cos(phi_1) + A2 cos(phi_2 + delta*w*t))

inline AmplitudePhase combine_amp_phase(const AmplitudePhase& lhs, const Drive& lhd, const AmplitudePhase& rhs, const Drive& rhd, double t) {
	double w1 = lhd.w;
	double deltaw = std::abs(lhd.w - rhd.w);
	double w2 = w1 + deltaw;
	double A1 = lhs.A;
	double A2 = rhs.A;
	double phi1 = lhs.phi;
	double phi2 = rhs.phi;

	double combA = std::sqrt(A1*A1 + A2*2 + 2 * A1 * A2 * std::cos(phi1 - phi2 - deltaw * t));
	double phi = atan2( (A1*std::sin(phi1) + A2 * std::sin(phi2 + deltaw*t)),(A1*std::cos(phi1) + A2 * std::cos(phi2 + deltaw*t)) );
	return {combA, phi};
}
}