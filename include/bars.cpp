#include "bars.hpp"

inline BarModes::BarModes(const Bar& bar_, BoundaryConditions left_, BoundaryConditions right_, double omega_max_)
	: bar(bar_), left(left_), right(right_), omega_max(omega_max_),
	  s2(bar_.get_s2()), u2(bar_.get_upsilon2()), c(bar_.get_wave_speed())
{
	using BC = BoundaryConditions;
	const double wbar_max = omega_max * bar.get_length() / c;

	auto is = [&](BC x, BC y) { return (left == x && right == y) || (left == y && right == x); };

	// Which base case, and whether we keep the symmetric (+1) or antisymmetric (-1) half of it.
	Base base = Base::FF;
	int half = 0;
	if      (is(BC::Free,    BC::Free))    { base = Base::FF; }
	else if (is(BC::Clamped, BC::Clamped)) { base = Base::CC; }
	else if (is(BC::Hinged,  BC::Hinged))  { base = Base::HH; }
	else if (is(BC::Clamped, BC::Free))    { base = Base::CF; }
	else if (is(BC::Sliding, BC::Sliding)) { base = Base::HH; }            // same frequencies as hinged-hinged
	else if (is(BC::Free,    BC::Sliding)) { base = Base::FF; half = +1; }
	else if (is(BC::Free,    BC::Hinged))  { base = Base::FF; half = -1; }
	else if (is(BC::Clamped, BC::Sliding)) { base = Base::CC; half = +1; }
	else if (is(BC::Clamped, BC::Hinged))  { base = Base::CC; half = -1; }
	else if (is(BC::Hinged,  BC::Sliding)) { base = Base::HH; half = +1; }

	std::vector<double> wbars;
	if (half == 0) {
		for (double root : scan_roots(base, s2, u2, wbar_max)) wbars.push_back(wbar_from_a(root, s2, u2));
	} else {
		// Full beam: twice as long, so slenderness 2s and wbar 2w, then keep modes of the right parity.
		const double s2f = 4.0 * s2;
		const BC endf = (base == Base::FF) ? BC::Free : (base == Base::CC) ? BC::Clamped : BC::Hinged;
		for (double af : scan_roots(base, s2f, u2, 2.0 * wbar_max)) {
			const Dispersion df = dispersion_at_a(af, s2f, u2);
			const auto Cf = null_vector(build_F(endf, endf, df, u2));
			double Wmax = 0.0, Pmax = 0.0;
			for (int i = 0; i <= 200; ++i) {
				const auto v = mode_eval(Cf, df, u2, i / 200.0);
				Wmax = std::max(Wmax, std::abs(v[0]));
				Pmax = std::max(Pmax, std::abs(v[1]));
			}
			const auto mid = mode_eval(Cf, df, u2, 0.5);
			const bool symmetric = std::abs(mid[1]) / Pmax < std::abs(mid[0]) / Wmax; // Psi(1/2) = 0
			if ((half > 0) == symmetric) wbars.push_back(0.5 * wbar_from_a(af, s2f, u2));
		}
	}

	std::sort(wbars.begin(), wbars.end());
	for (double wb : wbars) add_mode(wb);
}

// Mode shape for the requested ends at a known root
inline void BarModes::add_mode(double wbar_n) {
	const Dispersion d = dispersion(wbar_n, s2, u2);
	const Mat4 F = build_F(left, right, d, u2);
	auto Cn = null_vector(F);
	const double r = residual(F, Cn);

	const double m = simpson([&](double x) {
		const auto v = mode_eval(Cn, d, u2, x);
		return v[0] * v[0] + v[1] * v[1] / s2;
	}, 0.0, 1.0, std::size_t(400));
	for (double& cj : Cn) cj /= std::sqrt(m);

	a.push_back(d.a);
	b.push_back(d.b);
	above.push_back(d.above);
	w_bar.push_back(wbar_n);
	w.push_back(wbar_n * c / bar.get_length());
	res.push_back(r);
	C.push_back(Cn);
}