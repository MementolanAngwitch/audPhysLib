template<class F, class G>
requires std::invocable<F, float> && std::invocable<G, float>
class TransverseWaves { // Takes infinite string
public:
	TransverseWaves(String s_, F f_, G g_) : s(s_), f(f_), g(g_) {}
	float displacement(float x, float t) const {
		const float c = s.wavespeed();
		return f(x - c * t) + g(x + c * t);
	}
private:
	String s;
	F f; // Right-bound wave
	G g; // Left-bound wave
};

// Antiquated, keep for reference and solution checking
template<class Shape>
requires std::invocable<Shape, float>
class PluckedString { // Initial displacement, zero initial velocity
public:
	PluckedString(Air air_, String s_, std::size_t n_, Shape shape_) : air(air_), s(s_), n(n_), shape(shape_) {
		if (!(s.get_right() == BoundaryCondition::Fixed && s.get_left() == BoundaryCondition::Fixed))
			throw std::invalid_argument("PluckedString requires fixed ends");
		make_modes();
	}
	void make_modes() {
		const float L = s.get_length();
		for (std::size_t i = 0; i < n; ++i) {
			const float k = static_cast<float>(i + 1); // mode number 1..n
			B.push_back((2.0 / L) * simpson([&](float x) {
				return shape(x) * std::sin(k * std::numbers::pi * x / L);
			}, 0.0, L));
			A.push_back(0.0); // A_n = 0 for zero initial velocity
		}
	}
	float mode_displacement(std::size_t i, float x, float t) const {
		x = s.check_length(x);
		const float w_i = (i + 1) * std::numbers::pi * s.wavespeed() / s.get_length();
		return (A[i] * std::sin(w_i * t) + B[i] * std::cos(w_i * t) ) * std::sin( (w_i*x)/s.wavespeed());
	}
	float displacement(float x, float t) const {
		x = s.check_length(x);
		float sum = 0.0;
		for (std::size_t i = 0; i < n; ++i) sum += mode_displacement(i, x, t);
		return sum;
	}
private:
	Air air;
	String s;
	std::size_t n; // # of modes. TODO: currently an input
	Shape shape;   // Displacement function y(x,0)
	std::vector<float> A; // Mode contribution from y'
	std::vector<float> B; // Mode contribution from y
};

// Antiquated, keep for reference and solution checking
template<class Shape>
requires std::invocable<Shape, float>
class StruckString { // Zero initial displacement, initial velocity
public:
	StruckString(String s_, std::size_t n_, Shape shape_) : s(s_), n(n_), shape(shape_) {
		if (!(s.get_right() == BoundaryCondition::Fixed && s.get_left() == BoundaryCondition::Fixed))
			throw std::invalid_argument("StruckString requires fixed ends");
		make_modes();
	}
	void make_modes() {
		const float L = s.get_length();
		for (std::size_t i = 0; i < n; ++i) {
			const float k = static_cast<float>(i + 1); // mode number 1..n
			const float w_i = k * std::numbers::pi * s.wavespeed() / L;
			A.push_back((2.0 / (L * w_i)) * simpson([&](float x) {
				return shape(x) * std::sin(k * std::numbers::pi * x / L);
			}, 0.0, L));
			B.push_back(0.0); // B_n = 0 for zero initial displacement
		}
	}
	float mode_displacement(std::size_t i, float x, float t) const {
		x = s.check_length(x);
		const float w_i = (i + 1) * std::numbers::pi * s.wavespeed() / s.get_length();
		return (A[i] * std::sin(w_i * t) + B[i] * std::cos(w_i * t)) * std::sin(w_i * x / s.wavespeed());
	}
	float displacement(float x, float t) const {
		x = s.check_length(x);
		float sum = 0.0;
		for (std::size_t i = 0; i < n; ++i) sum += mode_displacement(i, x, t);
		return sum;
	}
private:
	String s;
	std::size_t n; // # of modes. TODO: currently an input
	Shape shape;   // Velocity function y'(x,0)
	std::vector<float> A; // Mode contribution from y'
	std::vector<float> B; // Mode contribution from y
};
/// Antiquated, keep for reference and solution checking
class HammerString {
public:
	HammerString(String s_, float v0_, float M_, float x0_) : s(s_), v0(v0_), M(M_), x0(x0_) {}
	float contact_displacement(float t) const;
	float displacement(float x , float t) const; // TODO: only true for long strings
private:
	String s;
	float v0; // Initial velocity [m/s]
	float M;  // Hammer mass [kg]
	float x0; // Strike point [m]
};

// class BowedString {
//     TODO
// };
