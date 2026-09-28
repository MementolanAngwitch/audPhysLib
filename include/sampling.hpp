#pragma once

namespace audphys{

template <class F>
std::vector<double> sample(F f, double sample_rate, std::size_t N) {
	std::vector<double> out(N);
	for (std::size_t n = 0; n < N; ++n) out[n] = f(n / sample_rate);
	return out;
}

}