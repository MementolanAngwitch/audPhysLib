#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <cstdint>
#include <cstddef>

namespace audphys{
template <class F>
std::vector<double> sample(F f, double sample_rate, std::size_t N) {
	std::vector<double> out(N);
	for (std::size_t n = 0; n < N; ++n) out[n] = f(n / sample_rate);
	return out;
}
void write_le(std::ofstream& file, std::int64_t value, int bytes);
bool write_mono_wav(const std::string fileName, int sample_rate, int bit_depth, const std::vector<double>& data);
}