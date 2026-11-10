#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <cstdint>
#include <cstddef>

namespace audphys{
template <class F>
std::vector<float> sample(F f, float sample_rate, std::size_t N) {
	std::vector<float> out(N);
	for (std::size_t n = 0; n < N; ++n) out[n] = f(n / sample_rate);
	return out;
}
void write_le(std::ofstream& file, std::int64_t value, int bytes);
bool write_mono_wav(const std::string fileName, int sample_rate, int bit_depth, const std::vector<float>& data);

}