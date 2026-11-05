#pragma once

namespace audphys{
template <class F>
std::vector<double> sample(F f, double sample_rate, std::size_t N);
void write_le(std::ofstream& file, std::int64_t value, int bytes);
bool write_mono_wav(const std::string fileName, int sample_rate, int bit_depth, const std::vector<double>& data)
}