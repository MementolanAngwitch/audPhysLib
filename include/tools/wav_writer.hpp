#pragma once
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

namespace audphys{
template <class F>
std::vector<double> sample(F f, double sample_rate, std::size_t N) {
	std::vector<double> out(N);
	for (std::size_t n = 0; n < N; ++n) out[n] = f(n / sample_rate);
	return out;
}

//Writes the lowest bytes of value, least significant byte first
inline void write_le(std::ofstream& file, std::int64_t value, int bytes) {
	const auto u = static_cast<std::uint64_t>(value);
	for (int i = 0; i < bytes; ++i) {
		file.put(static_cast<char>((u >>(8*i)) & 0xFF));
	}
}
//Writes mono PCM audio, bit_depth must be 16,24 or 32
inline bool write_mono_wav(const std::string fileName, int sample_rate, int bit_depth, const std::vector<double>& data) {
	//Standard: sr=48000, bit_depth = 16

	std::ofstream file(fileName, std::ios::binary);
	if (!file) return false;
	const int bytes_per_sample = bit_depth/8;

	// Header chunk
	file.write("RIFF",4);
	file.write("----",4); //file size - 8
	file.write("WAVE",4);

	//Format chunk
	file.write("fmt",4);
	write_le(file,16,4); //size
	write_le(file, 1,2); //compression code
	write_le(file, 1,2); //number of channels
	write_le(file, sample_rate,4); //sample rate
	write_le(file, sample_rate * bit_depth / 8, 4); //byte rate
	write_le(file, bit_depth/8, 2); //block depth
	write_le(file, bit_depth, 2); // bit depth

	//Data chunk
	file.write("data",4);
	file.write("----",4); //placeholder data size

	const std::int64_t pre = file.tellp();

	const double max_amplitude = std::pow(2.0, bit_depth-1) -1;
	for (double s: data) {
		s = std::clamp(s,-1.0,1.0);
		write_le(file, std::llround(s*max_amplitude), bytes_per_sample);
	}

	const std::int64_t post = file.tellp();

	file.seekp(pre - 4); //data size
	write_le(file, post-pre, 4);
	file.seekp(4, std::ios::beg); // file size
	write_le(file,post-8,4);

	return static_cast<bool>(file);
}


}