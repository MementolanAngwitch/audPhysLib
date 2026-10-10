#pragma once
#include <cmath>
#include <fstream>
namespace audphys{
template <class F>
std::vector<double> sample(F f, double sample_rate, std::size_t N) {
	std::vector<double> out(N);
	for (std::size_t n = 0; n < N; ++n) out[n] = f(n / sample_rate);
	return out;
}

inline void write_wav(const std::string fileName, int sample_rate = 48000, double dur, int bit_depth = 16, int max_amp, const std::vector<double>& data) {
	//Header chunk
	audioFile << "RIFF";
	audioFile << "----";
	audioFile << "WAVE";

	//Format chunk
	audioFile << "fmt";
	writeToFile(audioFile,16,4); //size
	writeToFile(audioFile, 1,2); //compression code
	writeToFile(audioFile, 1,2); //number of channels
	writeToFile(audioFile, sr,4); //sample rate
	writeToFile(audioFile, sr * bitDepth / 8, 4); //byte rate
	writeToFile(audioFile, bitDepth/8, 2); //block depth
	writeToFile(audioFile, bitDepth, 2); // bit depth

	//Data chunk
	audioFile << "data"; 
	audioFile << "----";

	int preAudioPosition = audioFile.tellp();
	auto maxAmplitude = pow(2, bitDepth - 1) - 1;
	for(int i =0; i < sr * duration; i++) {
		auto sample = data[i];
		int intSample = static_cast<int> (sample * maxAmplitude);
		writeToFile(audioFile, intSample, 2);
	}
	int postAudioPosition = audioFile.tellp();

	audioFile.seekp(preAudioPosition -4); //writing size of audio after "data "
	writeToFile(audioFile, postAudioPosition - preAudioPosition, 4);

	audioFile.seekp(4, std::ios::beg);
	writeToFile(audioFile, postAudioPosition - 8, 4 ); // Header size 

	audioFile.close();
}


};