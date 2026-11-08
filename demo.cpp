#include "include/oscillators.hpp"
#include "include/export.hpp"
#include "include/strings.hpp"
#include "include/tools.cpp"
#include "include/input.hpp"
#include <iostream>
// g++ demo.cpp -Iinclude -std=c++20 .\include\strings.cpp .\include\oscillators.cpp .\include\export.cpp .\include\input.cpp
using namespace audphys;

int main() {
    // Oscillators
    oscillators::Oscillator osc(0.01, 4000.0, 0.0);        // m [kg], K [N/m], R
    oscillators::Impulse hit(0.0, 0.0, 1.0);               // no force -> free oscillation
    oscillators::WaveMotion wave(osc, hit, 1.0, 0.0, 0.0); // x0 = 1, v0 = 0, a0 = 0

    std::cout << "x(0) = " << wave.displacement(0.0) << '\n';

    auto data = sample([&](double t) { return static_cast<float>(wave.displacement(t)); }, 48000, 48000);
    bool ok = write_mono_wav("out.wav", 48000, 16, data);
    std::cout << (ok ? "wrote out.wav\n" : "write failed\n");

    //String
    // using BC = strings::BoundaryCondition;
    // strings::String string(1.0,0.1,1.0,0.1,BC::Fixed,BC::Fixed); 
    // strings::PluckedString pluckedString(string, 10, [](double t){return t;});
    // return 0;
}