#include "include/oscillators.hpp"
#include "include/export.hpp"
#include "include/strings.hpp"
#include "include/tools.cpp"
#include "include/input.hpp"
#include <iostream>
// g++ demo.cpp -Iinclude -std=c++20 .\include\strings.cpp .\include\oscillators.cpp .\include\export.cpp .\include\input.cpp

using namespace audphys;
using namespace audphys::strings;
 
// a steel low-E guitar string
constexpr float L    = 0.65f;     // length [m]
constexpr float r    = 0.00053f;  // radius [m]
constexpr float rho  = 7850.0f;   // steel density [kg/m^3]
constexpr float T    = 70.0f;     // tension [N]
constexpr float E1   = 2.0e11f;   // Young's modulus, elastic part [N/m^2]
constexpr float loss = 1.0e-4f;   // E2 / E1: internal friction of the steel
constexpr float G    = 1.0e-3f;   // bridge conductance [s/kg], a stand-in value
 
constexpr int   sr    = 48000;    // sample rate [Hz]
constexpr float dt    = 1.0f / sr;
constexpr float secs  = 3.0f;
constexpr float f_max = 5000.0f;  // highest mode kept [Hz]; keeps w*dt well below 2
 

void save(const char* path, std::vector<float> x) {
    float peak = 0.0f;
    for (float v : x) peak = std::max(peak, std::abs(v));
    if (peak > 0.0f) for (float& v : x) v *= 0.9f / peak;
    std::printf("%s %s\n", write_mono_wav(path, sr, 16, x) ? "wrote" : "FAILED to write", path);
}
 
int main() {
    const float mu = rho * std::numbers::pi_v<float> * r * r; // linear density [kg/m]
    const String s(L, mu, T, r, {E1, E1 * loss});
    const Air air(1.2f, 1.5e-5f);
 
    const float w1 = std::numbers::pi_v<float> * s.wavespeed() / L; // fundamental [rad/s]
    const std::size_t n = static_cast<std::size_t>(f_max / (w1 / (2.0f * std::numbers::pi_v<float>)));
    const std::size_t N = static_cast<std::size_t>(secs * sr);
    std::printf("fundamental %.1f Hz, %zu modes, fundamental decay time %.2f s\n",
                w1 / (2.0f * std::numbers::pi_v<float>), n, tau_total(air, s, w1, G));
 
    const float x_listen = 0.9f * L; // a pickup near the bridge end
    const float x0       = 0.2f * L; // where every gesture touches the string
 
    // 1. Pluck: pull the string 3 mm and let go.
    {
        PlayString str(s, air, G, n, x_listen);
        str.pluck(x0, 0.003f);
        std::vector<float> out(N);
        for (std::size_t k = 0; k < N; ++k) { str.step(dt); out[k] = str.output(); }
        save("pluck.wav", out);
    }
 
    // 2. Strike: an instant 2 mN·s tap
    {
        PlayString str(s, air, G, n, x_listen);
        str.strike(x0, 0.002f);
        std::vector<float> out(N);
        for (std::size_t k = 0; k < N; ++k) { str.step(dt); out[k] = str.output(); }
        save("strike.wav", out);
    }
 
    // 3. Hammer: a 3 g felt hammer at 2 m/s. 
    {
        PlayString str(s, air, G, n, x_listen);
        str.set_contact(x0);
        Hammer hammer(0.003f, 1.0e9f, 2.5f, 2.0f);
        std::vector<float> out(N);
        std::size_t contacts = 0, last = 0;
        bool touching = false;
        for (std::size_t k = 0; k < N; ++k) {
            const float F = hammer.step(str, dt);
            if (F > 0.0f) { last = k; if (!touching) ++contacts; }
            touching = F > 0.0f;
            out[k] = str.output();
        }
        std::printf("hammer: %zu contact(s), last at %.2f ms; reflection returns after %.2f ms\n",
                    contacts, 1e3f * last * dt, 1e3f * 2.0f * x0 / s.wavespeed());
        save("hammer.wav", out);
    }
}
 