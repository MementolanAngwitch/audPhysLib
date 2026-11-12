#include <cstdio>
#include <cmath>
#include <algorithm>
#include "strings.hpp"
#include "export.hpp"
#include "analytic.hpp"
 
using namespace audphys;
using namespace audphys::strings;
 
// g++ demo.cpp -Iinclude -std=c++20 .\include\strings.cpp .\include\oscillators.cpp .\include\export.cpp .\include\input.cpp
#if 1
// Tests generated with Claude 
// A steel guitar high-E string (0.25 mm diameter, 650 mm scale, tuned to ~329.6 Hz)
constexpr float L   = 0.65f;
constexpr float r   = 0.000125f;
constexpr float rho = 7850.0f;                                   // steel [kg/m^3]
const     float mu  = rho * std::numbers::pi_v<float> * r * r;   // [kg/m]
const     float T   = mu * std::pow(2.0f * L * 329.6f, 2.0f);    // T = mu c^2, c = 2 L f1
 
constexpr float fs = 48000.0f;
constexpr float dt = 1.0f / fs;
 
static int failures = 0;
static void check(bool ok, const char* what) {
    std::printf("  [%s] %s\n", ok ? "PASS" : "FAIL", what);
    if (!ok) ++failures;
}
 
static void normalize(std::vector<float>& v) {
    float peak = 0.0f;
    for (float x : v) peak = std::max(peak, std::abs(x));
    if (peak > 0.0f) for (float& x : v) x *= 0.9f / peak;
}
 
// ---------------------------------------------------------------------------
// 1. The String itself: derived quantities should agree with each other.
static void test_string_basics() {
    std::puts("\n1. String basics");
    String s(L, mu, T, r, {2e11f, 0.0f});
    const float f1 = s.fundamental() / (2.0f * std::numbers::pi_v<float>);
    std::printf("  mu = %.3e kg/m, T = %.1f N, c = %.1f m/s, f1 = %.2f Hz\n", mu, T, s.wavespeed(), f1);
    check(std::abs(f1 - 329.6f) < 0.1f, "fundamental() gives the tuned pitch");
    check(std::abs(s.get_rho() - rho) / rho < 1e-4f, "get_rho() recovers volumetric density from mu and r");
    check(std::abs(s.check_length(-0.1f) - 0.1f) < 1e-6f, "check_length folds x = -0.1 to 0.1");
    check(std::abs(s.check_length(0.7f) - 0.6f) < 1e-6f, "check_length folds x = L+0.05 to L-0.05");
    // Two lengths past the end: a fold should land back inside [0, L]
    const float far = s.check_length(1.5f);
    check(far >= 0.0f && far <= L, "check_length(1.5) lands inside [0, L]");
}
 
// ---------------------------------------------------------------------------
// 2. Undamped pluck vs the analytic (Fourier series) PluckedString.
//    Removing every loss: zero air density, real Young's modulus, zero support conductance.
//    Both use the same modes, so the only difference left is the time stepping.
static void test_pluck_vs_analytic() {
    std::puts("\n2. Undamped pluck vs analytic PluckedString");
    const Air vacuum(0.0f, 1.5e-5f);
    String s(L, mu, T, r, {2e11f, 0.0f});
    const float x0 = 0.2f * L, h = 0.002f, xl = 0.37f * L;
    auto triangle = [&](float x) { return x < x0 ? h * x / x0 : h * (L - x) / (L - x0); };
 
    // The recurrence is exact, so the error should be rounding-level at any step size.
    float err_at[2] = {0.0f, 0.0f};
    const int subs[2] = {1, 8};
    for (int k2 = 0; k2 < 2; ++k2) {
        const float h_dt = dt / subs[k2];
        PlayString ps(s, vacuum, 0.0f, 5000.0f, h_dt, xl);
        PluckedString ref(vacuum, s, ps.modes(), triangle);
        ps.pluck(x0, h);
        float max_err = 0.0f, max_ref = 0.0f;
        const std::size_t N = static_cast<std::size_t>(0.05f / h_dt);  // 50 ms
        for (std::size_t k = 0; k < N; ++k) {
            const float a = ps.output(), b = ref.displacement(xl, k * h_dt);
            max_err = std::max(max_err, std::abs(a - b));
            max_ref = std::max(max_ref, std::abs(b));
            ps.step(h_dt);
        }
        err_at[k2] = max_err / max_ref;
        std::printf("  dt = 1/%-6.0f %zu modes, max error over 50 ms = %.2f%% of peak\n",
                    fs * subs[k2], ps.modes(), 100.0f * err_at[k2]);
    }
    check(err_at[1] < 0.001f, "pluck matches the analytic series at a fine step");
    check(err_at[0] < 0.001f, "stepper is accurate at the audio rate");
}
 
// ---------------------------------------------------------------------------
// 3. Damping times per mode. Prints each mechanism separately.
static void test_damping_table() {
    std::puts("\n3. Decay times per mode [s]");
    const Air air;
    String s(L, mu, T, r, {2e11f, 2e11f * 2e-4f});
    const float G = 1e-3f;
    std::puts("   n    f [Hz]   tau_air  tau_int  tau_sup  tau_total");
    float sup_first = 0.0f, sup_last = 0.0f;
    for (int n : {1, 2, 4, 8, 16}) {
        const float w = n * s.fundamental();
        const float sup = tau_support(s, G);
        if (n == 1) sup_first = sup;
        sup_last = sup;
        std::printf("  %2d  %8.1f  %7.3f  %7.3f  %7.3f  %7.3f\n", n, w / (2 * std::numbers::pi_v<float>),
                    tau_air(air, s, w), tau_internal(s, w), sup, tau_total(air, s, w, G));
    }
    // With constant G, support losses are the same for every mode.
    check(sup_first == sup_last, "tau_support is mode-independent for constant G");
}
 
// ---------------------------------------------------------------------------
// 4. Stability of the time stepper. The old semi-implicit Euler blew up above fs / pi
//    (~15.3 kHz at 48 kHz). The exact recurrence should stay bounded for any f_max.
static void test_stability() {
    std::puts("\n4. Stepper stability vs f_max");
    const Air vacuum(0.0f, 1.5e-5f);
    String s(L, mu, T, r, {2e11f, 0.0f});
    for (float f_max : {10000.0f, 15000.0f, 16000.0f, 20000.0f}) {
        PlayString ps(s, vacuum, 0.0f, f_max, dt, 0.37f * L);
        ps.pluck(0.2f * L, 0.002f);
        float peak = 0.0f;
        for (std::size_t k = 0; k < static_cast<std::size_t>(0.2f * fs); ++k) {
            ps.step(dt);
            peak = std::max(peak, std::abs(ps.output()));
        }
        std::printf("  f_max = %5.0f Hz -> %3zu modes, peak |y| = %.3e m %s\n",
                    f_max, ps.modes(), peak, (peak > 0.01f || !std::isfinite(peak)) ? "<-- blew up" : "");
    }
}
 
// ---------------------------------------------------------------------------
// 5. Strike: an impulse J sets modal velocities; for a short while the contact
//    point should move at roughly the same speed it was given, then reflections arrive.
static void test_strike() {
    std::puts("\n5. Strike");
    const Air air;
    String s(L, mu, T, r, {2e11f, 2e11f * 2e-4f});
    PlayString ps(s, air, 1e-3f, 12000.0f, dt, 0.37f * L);
    ps.strike(0.12f * L, 1e-4f);  // 0.1 mN*s
    check(std::abs(ps.contact_displacement()) < 1e-12, "string is still flat right after the impulse");
    ps.step(dt);
    check(ps.contact_displacement() > 0.0, "contact point moves in the direction of the impulse");
}
 
// ---------------------------------------------------------------------------
// 6. Rendered audio: a pluck and a felt-hammer strike, written to WAV.
static void render_audio() {
    std::puts("\n6. Audio");
    const Air air;
    String s(L, mu, T, r, {2e11f, 2e11f * 2e-4f});
    const std::size_t N = static_cast<std::size_t>(3.0f * fs);
 
    PlayString plucked(s, air, 1e-3f, 12000.0f, dt, 0.37f * L);
    plucked.pluck(0.2f * L, 0.002f);
    auto pl = render(plucked, N, dt);
    normalize(pl);  // output is in metres (~mm); the WAV writer clamps to [-1, 1]
    check(write_mono_wav("pluck.wav", static_cast<int>(fs), 16, pl), "wrote pluck.wav");
 
    PlayString hammered(s, air, 1e-3f, 12000.0f, dt, 0.37f * L);
    hammered.set_contact(0.12f * L);  // Hammer never sets its own strike point (see notes)
    Hammer hammer(0.002f, 1e9f, 2.5f, 1.0f);
    std::vector<float> hm(N);
    float peakF = 0.0f; std::size_t contact_steps = 0;
    for (std::size_t k = 0; k < N; ++k) {
        const float F = hammer.step(hammered, dt);  // NB: this also steps the string
        peakF = std::max(peakF, F);
        if (F > 0.0f) ++contact_steps;
        hm[k] = hammered.output();
    }
    std::printf("  hammer: peak force %.2f N, in contact for %zu samples (%.2f ms)\n",
                peakF, contact_steps, 1000.0f * contact_steps * dt);
    normalize(hm);
    check(write_mono_wav("hammer.wav", static_cast<int>(fs), 16, hm), "wrote hammer.wav");
}
 
int main() {
    test_string_basics();
    test_pluck_vs_analytic();
    test_damping_table();
    test_stability();
    test_strike();
    render_audio();
    std::printf("\n%d check(s) failed\n", failures);
    return failures == 0 ? 0 : 1;
}
#endif