#include <cstdio>
#include <cmath>
#include <algorithm>
#include "strings.hpp"
#include "export.hpp"
#include "analytic.hpp"
 
using namespace audphys;
using namespace audphys::strings;
 
// g++ demo.cpp -Iinclude -std=c++20 .\include\strings.cpp .\include\oscillators.cpp .\include\export.cpp .\include\input.cpp
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
        PlayString ps(s, vacuum, 2.0f * std::numbers::pi_v<float> * 5000.0f, h_dt, xl);
        PluckedString ref(vacuum, s, ps.modes(), triangle);
        ps.pluck(x0, h);
        float max_err = 0.0f, max_ref = 0.0f;
        const std::size_t N = static_cast<std::size_t>(0.05f / h_dt);  // 50 ms
        for (std::size_t k = 0; k < N; ++k) {
            const float a = ps.output(), b = ref.displacement(xl, k * h_dt);
            max_err = std::max(max_err, std::abs(a - b));
            max_ref = std::max(max_ref, std::abs(b));
            ps.step();
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
    std::puts("\n3. Decay times per mode [s] (resistive bridge at x = L)");
    const Air air;
    const float G = 1e-3f;  // conductance at the bridge [m/(N s)]
    String s(L, mu, T, r, {2e11f, 2e11f * 2e-4f}, Support::rigid(), Support(0.0f, 0.0f, 1.0f / G));
    std::puts("   n    f [Hz]   tau_air  tau_int  tau_sup  tau_total");
    for (std::size_t n : {1, 2, 4, 8, 16}) {
        const float w = mode_frequency(s, n);
        std::printf("  %2zu  %8.1f  %7.3f  %7.3f  %7.3f  %7.3f\n", n, w / (2 * std::numbers::pi_v<float>),
                    tau_air(air, s, w), tau_internal(s, w), tau_support(s, n), tau_total(air, s, n));
    }
    // One resistive end: tau = L / (T G)
    check(std::abs(tau_support(s, 3) - L / (T * G)) / (L / (T * G)) < 1e-3f, "resistive end gives tau = L / (T G)");
}

// ---------------------------------------------------------------------------
// 3b. Mass-loaded end vs the exact equation cot(kL) = (m / m_s) kL, solved with brent.
static void test_mass_loaded() {
    std::puts("\n3b. Mass-loaded end vs exact cot(kL) = (m/m_s) kL");
    const float m = 0.01f, ms = mu * L, ratio = m / ms;
    String s(L, mu, T, r, {2e11f, 0.0f}, Support::rigid(), Support(m, 0.0f, 0.0f));
    auto f = [ratio](float x) { return std::cos(x) - ratio * x * std::sin(x); };  // no poles
    float worst = 0.0f;
    for (std::size_t n : {1, 2, 3, 5, 10}) {
        const float pi = std::numbers::pi_v<float>;
        const float exact = s.wavespeed() * brent(f, n * pi, n * pi + 0.5f * pi, 1e-7f, 100) / L;
        const float model = mode_frequency(s, n);
        const float rigid = n * s.fundamental();
        const float cents_model = 1200.0f * std::log2(model / rigid), cents_exact = 1200.0f * std::log2(exact / rigid);
        std::printf("  n = %2zu  shift vs rigid: model %+6.3f cents, exact %+6.3f cents\n", n, cents_model, cents_exact);
        worst = std::max(worst, std::abs(cents_model - cents_exact));
    }
    check(worst < 0.1f, "perturbation model within 0.1 cent of the exact equation");
    bool threw = false;
    try { String soft(L, mu, T, r, {2e11f, 0.0f}, Support::rigid(), Support(1e-5f, 0.0f, 0.0f)); mode_frequency(soft, 1); }
    catch (const std::invalid_argument&) { threw = true; }
    check(threw, "a support lighter than the string is rejected");
}

// ---------------------------------------------------------------------------
// 4. Stability of the time stepper. The exact recurrence should stay bounded for any
//    w_max; modes are also capped at Nyquist.
static void test_stability() {
    std::puts("\n4. Stepper stability vs f_max");
    const Air vacuum(0.0f, 1.5e-5f);
    String s(L, mu, T, r, {2e11f, 0.0f});
    for (float f_max : {10000.0f, 15000.0f, 16000.0f, 20000.0f}) {
        PlayString ps(s, vacuum, 2.0f * std::numbers::pi_v<float> * f_max, dt, 0.37f * L);
        ps.pluck(0.2f * L, 0.002f);
        float peak = 0.0f;
        for (std::size_t k = 0; k < static_cast<std::size_t>(0.2f * fs); ++k) {
            ps.step();
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
    PlayString ps(s, air, 2.0f * std::numbers::pi_v<float> * 12000.0f, dt, 0.37f * L);
    ps.strike(0.12f * L, 1e-4f);  // 0.1 mN*s
    check(std::abs(ps.contact_displacement()) < 1e-12, "string is still flat right after the impulse");
    ps.step();
    check(ps.contact_displacement() > 0.0, "contact point moves in the direction of the impulse");
}

// ---------------------------------------------------------------------------
// 6. Rendered audio: a pluck and a felt-hammer strike, written to WAV.
static void render_audio() {
    std::puts("\n6. Audio");
    const Air air;
    String s(L, mu, T, r, {2e11f, 2e11f * 2e-4f});
    const std::size_t N = static_cast<std::size_t>(3.0f * fs);

    PlayString plucked(s, air, 2.0f * std::numbers::pi_v<float> * 12000.0f, dt, 0.37f * L);
    plucked.pluck(0.2f * L, 0.002f);
    auto pl = render(plucked, N);
    normalize(pl);  // output is in metres (~mm); the WAV writer clamps to [-1, 1]
    check(write_mono_wav("pluck.wav", static_cast<int>(fs), 16, pl), "wrote pluck.wav");

    PlayString hammered(s, air, 2.0f * std::numbers::pi_v<float> * 12000.0f, dt, 0.37f * L);
        Hammer hammer(0.002f, 1e9f, 2.5f, 1.0f, 0.12f * L);
    std::vector<float> hm(N);
    float peakF = 0.0f; std::size_t contact_steps = 0;
    for (std::size_t k = 0; k < N; ++k) {
        const float F = hammer.step(hammered);  // NB: this also steps the string
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
    test_mass_loaded();
    test_stability();
    test_strike();
    render_audio();
    std::printf("\n%d check(s) failed\n", failures);
    return failures == 0 ? 0 : 1;
}