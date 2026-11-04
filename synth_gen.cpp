// gen: builds the mode tables for the audPhysLib synth.
// For every instrument, note and strike position it records each vibration mode as

#include <cstdio>
#include <cmath>
#include <numbers>
#include <string>
#include <vector>
#include "audphyslib.hpp"

namespace ap = audphys;
constexpr double PI = std::numbers::pi;
constexpr double W_MAX = 2 * PI * 16000;        
struct M { double w, beta, c, s; };

static double midi_w(int note) { return 2 * PI * 440.0 * std::pow(2.0, (note - 69) / 12.0); }

static std::string json(const std::vector<M>& ms) {
    std::string out = "[";
    char buf[96];
    for (std::size_t i = 0; i < ms.size(); ++i) {
        std::snprintf(buf, sizeof buf, "%s[%.6g,%.4g,%.5g,%.5g]", i ? "," : "", ms[i].w, ms[i].beta, ms[i].c, ms[i].s);
        out += buf;
    }
    return out + "]";
}

// ---- plucked string: pickup is the force on the bridge, proportional to the slope at x = 0
static std::vector<M> string_note(double w, double pos) {
    const double L = 0.65, rho = 0.0006;
    const ap::strings::String s(L, ap::strings::tension_for_ang_freq(L, rho, w), rho);
    const double x0 = pos * L, h = 0.003;
    auto shape = [&](double x) { return x < x0 ? h * x / x0 : h * (L - x) / (L - x0); };
    const int N = std::min(40, int(W_MAX / w));
    auto nm = ap::strings::pluck_modes(s, shape, N, 0.0);
    std::vector<M> out;
    for (auto& md : nm.modes) {
        const double wn = ap::strings::mode_ang_freq(s, md.n), k = ap::strings::mode_wave_number(s, md.n);
        const double beta = 0.7 + 4e-5 * wn;
        out.push_back({wn, beta, k * md.A, k * beta * md.A / wn});   // released at rest: B = beta A / w
    }
    return out;
}

// ---- membranes: struck, listened to at a point off the centre
static std::vector<M> membrane_modes_out(std::vector<ap::membranes::CircularMode>& mds, double rl, double thl) {
    std::vector<M> out;
    for (auto& md : mds) {
        const double wd = ap::membranes::damped_angular_frequency(md);
        if (wd > W_MAX) continue;
        out.push_back({wd, md.beta, 0.0, md.A * ap::membranes::mode_shape_polar(md, rl, thl)}); // phi = -pi/2: pure sine
    }
    return out;
}

static std::vector<M> drum_note(double w, double pos) {
    const double a = 0.2, ps = 0.2;
    ap::membranes::Membrane mem{ps, ap::membranes::circular_tension_for_ang_freq(ps, a, w, 0, 1)};
    auto mds = ap::membranes::circular_modes(mem, a, 5, 3);
    for (auto& md : mds) ap::membranes::set_beta(md, 9.0 + 0.004 * md.w);
    ap::membranes::strike(mds, 0.02, pos * a, 0.0);
    return membrane_modes_out(mds, 0.42 * a, 0.5);
}

static std::vector<M> timpani_note(double w, double pos) {
    // pitched by the (1,1) mode, as real timpani are; the kettle stiffens the (0,1) mode
    const double a = 0.33, ps = 0.26, rl = 0.3 * a, thl = 0.5;
    ap::membranes::Membrane mem{ps, ap::membranes::circular_tension_for_ang_freq(ps, a, w, 1, 1)};
    auto mds = ap::membranes::circular_modes(mem, a, 5, 3);
    for (auto& md : mds) ap::membranes::set_beta(md, 2.5 + 0.002 * md.w);
    ap::membranes::strike(mds, 0.02, pos * a, 0.0);
    auto out = std::vector<M>{};
    for (auto& md : mds) {
        if (md.m == 0 && md.n == 1) {
            ap::membranes::AirVolume av{101325.0, 0.12, 1.4, md};
            auto air = ap::membranes::air_loaded_mode(av);
            ap::membranes::set_beta(air, 3.0 * md.beta);   // this mode pumps air and dies quickly
            ap::membranes::strike_air_loaded(air, 0.02, pos * a);
            out.push_back({ap::membranes::damped_angular_frequency(air), air.beta, 0.0, ap::membranes::symmetric_normal_mode(air, rl)});
            continue;
        }
        const double wd = ap::membranes::damped_angular_frequency(md);
        if (wd > W_MAX) continue;
        out.push_back({wd, md.beta, 0.0, md.A * ap::membranes::mode_shape_polar(md, rl, thl)});
    }
    return out;
}

// ---- free-free aluminium bar (glockenspiel), struck at x_s, heard near one end
static std::vector<M> bar_note(double w, double pos) {
    using enum ap::bars::BendingEnd;
    ap::bars::Bar b{1.0, ap::bars::rectangular_section(0.032, 0.008), 6.9e10, 2700, 2.6e10};
    b.length = ap::bars::length_for_ang_freq(b, Free, Free, w, 1);
    const double kc = b.crossSection.kappa * ap::bars::longitudinal_phase_speed(b);
    std::vector<M> out;
    for (double wn : ap::bars::bending_natural_frequency(b, Free, Free, 6)) {
        if (wn > W_MAX) break;
        const double A = ap::bars::strike_amplitude(b, Free, Free, wn, 0.02, pos * b.length);
        const double phi_l = ap::bars::mode_shape(Free, Free, std::sqrt(wn / kc), b.length, 0.03 * b.length);
        out.push_back({wn, 0.9 + 3e-4 * wn, 0.0, A * phi_l});
    }
    return out;
}

// ---- clamped steel plate, struck at r0, heard at r = 0.3a (axisymmetric modes)
static std::vector<M> plate_note(double w, double pos) {
    ap::plates::CircularThinPlate p{7800, 0.3, 2.0e11, 0.0015, 0.2};
    p.a = ap::plates::radius_for_ang_freq(p, w, 1);
    std::vector<ap::plates::PlateMode> mds;
    for (std::size_t n = 1; n <= 8; ++n) {
        auto md = ap::plates::make_mode(p, n);
        if (md.w > W_MAX) break;
        mds.push_back(md);
    }
    ap::plates::strike(mds, 0.02, pos * p.a, p);
    std::vector<M> out;
    for (std::size_t i = 0; i < mds.size(); ++i) {
        const auto& md = mds[i];   // displacement() at t = T/4 with phi = -pi/2 gives A * shape(r)
        const double A_shape = ap::plates::displacement(md, i + 1, 0.3 * p.a, PI / (2.0 * md.w));
        out.push_back({md.w, 1.6 + 5e-4 * md.w, 0.0, A_shape});
    }
    return out;
}

int main() {
    struct Inst { const char* id; const char* name; int base; std::vector<double> pos; std::vector<M> (*fn)(double, double); };
    std::vector<Inst> insts = {
        {"string",  "Plucked string",   48, {0.06, 0.12, 0.2, 0.32, 0.5},  string_note},
        {"drum",    "Hand drum",        43, {0.05, 0.25, 0.45, 0.65, 0.85}, drum_note},
        {"timpani", "Timpani",          40, {0.15, 0.35, 0.55, 0.75, 0.9},  timpani_note},
        {"bar",     "Glockenspiel bar", 72, {0.04, 0.12, 0.22, 0.34, 0.5},  bar_note},
        {"plate",   "Steel plate",      60, {0.0, 0.2, 0.4, 0.6, 0.8},     plate_note},
    };
    std::printf("{\"instruments\":[");
    for (std::size_t i = 0; i < insts.size(); ++i) {
        const auto& in = insts[i];
        std::printf("%s{\"id\":\"%s\",\"name\":\"%s\",\"base\":%d,\"pos\":[", i ? "," : "", in.id, in.name, in.base);
        for (std::size_t j = 0; j < in.pos.size(); ++j) std::printf("%s%g", j ? "," : "", in.pos[j]);
        std::printf("],\"notes\":[");
        for (int k = 0; k < 25; ++k) {
            std::printf("%s[", k ? "," : "");
            for (std::size_t j = 0; j < in.pos.size(); ++j)
                std::printf("%s%s", j ? "," : "", json(in.fn(midi_w(in.base + k), in.pos[j])).c_str());
            std::printf("]");
        }
        std::printf("]}");
    }
    std::printf("]}\n");
}
