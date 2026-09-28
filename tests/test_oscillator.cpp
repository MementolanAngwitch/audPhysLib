#include "oscillators.hpp"
#include "test_util.hpp"
#include <cassert>

int main() {
    // s = 400 N/m, m = 1 kg  →  ω0 = 20 rad/s
    assert(close(audphys::natural_angular_freq(400.0, 1.0), 20.0));
    // Check conversions between omega and frequency
}