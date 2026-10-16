#pragma once
namespace audphys {
// ∫[a,b] f(x) dx ≈ h/3 [ f(x0) + 4·Σ_{j=1..n} f(x_{2j-1}) + 2·Σ_{j=1..n-1} f(x_{2j}) + f(x_{2n}) ]
// with h = (b-a)/(2n), x_j = a + j·h
template<class F>
inline double simpson(F f, double a, double b, int n) {
    double h = (b - a) / (2.0 * n);
    double oddSum = 0.0, evenSum = 0.0;
    for (int j = 1; j <= n; ++j) {
        oddSum += f(a + (2.0 * j - 1.0) * h);           // x_1, x_3, ..., x_{2n-1}
        if (j < n) evenSum += f(a + (2.0 * j) * h);     // x_2, x_4, ..., x_{2n-2}
    }
    return h / 3.0 * (f(a) + 4.0 * oddSum + 2.0 * evenSum + f(b));
}
}