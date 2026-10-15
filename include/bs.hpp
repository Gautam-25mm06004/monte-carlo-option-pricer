#pragma once
// Black-Scholes closed forms and the Kemna-Vorst geometric-average Asian option (used as a control variate).
#include <cmath>

namespace mc {

inline double ncdf(double x) { return 0.5 * std::erfc(-x / std::sqrt(2.0)); }
inline double npdf(double x) { return std::exp(-0.5 * x * x) / std::sqrt(2.0 * M_PI); }

struct Market { double S0, K, r, sigma, T; };

inline double bs_call(const Market& m) {
    double sT = m.sigma * std::sqrt(m.T);
    double d1 = (std::log(m.S0 / m.K) + (m.r + 0.5 * m.sigma * m.sigma) * m.T) / sT, d2 = d1 - sT;
    return m.S0 * ncdf(d1) - m.K * std::exp(-m.r * m.T) * ncdf(d2);
}
inline double bs_put(const Market& m) { return bs_call(m) - m.S0 + m.K * std::exp(-m.r * m.T); }
inline double bs_delta(const Market& m) {
    double d1 = (std::log(m.S0 / m.K) + (m.r + 0.5 * m.sigma * m.sigma) * m.T) / (m.sigma * std::sqrt(m.T));
    return ncdf(d1);
}
inline double bs_vega(const Market& m) {
    double d1 = (std::log(m.S0 / m.K) + (m.r + 0.5 * m.sigma * m.sigma) * m.T) / (m.sigma * std::sqrt(m.T));
    return m.S0 * npdf(d1) * std::sqrt(m.T);
}
inline double bs_gamma(const Market& m) {
    double sT = m.sigma * std::sqrt(m.T);
    double d1 = (std::log(m.S0 / m.K) + (m.r + 0.5 * m.sigma * m.sigma) * m.T) / sT;
    return npdf(d1) / (m.S0 * sT);
}

// Geometric average over n equally spaced fixing dates t_i = iT/n (i = 1..n): log G is normal.
inline double geo_asian_call(const Market& m, int n) {
    double mu = std::log(m.S0) + (m.r - 0.5 * m.sigma * m.sigma) * m.T * (n + 1) / (2.0 * n);
    double var = m.sigma * m.sigma * m.T * (n + 1) * (2.0 * n + 1) / (6.0 * n * n);
    double sd = std::sqrt(var);
    double d1 = (mu - std::log(m.K) + var) / sd, d2 = d1 - sd;
    return std::exp(-m.r * m.T) * (std::exp(mu + 0.5 * var) * ncdf(d1) - m.K * ncdf(d2));
}

}  // namespace mc
