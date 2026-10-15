#pragma once
// Multithreaded Monte Carlo pricer with antithetic variates, control variates and pathwise / likelihood-ratio Greeks.
#include <cmath>
#include <cstdint>
#include <random>
#include <thread>
#include <vector>
#include "bs.hpp"

namespace mc {

struct Estimate { double price, stderr_; std::uint64_t paths; };

// Running sums for mean / variance / covariance of (payoff Y, control X).
struct Acc {
    double sy = 0, syy = 0, sx = 0, sxx = 0, sxy = 0; std::uint64_t n = 0;
    void add(double y, double x = 0) { sy += y; syy += y * y; sx += x; sxx += x * x; sxy += x * y; ++n; }
    void merge(const Acc& o) { sy += o.sy; syy += o.syy; sx += o.sx; sxx += o.sxx; sxy += o.sxy; n += o.n; }
};

inline Estimate plain(const Acc& a) {
    double m = a.sy / a.n, v = a.syy / a.n - m * m;
    return {m, std::sqrt(v * a.n / (a.n - 1) / a.n), a.n};
}
// Control variate: Y_cv = Y - b (X - E[X]) with b = Cov(X,Y)/Var(X) estimated from the same paths.
inline Estimate controlled(const Acc& a, double EX) {
    double my = a.sy / a.n, mx = a.sx / a.n;
    double vx = a.sxx / a.n - mx * mx, vy = a.syy / a.n - my * my, cxy = a.sxy / a.n - mx * my;
    double b = cxy / vx;
    double v = vy - cxy * cxy / vx;
    return {my - b * (mx - EX), std::sqrt(v / (a.n - 1)), a.n};
}

template <class PathFn>
Acc run_parallel(std::uint64_t n_paths, unsigned threads, std::uint64_t seed, PathFn fn) {
    std::vector<Acc> parts(threads);
    std::vector<std::thread> pool;
    for (unsigned t = 0; t < threads; ++t) {
        pool.emplace_back([&, t] {
            std::mt19937_64 rng(seed * 1000003ull + t);           // independent stream per thread
            std::normal_distribution<double> N01;
            std::uint64_t mine = n_paths / threads + (t < n_paths % threads ? 1 : 0);
            for (std::uint64_t i = 0; i < mine; ++i) fn(rng, N01, parts[t]);
        });
    }
    for (auto& th : pool) th.join();
    Acc all;
    for (auto& p : parts) all.merge(p);
    return all;
}

enum class Method { Plain, Antithetic, Control, AntitheticControl };

// European call. Control variate: discounted S_T, whose expectation is S0.
inline Estimate european_call(const Market& m, std::uint64_t n, Method meth, unsigned threads = 1, std::uint64_t seed = 1) {
    const double drift = (m.r - 0.5 * m.sigma * m.sigma) * m.T, vol = m.sigma * std::sqrt(m.T), df = std::exp(-m.r * m.T);
    bool anti = meth == Method::Antithetic || meth == Method::AntitheticControl;
    bool cv = meth == Method::Control || meth == Method::AntitheticControl;
    std::uint64_t draws = anti ? n / 2 : n;   // same number of payoff evaluations for a fair comparison
    Acc a = run_parallel(draws, threads, seed, [&](std::mt19937_64& g, std::normal_distribution<double>& N, Acc& acc) {
        double z = N(g);
        double s1 = m.S0 * std::exp(drift + vol * z);
        double y = df * std::max(s1 - m.K, 0.0), x = df * s1;
        if (anti) {
            double s2 = m.S0 * std::exp(drift - vol * z);
            y = 0.5 * (y + df * std::max(s2 - m.K, 0.0));
            x = 0.5 * (x + df * s2);
        }
        acc.add(y, x);
    });
    return cv ? controlled(a, m.S0) : plain(a);
}

// Greeks for the European call: pathwise delta & vega, likelihood-ratio gamma.
struct Greeks { double delta, delta_se, vega, vega_se, gamma, gamma_se; };
inline Greeks european_greeks(const Market& m, std::uint64_t n, unsigned threads = 1, std::uint64_t seed = 2) {
    const double drift = (m.r - 0.5 * m.sigma * m.sigma) * m.T, vol = m.sigma * std::sqrt(m.T), df = std::exp(-m.r * m.T);
    std::vector<Acc> parts(3);
    Acc d = run_parallel(n, threads, seed, [&](auto& g, auto& N, Acc& acc) {
        double z = N(g), sT = m.S0 * std::exp(drift + vol * z);
        acc.add(sT > m.K ? df * sT / m.S0 : 0.0);
    });
    Acc v = run_parallel(n, threads, seed, [&](auto& g, auto& N, Acc& acc) {
        double z = N(g), sT = m.S0 * std::exp(drift + vol * z);
        acc.add(sT > m.K ? df * sT * (std::log(sT / m.S0) - (m.r + 0.5 * m.sigma * m.sigma) * m.T) / m.sigma : 0.0);
    });
    Acc ga = run_parallel(n, threads, seed, [&](auto& g, auto& N, Acc& acc) {
        double z = N(g), sT = m.S0 * std::exp(drift + vol * z);
        double w = ((z * z - 1.0) / (m.S0 * m.S0 * m.sigma * m.sigma * m.T)) - z / (m.S0 * m.S0 * vol);   // LR weight for d2/dS0^2
        acc.add(df * std::max(sT - m.K, 0.0) * w);
    });
    Estimate D = plain(d), V = plain(v), G = plain(ga);
    return {D.price, D.stderr_, V.price, V.stderr_, G.price, G.stderr_};
}

// Arithmetic-average Asian call with n_fix monitoring dates; optional geometric-Asian control variate.
inline Estimate asian_call(const Market& m, int n_fix, std::uint64_t n, bool use_cv, unsigned threads = 1, std::uint64_t seed = 3) {
    const double dt = m.T / n_fix, drift = (m.r - 0.5 * m.sigma * m.sigma) * dt, vol = m.sigma * std::sqrt(dt), df = std::exp(-m.r * m.T);
    Acc a = run_parallel(n, threads, seed, [&](auto& g, auto& N, Acc& acc) {
        double logS = std::log(m.S0), sum = 0, sumlog = 0;
        for (int i = 0; i < n_fix; ++i) { logS += drift + vol * N(g); sum += std::exp(logS); sumlog += logS; }
        double arith = sum / n_fix, geo = std::exp(sumlog / n_fix);
        acc.add(df * std::max(arith - m.K, 0.0), df * std::max(geo - m.K, 0.0));
    });
    return use_cv ? controlled(a, geo_asian_call(m, n_fix)) : plain(a);
}

}  // namespace mc
