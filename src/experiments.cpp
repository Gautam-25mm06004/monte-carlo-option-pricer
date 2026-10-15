// Convergence and variance-reduction study; writes results/results.json and prints tables.
#include <chrono>
#include <cstdio>
#include "mc.hpp"

using namespace mc;
using clk = std::chrono::steady_clock;

int main() {
    Market m{100, 100, 0.05, 0.2, 1.0};
    const double exact = bs_call(m);
    FILE* js = std::fopen("results/results.json", "w");
    std::fprintf(js, "{\n \"market\": {\"S0\": 100, \"K\": 100, \"r\": 0.05, \"sigma\": 0.2, \"T\": 1},\n \"bs_call\": %.6f,\n", exact);
    std::printf("Black-Scholes call = %.6f\n\n", exact);

    // 1) Variance reduction at fixed N, plus efficiency (variance x time)
    const char* names[4] = {"plain", "antithetic", "control (S_T)", "antithetic + control"};
    Method ms[4] = {Method::Plain, Method::Antithetic, Method::Control, Method::AntitheticControl};
    const std::uint64_t N = 4000000;
    double base_var = 0, base_t = 0;
    std::printf("| method | price | std error | |error|/SE | variance reduction | time (s) | efficiency gain |\n|---|---|---|---|---|---|---|\n");
    std::fprintf(js, " \"european_1M\": [\n");
    for (int k = 0; k < 4; ++k) {
        auto t0 = clk::now();
        Estimate e = european_call(m, N, ms[k], 2, 7);
        double t = std::chrono::duration<double>(clk::now() - t0).count();
        double var = e.stderr_ * e.stderr_;
        if (k == 0) { base_var = var; base_t = t; }
        std::printf("| %s | %.4f | %.5f | %.2f | %.1fx | %.3f | %.1fx |\n", names[k], e.price, e.stderr_, std::fabs(e.price - exact) / e.stderr_, base_var / var, t, (base_var * base_t) / (var * t));
        std::fprintf(js, "  {\"method\": \"%s\", \"price\": %.6f, \"se\": %.6f, \"z\": %.3f, \"var_reduction\": %.2f, \"time_s\": %.4f, \"efficiency_gain\": %.2f}%s\n",
                     names[k], e.price, e.stderr_, (e.price - exact) / e.stderr_, base_var / var, t, (base_var * base_t) / (var * t), k < 3 ? "," : "");
    }
    std::fprintf(js, " ],\n");

    // 2) Convergence: error vs N should fall like 1/sqrt(N)
    std::printf("\n| paths | plain price | plain SE | anti+control price | anti+control SE |\n|---|---|---|---|---|\n");
    std::fprintf(js, " \"convergence\": [\n");
    std::uint64_t Ns[6] = {1000, 10000, 100000, 1000000, 10000000, 40000000};
    for (int i = 0; i < 6; ++i) {
        Estimate a = european_call(m, Ns[i], Method::Plain, 2, 11), b = european_call(m, Ns[i], Method::AntitheticControl, 2, 11);
        std::printf("| %llu | %.4f | %.4f | %.5f | %.5f |\n", (unsigned long long)Ns[i], a.price, a.stderr_, b.price, b.stderr_);
        std::fprintf(js, "  {\"paths\": %llu, \"plain\": %.6f, \"plain_se\": %.6f, \"ac\": %.6f, \"ac_se\": %.6f}%s\n", (unsigned long long)Ns[i], a.price, a.stderr_, b.price, b.stderr_, i < 5 ? "," : "");
    }
    std::fprintf(js, " ],\n");

    // 3) Greeks vs closed form
    Greeks g = european_greeks(m, 4000000, 2);
    std::printf("\n| Greek | Monte Carlo | std error | Black-Scholes |\n|---|---|---|---|\n");
    std::printf("| delta (pathwise) | %.5f | %.5f | %.5f |\n| vega (pathwise) | %.4f | %.4f | %.4f |\n| gamma (likelihood ratio) | %.5f | %.5f | %.5f |\n",
                g.delta, g.delta_se, bs_delta(m), g.vega, g.vega_se, bs_vega(m), g.gamma, g.gamma_se, bs_gamma(m));
    std::fprintf(js, " \"greeks\": {\"delta\": [%.6f, %.6f, %.6f], \"vega\": [%.6f, %.6f, %.6f], \"gamma\": [%.6f, %.6f, %.6f]},\n",
                 g.delta, g.delta_se, bs_delta(m), g.vega, g.vega_se, bs_vega(m), g.gamma, g.gamma_se, bs_gamma(m));

    // 4) Arithmetic Asian with geometric control variate
    std::printf("\n| Asian (52 weekly fixings) | price | std error | variance reduction |\n|---|---|---|---|\n");
    Estimate p0 = asian_call(m, 52, 1000000, false, 2), p1 = asian_call(m, 52, 1000000, true, 2);
    double vr = (p0.stderr_ * p0.stderr_) / (p1.stderr_ * p1.stderr_);
    std::printf("| plain | %.5f | %.5f | 1x |\n| geometric control variate | %.5f | %.5f | %.0fx |\n", p0.price, p0.stderr_, p1.price, p1.stderr_, vr);
    std::printf("(geometric Asian closed form = %.5f)\n", geo_asian_call(m, 52));
    std::fprintf(js, " \"asian\": {\"fixings\": 52, \"plain\": [%.6f, %.6f], \"cv\": [%.6f, %.6f], \"var_reduction\": %.1f, \"geo_closed_form\": %.6f},\n",
                 p0.price, p0.stderr_, p1.price, p1.stderr_, vr, geo_asian_call(m, 52));

    // 5) Thread scaling
    std::printf("\n| threads | time for 20M paths (s) | speed-up |\n|---|---|---|\n");
    std::fprintf(js, " \"scaling\": [\n");
    double t1 = 0;
    unsigned hw = std::thread::hardware_concurrency();
    for (unsigned th = 1; th <= hw; ++th) {
        auto t0 = clk::now();
        european_call(m, 20000000, Method::Plain, th, 5);
        double t = std::chrono::duration<double>(clk::now() - t0).count();
        if (th == 1) t1 = t;
        std::printf("| %u | %.3f | %.2fx |\n", th, t, t1 / t);
        std::fprintf(js, "  {\"threads\": %u, \"time_s\": %.4f, \"speedup\": %.3f}%s\n", th, t, t1 / t, th < hw ? "," : "");
    }
    std::fprintf(js, " ]\n}\n");
    std::fclose(js);
}
