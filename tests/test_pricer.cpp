#include <cmath>
#include <cstdio>
#include "mc.hpp"
using namespace mc;
static int failures = 0;
#define CHECK(c) do { if (!(c)) { std::printf("FAIL line %d: %s\n", __LINE__, #c); ++failures; } } while (0)

int main() {
    Market m{100, 100, 0.05, 0.2, 1.0};
    CHECK(std::fabs(bs_call(m) - 10.450584) < 1e-5);                       // textbook value
    CHECK(std::fabs(bs_call(m) - bs_put(m) - (m.S0 - m.K * std::exp(-m.r * m.T))) < 1e-12);  // put-call parity
    // finite-difference checks of closed-form Greeks
    double h = 1e-3;
    Market up = m, dn = m; up.S0 += h; dn.S0 -= h;
    CHECK(std::fabs((bs_call(up) - bs_call(dn)) / (2 * h) - bs_delta(m)) < 1e-6);
    CHECK(std::fabs((bs_call(up) - 2 * bs_call(m) + bs_call(dn)) / (h * h) - bs_gamma(m)) < 1e-4);
    Market vu = m, vd = m; vu.sigma += h; vd.sigma -= h;
    CHECK(std::fabs((bs_call(vu) - bs_call(vd)) / (2 * h) - bs_vega(m)) < 1e-4);
    // geometric Asian with 1 fixing = European
    CHECK(std::fabs(geo_asian_call(m, 1) - bs_call(m)) < 1e-10);
    // every MC estimator within 4 standard errors of the exact price
    for (Method k : {Method::Plain, Method::Antithetic, Method::Control, Method::AntitheticControl}) {
        Estimate e = european_call(m, 400000, k, 2, 99);
        CHECK(std::fabs(e.price - bs_call(m)) < 4 * e.stderr_);
    }
    // variance reduction actually reduces the standard error
    CHECK(european_call(m, 400000, Method::AntitheticControl, 2).stderr_ < 0.5 * european_call(m, 400000, Method::Plain, 2).stderr_);
    // Greeks within 4 SE
    Greeks g = european_greeks(m, 400000, 2);
    CHECK(std::fabs(g.delta - bs_delta(m)) < 4 * g.delta_se);
    CHECK(std::fabs(g.vega - bs_vega(m)) < 4 * g.vega_se);
    CHECK(std::fabs(g.gamma - bs_gamma(m)) < 4 * g.gamma_se);
    // Monte Carlo geometric Asian matches its closed form (control-variate expectation is correct)
    Estimate a = asian_call(m, 12, 400000, false, 2), b = asian_call(m, 12, 400000, true, 2);
    CHECK(std::fabs(a.price - b.price) < 4 * a.stderr_);
    CHECK(b.stderr_ < 0.2 * a.stderr_);
    // same seed + same thread count -> identical result
    CHECK(european_call(m, 10000, Method::Plain, 2, 5).price == european_call(m, 10000, Method::Plain, 2, 5).price);
    if (failures) { std::printf("%d FAILURES\n", failures); return 1; }
    std::printf("ALL TESTS PASSED\n");
}
