# Multithreaded Monte Carlo Option Pricer with Variance Reduction (C++17)

Prices European and arithmetic-average Asian calls by Monte Carlo, measures how much each variance-reduction technique
helps, and checks Monte Carlo Greeks against Black–Scholes. Header-only, no dependencies beyond the standard library.

**Techniques:** antithetic variates · control variates (discounted terminal price for Europeans; the closed-form
geometric Asian of Kemna & Vorst (1990) for arithmetic Asians) · pathwise delta and vega · likelihood-ratio gamma ·
per-thread independent random streams.

## Results (S₀ = K = 100, r = 5%, σ = 20%, T = 1; exact Black–Scholes price 10.4506)

**European call, 4M payoff evaluations each**

| Method | Price | Std error | Variance reduction | Efficiency gain (variance × time) |
|---|---|---|---|---|
| Plain | 10.4620 | 0.00736 | 1× | 1× |
| Antithetic | 10.4636 | 0.00520 | 2.0× | 2.8× |
| Control variate (S_T) | 10.4567 | 0.00281 | 6.9× | 6.7× |
| **Antithetic + control** | **10.4525** | **0.00138** | **28.6×** | **44.4×** |

**Arithmetic Asian call, 52 weekly fixings, 1M paths:** the geometric-Asian control variate cuts the variance **1,315×**
(std error 0.0081 → 0.00022; price 5.854).

**Greeks (4M paths) vs closed form:** delta 0.63694 ± 0.00029 (BS 0.63683) · vega 37.54 ± 0.04 (BS 37.52) ·
gamma 0.01876 ± 0.00007 (BS 0.01876).

**Convergence:** the standard error falls as 1/√N, from 0.46 at 1,000 paths to 0.0023 at 40M (plain), and 0.097 → 0.00043
with antithetic + control.

**Threads:** 2 threads give a 1.32× speed-up on the 2-vCPU test VM (hyper-threaded cores share execution units).

Full tables: [`results/experiments.md`](results/experiments.md) · JSON: [`results/results.json`](results/results.json)

## Build and run
```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j
./build/test_pricer      # closed forms, parity, finite-difference Greeks, MC within 4 SE, reproducibility
./build/experiments      # writes results/results.json
```
