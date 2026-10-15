Black-Scholes call = 10.450584

| method | price | std error | |error|/SE | variance reduction | time (s) | efficiency gain |
|---|---|---|---|---|---|---|
| plain | 10.4620 | 0.00736 | 1.54 | 1.0x | 0.110 | 1.0x |
| antithetic | 10.4636 | 0.00520 | 2.49 | 2.0x | 0.077 | 2.8x |
| control (S_T) | 10.4567 | 0.00281 | 2.19 | 6.9x | 0.113 | 6.7x |
| antithetic + control | 10.4525 | 0.00138 | 1.39 | 28.6x | 0.071 | 44.4x |

| paths | plain price | plain SE | anti+control price | anti+control SE |
|---|---|---|---|---|
| 1000 | 10.2009 | 0.4615 | 10.50222 | 0.09687 |
| 10000 | 10.6430 | 0.1513 | 10.46339 | 0.02802 |
| 100000 | 10.4660 | 0.0466 | 10.44379 | 0.00860 |
| 1000000 | 10.4620 | 0.0147 | 10.44873 | 0.00275 |
| 10000000 | 10.4561 | 0.0047 | 10.45099 | 0.00087 |
| 40000000 | 10.4544 | 0.0023 | 10.45120 | 0.00043 |

| Greek | Monte Carlo | std error | Black-Scholes |
|---|---|---|---|
| delta (pathwise) | 0.63694 | 0.00029 | 0.63683 |
| vega (pathwise) | 37.5418 | 0.0378 | 37.5240 |
| gamma (likelihood ratio) | 0.01876 | 0.00007 | 0.01876 |

| Asian (52 weekly fixings) | price | std error | variance reduction |
|---|---|---|---|
| plain | 5.87649 | 0.00813 | 1x |
| geometric control variate | 5.85399 | 0.00022 | 1315x |
(geometric Asian closed form = 5.63743)

| threads | time for 20M paths (s) | speed-up |
|---|---|---|
| 1 | 0.745 | 1.00x |
| 2 | 0.563 | 1.32x |
