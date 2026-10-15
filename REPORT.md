# Variance Reduction for Monte Carlo Option Pricing: A Quantitative Comparison

**Gautam Jain** — IIT Bhubaneswar

## Abstract
We implement a multithreaded Monte Carlo pricer in C++17 and measure the variance and efficiency gains of antithetic and
control variates for a European call and an arithmetic Asian call under Black–Scholes dynamics. Combining antithetic
sampling with a terminal-price control variate reduces variance 28.6× for the European call (44.4× after accounting for
run time); a geometric-Asian control variate reduces the Asian option's variance 1,315×. Pathwise and likelihood-ratio
Greek estimators agree with closed forms within their standard errors.

## 1. Methods
Under the risk-neutral measure, `S_T = S₀ exp((r − σ²/2)T + σ√T Z)`.
- **Antithetic**: average the payoffs for Z and −Z; the pair has negative correlation for monotone payoffs.
- **Control variate**: `Y − b(X − E[X])` with `b = Cov(X,Y)/Var(X)` estimated from the same paths. For the European call,
  `X = e^{−rT} S_T` with `E[X] = S₀`. For the arithmetic Asian, X is the discounted geometric-average payoff, whose price is
  known in closed form because the log of a geometric average of lognormals is normal (Kemna & Vorst 1990).
- **Greeks**: pathwise delta `e^{−rT} 1{S_T > K} S_T/S₀` and vega `e^{−rT} 1{S_T > K} S_T(−σT + √T Z)`; likelihood-ratio
  gamma with score weight `(Z² − 1)/(S₀²σ²T) − Z/(S₀²σ√T)` (the payoff's second derivative does not exist pathwise).
- **Parallelism**: each thread owns an `mt19937_64` stream seeded from (seed, thread id) and its own accumulators, merged at
  the end; no shared state in the hot loop.

## 2. Verification
Unit tests check the closed form against the textbook value 10.450584, put–call parity, closed-form Greeks against finite
differences, that the geometric Asian with one fixing equals Black–Scholes, that every estimator lies within 4 standard
errors of the exact price, and reproducibility for a fixed seed and thread count.

## 3. Results
See README tables. Two observations:
1. **The best single technique depends on the payoff.** For the at-the-money European call, the terminal-price control
   variate (correlation ≈ 0.93 with the payoff) beats antithetic sampling; combining them works best because antithetic
   pairs already remove the linear part of the noise that the control variate does not capture.
2. **A closely related closed-form instrument is an extremely strong control.** Arithmetic and geometric averages of the same
   path are almost perfectly correlated, which explains the 1,315× reduction.

Because all European methods reuse the same random stream, their errors are correlated (all estimates sit slightly above
the exact price); the z-scores range from 1.4 to 2.5.

## 4. Limitations
Black–Scholes dynamics only; pseudo-random rather than quasi-random (Sobol) numbers; no path-dependent Greeks for Asians;
thread scaling measured on a 2-vCPU virtual machine.

## References
1. P. Boyle. Options: a Monte Carlo approach. *Journal of Financial Economics* 4(3), 323–338, 1977.
2. P. Glasserman. *Monte Carlo Methods in Financial Engineering*. Springer, 2003.
3. A. Kemna, A. Vorst. A pricing method for options based on average asset values. *Journal of Banking & Finance* 14(1), 113–129, 1990.
4. M. Broadie, P. Glasserman. Estimating security price derivatives using simulation. *Management Science* 42(2), 269–285, 1996.
