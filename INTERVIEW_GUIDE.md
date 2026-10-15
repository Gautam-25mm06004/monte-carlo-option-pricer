# Interview guide — Monte Carlo option pricer

## 30-second version
"I wrote a multithreaded Monte Carlo pricer in C++ and measured variance reduction properly — at a fixed number of payoff
evaluations and also per unit of compute time. Antithetic plus a terminal-price control variate gives 28× lower variance on
a European call; for an arithmetic Asian, using the closed-form geometric Asian as a control gives over 1,000×. I also
compute delta and vega pathwise and gamma by likelihood ratio, and they match Black–Scholes within their standard errors."

## Explain these
- **Standard error** = sample std / √N → to halve the error you need 4× the paths. Variance reduction is the cheap way.
- **Antithetic**: use Z and −Z. Works when the payoff is monotone in Z (calls are).
- **Control variate**: subtract something with known mean that is correlated with the payoff; optimal b = Cov/Var.
  Variance falls by a factor 1/(1 − ρ²).
- **Why geometric Asian works so well**: geometric and arithmetic averages of the same path are almost the same number;
  the geometric one has a closed-form price.
- **Pathwise vs likelihood ratio**: pathwise differentiates the payoff (needs it to be continuous — fine for delta/vega
  of a call); LR differentiates the density (works for gamma or digitals, but noisier).
- **Efficiency** = 1/(variance × time). A method that reduces variance but costs more time can lose.

## Likely questions
1. *How do you parallelise MC correctly?* Independent RNG streams per thread (different seeds or skip-ahead), local
   accumulators, merge at the end; never share an RNG.
2. *Why only 1.32× with 2 threads?* The VM's 2 vCPUs are hyper-threads of one core sharing floating-point units.
3. *What would you add?* Sobol quasi-random numbers with Brownian bridge, importance sampling for deep out-of-the-money
   options, Heston or local-volatility paths, GPU.
4. *How do you know your price is right?* Compare to closed form within a few standard errors, put–call parity, and limits
   (one-fixing Asian = European).
