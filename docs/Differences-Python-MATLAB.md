# Differences between the Python and MATLAB implementations

The two implementations follow the same method and the same sequence of operations. This page lists where they differ, as found by reading the code of both.
They are **not** bit-for-bit equivalent: different random number generators, solvers and tolerances give different chains for the same seed.

## Summary

| Aspect | Python | MATLAB |
|---|---|---|
| Entry point | `python/examples/CHRJ_3comp/example_hefa.py` | `matlab/examples/CHRJ_3comp/example_HEFA.m` |
| Sampler arguments | keyword arguments with defaults, more settings exposed | 29 positional arguments; several settings fixed inside the file |
| Seed | argument `seed` (default 1234) | fixed `baseSeed = 1234`, not an argument |
| Random numbers | `RandomState` (MT19937) with a ziggurat built to follow MATLAB's `randn` stream | `rng(1234, 'twister')`, workers seeded `1234 + spmdIndex` if a pool exists |
| Gibbs move | enumeration of all candidates (`gibbs_tries=None`) or multiple-try | multiple-try only; `gibbs_tries` required |
| `scaling_factor_nc` | scalar or array | scalar |
| Log file | argument `log_file` (default `MCMC.txt`; the example uses `demc_log.txt`) | always `MCMC.txt` in the current folder |
| Progress display | text line every 100 iterations | waitbar with *Cancel* button (desktop) or text line every 100 iterations |
| Parallelism | `multiprocessing` pool with shared memory; optional numba | `parfor` (Parallel Computing Toolbox) |
| Distillation root finder | Brent's method (tolerance 0.01 K) or compiled numba version | `fzero` (default tolerance) |
| Boiling temperature | Brent's method in `[200, 1200]` K | `fsolve` from 500 K (Optimization Toolbox) |
| Where the data are found | `data_dir` argument (path) | MATLAB path (`addpath`); file names without extension |
| Species database | `database/` located from the source file; editable install | folders of `database/` added to the path by the example |
| `deltaT_dist` and the distillation cache | recomputes its own curves | shares the curves with `distillation` through a `persistent` cache, when the batch is the same |
| Mono-/bis-substituted isoparaffins | removed | not present |
| Plots | `visualize_results` (also `AcceptanceRate` and `Rhat`) | inside `props_pushforward` (no `AcceptanceRate` and `Rhat`) |
| Pushforward and plots | two functions: `props_pushforward`, `visualize_results` | one function: `props_pushforward` |
| Sobol' analysis | `seed` argument (default 0) | local stream with seed 0; the global generator is untouched |
| Sobol' samples | evaluates the model on the original (unrounded) distinct samples | evaluates on the distinct samples rounded to 6 decimals |
| Output of the example | `MAP.txt`, `demc_log.txt`, `Figures/` | `MAP.txt`, `MCMC.txt`, `Figures/`, `CHRJ_3comp.mat` |

## What is the same

- the parametrisation `(x, n_C, eta)`, the decoding and the folding at the boundaries;
- the prior, the isomer correction and the Gaussian likelihood, with the same mixing rules and units;
- the distillation curve algorithm in 100 flash steps and the cheap variant;
- the two-stage delayed acceptance, the three-phase burn-in with tempering and outlier resets, the production ladder, the snooker and Gibbs moves, the swap rules;
- the stopping rule (R-hat only) and the effective sample size as a diagnostic;
- the MAP (highest log-posterior sample), the 25-point temperature grids and the percentile bands, the grouped Sobol' indices;
- the family keys (`alkylnaphthalenes`), the detection of the property from the name of the data file (first match, error if none) and the treatment of flash points at or below 0 degF (`-Inf` likelihood, `NaN` in `MAP.txt`).

## Practical consequences

- Do not expect identical numbers between the two versions. Compare the posterior summaries (MAP, credible bands), not the chains.
- A setting that you can change from the call in Python (seed, `rhat_spacing`, `outlier_every`, `burnin_cap`, ...) may require an edit of
  `DifferentialEvolutionMarkovChain.m` in MATLAB.
- Check the known limitations of each implementation ([Python](python/Known-limitations.md), [MATLAB](matlab/Known-limitations.md)) before moving a setup from one to the other.
