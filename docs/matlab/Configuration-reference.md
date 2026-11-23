# Configuration reference (MATLAB)

Every setting of `example_HEFA.m` and every argument of `DifferentialEvolutionMarkovChain`. For the meaning of each parameter in the algorithm, see [Sampler](../Sampler.md).

## Real fuel and data

| Variable | Meaning | Example |
|---|---|---|
| `Dataset` | cell array of file names (no extension), one per property; the files must be on the MATLAB path | `{'rho','distillation','molWeight','HC','DCN','LHV','flash'}` |
| `fuel_name` | label in the figures and in `MAP.txt` | `'CHRJ POSF-6152'` |
| `Tbubble_realFuel` | bubble temperature of the real fuel [K]; `1e+18` if unknown | 445.00 |

## Surrogate definition

| Variable | Meaning | Example |
|---|---|---|
| `numComponents` | number of components | 3 |
| `families` | cell array with the family of each component (see [Input data](../Input-data.md#species-database-database)) | `{'nparaffins','isoparaffins','isoparaffins'}` |
| `n_ranges` | cell array with the carbon numbers of each component | `{8:15, 9:12, 13:16}` |
| `LowerBound_eta_B_star`, `UpperBound_eta_B_star` | bounds on `eta_B_star_norm`, used by the prior | `[0 0 0]`, `[1 1 1]` |
| `LowerBound_molFrac`, `UpperBound_molFrac` | bounds on each molar fraction, including the last | `[0.1 0.1 0.1]`, `[0.75 0.75 0.75]` |
| `alpha` | Dirichlet concentration | `ones(1, numComponents)` |
| `classes` | `classes{i} = Hydrocarbons(families{i}, n_ranges{i})` | |

## Arguments of `DifferentialEvolutionMarkovChain`

The arguments are **positional**; the order is the one in the table.

| # | Argument | Meaning | Example |
|---|---|---|---|
| 1 | `families` | families of the components | |
| 2 | `Posterior_PDF` | handle of the full posterior (`pdf_build_speed`) | `posterior` |
| 3 | `Posterior_PDF_cheap` | handle of the cheap posterior (`pdf_build_speed_cheap`) | `posterior_cheap` |
| 4 | `classes` | species of each component | |
| 5, 6 | `LowerBound_x`, `UpperBound_x` | bounds of the molar fractions | |
| 7 | `n_ranges` | carbon numbers of each component | |
| 8, 9 | `LowerBound_eta_B_star`, `UpperBound_eta_B_star` | bounds of `eta` (not used by the sampler; they act through the prior in the posterior handles) | |
| 10 | `maxIterations` | hard stop per chain | 20000 |
| 11 | `t_burnin` | length of the tempered phase of the burn-in; the real burn-in is at least `2*t_burnin` | 500 |
| 12 | `scaling_factor_X` | jump scaling, molar fractions | 0.1 |
| 13 | `scaling_factor_nc` | jump scaling, carbon numbers (a scalar) | 1.0 |
| 14 | `scaling_factor_eta` | jump scaling, topochemical indices | 1.0 |
| 15 | `noise_X` | noise of the jumps, molar fractions | 1e-06 |
| 16 | `noise_nc` | noise, carbon numbers | 1e-06 |
| 17 | `noise_eta` | noise, topochemical indices | 1e-06 |
| 18 | `N_chains` | number of chains; at least 8 | `2*(3*numComponents-1)` = 16 |
| 19 | `beta_min` | lowest inverse temperature in the burn-in ladder | 0.0001 |
| 20 | `T_ladder` | `'linear'` or `'geometric'` | `'geometric'` |
| 21 | `swap_freq` | iterations between swaps in the burn-in | 20 |
| 22 | `n_pairs` | maximum number of chain pairs in the DE jump | 1 |
| 23 | `p_gibbs` | probability of the Gibbs move | 0.3 |
| 24 | `p_snooker` | probability of the snooker move once the fractions block is chosen | 0.1 |
| 25 | `gibbs_tries` | tries of the multiple-try Gibbs move (required, no `None` option) | 32 |
| 26 | `n_cold` | chains at beta = 1 after the burn-in; at least 4, and `N_chains - n_cold >= 4` | `N_chains/2` |
| 27 | `rhat_stop` | R-hat threshold | 1.1 |
| 28 | `rhat_min_samples` | minimum iterations after the burn-in | 2000 |
| 29 | `ess_min` | **diagnostic only**: warning if the minimum bulk ESS of the molar fractions is below it | 400 |

## Settings fixed inside the sampler

These are variables at the beginning of `DifferentialEvolutionMarkovChain.m`, not arguments:

| Variable | Value |
|---|---|
| `baseSeed` | 1234 |
| `beta_min_sampling` | 0.05 |
| `swap_every_sampling` | 1 |
| `outlier_every` | `max(50, floor(t_burnin/4))` |
| `burnin_cap` | 0.75 |
| `rhat_checks` | 3 |
| `rhat_spacing` | 500 |

Edit the file to change them. Keep `rhat_spacing` a multiple of 100.

## Post-processing

| Variable | Meaning | Example |
|---|---|---|
| `band_percentiles` | `'True'` colours the credible band with the percentiles, `'False'` does not | `'True'` |
| `confidence_width` | width of the credible interval (0.95 gives the 2.5-97.5 % band) | 0.95 |

## Choosing values

- `maxIterations` is a ceiling: the run stops earlier if the R-hat criterion is met. Increase it first if it does not converge.
- With the parallel pool, the cost per iteration is dominated by the distillation curves of the proposals that pass the first stage.
- After a run, look at `info.swap_rates`, `info.n_resets`, `AR(end,:)`, `R_hat` and the log `MCMC.txt`.
