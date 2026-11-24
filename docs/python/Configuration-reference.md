# Configuration reference

Every parameter of the CHRJ example (`python/examples/CHRJ_3comp/example_hefa.py`) and of
`run_demc` (`mcmc/demc.py`). The "Example" column gives the value used in the example. See
[Sampler](../Sampler.md) for the meaning of each parameter in the algorithm.

## Real fuel and data

| Name | Meaning | Example |
|---|---|---|
| `DATASET` | properties to fit; each must be a file `<name>.csv` in `DATA_DIR` | `rho, distillation, molWeight, HC, DCN, LHV, flash` |
| `DATA_DIR` | folder with the CSV files, relative to the script | `../../../exp/CHRJ_POSF6152` |
| `FUEL_NAME` | label in figures and `MAP.txt` | `CHRJ POSF-6152` |
| `T_BUBBLE` | bubble temperature of the real fuel [K]; `1e18` if unknown | 445.00 |

## Surrogate definition

| Name | Meaning | Example |
|---|---|---|
| `NUM_COMPONENTS` | number of components `N_s` | 3 |
| `FAMILIES` | family of each component (see [Input data](../Input-data.md#species-database-database)) | `nparaffins, isoparaffins, isoparaffins` |
| `N_RANGES` | carbon numbers of each component | C8-C15, C9-C12, C13-C16 |
| `LOWER_ETA`, `UPPER_ETA` | bounds on `eta_B_star_norm`; act only in the prior | 0, 1 |
| `LOWER_MOL_FRAC`, `UPPER_MOL_FRAC` | bounds on each molar fraction, including the last | 0.1, 0.75 |
| `ALPHA_DIRICHLET` | Dirichlet concentration of the prior | `ones(N_s)` |

`N_RANGES` is used to load the species (`load_hydrocarbons`) and in the discrete-uniform prior. The sampler takes the
carbon numbers it can reach from the loaded species lists, not from `N_RANGES` directly.
The `n_ranges`, `lower_bound_eta` and `upper_bound_eta` arguments of `run_demc` are accepted for
compatibility with the MATLAB signature but the standard layout does not use them; the bounds on `η*_B`
take effect only through the prior in `build_pdfs`.

## Chains

| Name (`run_demc` argument) | Meaning | Default | Example |
|---|---|---|---|
| `N_CHAINS` (`N_chains`) | number of chains; at least 8 | required | `2*(3*N_s-1)` = 16 |
| `MAX_ITERATIONS` (`max_iterations`) | hard stop, per chain | required | 20000 |
| `SEED` (`seed`) | random seed | 1234 | 42 |
| `x_init` | normalised initial population, shape `(N_chains, 3*N_s-1)`; bypasses the random initialisation | `None` | not used |
| `layout` | `ParameterLayout`; `None` builds `SurrogateLayout` | `None` | not used |
| `log_file` | text log of the run; `None` disables it | `"MCMC.txt"` | `demc_log.txt` |
| `info` | dict filled with run information | `None` | `{}` |

## Moves

| Name | Meaning | Default | Example |
|---|---|---|---|
| `SCALING_X`, `SCALING_NC`, `SCALING_ETA` | jump scaling for the three blocks; `scaling_factor_nc` may be an array | required | 0.1, 1.0, 1.0 |
| `NOISE_X`, `NOISE_NC`, `NOISE_ETA` | std of the Gaussian noise added to the jumps | required | 1e-6 |
| `N_PAIRS` (`n_pairs`) | maximum number of chain pairs in the DE jump; `D` is drawn in `1..n_pairs` | 1 | 1 |
| `P_GIBBS` (`p_gibbs`) | probability of the Gibbs-type move per chain and iteration | 0.0 | 0.3 |
| `P_SNOOKER` (`p_snooker`) | probability of the snooker move once the fractions block is selected | 0.0 | 0.1 |
| `GIBBS_TRIES` (`gibbs_tries`) | number of tries of the multiple-try Gibbs move; `None` enumerates all candidates | `None` | 32 |

## Burn-in and parallel tempering

| Name | Meaning | Default | Example |
|---|---|---|---|
| `T_BURNIN` (`t_burnin`) | length of the tempered phase A; the real burn-in is at least `2*t_burnin` | required | 500 |
| `BETA_MIN` (`beta_min`) | lowest inverse temperature of the burn-in ladder | 1e-4 | 1e-4 |
| `T_LADDER` (`T_ladder`) | `"geometric"` or `"linear"` | `"geometric"` | `geometric` |
| `SWAP_FREQ` (`swap_freq`) | iterations between swaps in phase A | 20 | 20 |
| `outlier_every` | iterations between outlier checks in phase B | `max(50, t_burnin//4)` | not set |
| `burnin_cap` | burn-in is forced to end at `burnin_cap * max_iterations` | 0.75 | not set |
| `N_COLD` (`n_cold`) | number of cold chains in production; at least 4, and `N_chains - n_cold >= 4` | `N_chains // 2` | `N_CHAINS // 2` |
| `beta_min_sampling` | lowest `beta` of the production ladder | 0.05 | not set |
| `swap_every_sampling` | iterations between swaps in production | 1 | not set |

## Convergence

| Name | Meaning | Default | Example |
|---|---|---|---|
| `RHAT_STOP` (`rhat_stop`) | R-hat threshold | 1.1 | 1.1 |
| `rhat_checks` | consecutive checks below the threshold | 3 | not set |
| `rhat_spacing` | iterations between checks; **must be a positive multiple of 100** (`ValueError` otherwise) | 500 | not set |
| `RHAT_MIN_SAMPLES` (`rhat_min_samples`) | minimum post-burn-in iterations | 2000 | 2000 |
| `ESS_MIN` (`ess_min`) | **diagnostic only**: a warning is written in the log if the minimum bulk ESS of the molar fractions is below it; it does not affect the stopping | 400 | 400 |

## Extension points of `run_demc`

These arguments serve layouts other than the standard one (for example the hierarchical model of the `bayesaf-hierarchical` repository). They have no effect with the standard layout and are not exposed in the example.

| Name | Meaning | Default |
|---|---|---|
| `layout` | `ParameterLayout` describing the state vector (units, blocks, posterior arguments) | standard `SurrogateLayout` |
| `scaling_factor_cont`, `noise_cont` | jump scale and noise of the `"continuous"` blocks of a layout | 1.0, 1e-6 |
| `crossover_every` | iterations between exchanges of one block of coordinates between two chains of the same temperature group (layouts with `crossover_blocks`) | `None` (off) |
| `block_proposals`, `p_block_proposal`, `block_proposal_rounds` | independence proposals of whole blocks drawn from an object with `blocks`, `sample` and `logq`; probability per chain and round; rounds per iteration | `None`, 0.0, 1 |
| `aux_moves` | object with `step(X, phys, p_X, beta, rng, adapt)`, called once per iteration after the other moves | `None` |
| `de_skip_units` | units whose discrete blocks the DE move never draws | `None` |
| `checkpoint_every`, `checkpoint_path` | write a pickle of the cold chains every that many iterations | `None` (off) |
| `burnin_resets` | if `False`, the chains below the others are not reset during the burn-in, which ends at `2*t_burnin` | `True` |

## Post-processing

| Name | Meaning | Example |
|---|---|---|
| `CONFIDENCE_WIDTH` | width of the credible interval in the pushforward and figures (0.95 gives the 2.5-97.5% band) | 0.95 |
| `BAND_PERCENTILES` | colour the band by percentiles | `True` |

## Choosing values

- **Chains.** The example uses `2*(3*N_s-1)` chains, half cold and half hot after the burn-in. For `N_s = 3` this is 16.
- **Iterations.** The run stops early when the convergence criterion is met, so `MAX_ITERATIONS` is a ceiling. Increase it first if
  the run reports no convergence.
- **Cost.** Per iteration, the cost is dominated by distillation curves of the proposals that pass stage 1. A smaller
  `P_GIBBS` with `GIBBS_TRIES = None` can be costlier than a multiple-try move with few tries, because it enumerates all candidates.
- **Diagnostics to check** after a run: `info["swap_rates"]` (not near 0 or 1), `info["n_resets"]` (outlier resets in
  the burn-in), the acceptance rate `AR` and the R-hat trace.
