# Post-processing and outputs

## Return values of `run_demc`

```python
(chain, posterior_pdf, chain_reshaped, posterior_reshaped,
 AR, R_hat, t_convergence) = run_demc(...)
```

Seven values are returned (the type annotation in the source mentions eight).

| Name | Shape | Content |
|---|---|---|
| `chain` | `(t_convergence, 3*Nc-1, n_cold)` | **physical** parameters of the cold chains, columns as in the layout. Rows after the iteration at which the run stopped are not included |
| `posterior_pdf` | `(t_convergence, n_cold)` | log-posterior of the corrected posterior (with the isomer correction) |
| `chain_reshaped` | `(N_post * n_cold, 3*Nc-1)` | post-burn-in samples of all cold chains stacked |
| `posterior_reshaped` | `(N_post * n_cold,)` | corresponding log-posterior |
| `AR` | `(t_convergence, n_cold)` | cumulative acceptance rate, in percent |
| `R_hat` | `(n_checks, 3*Nc-1)` | R-hat history, one row every 2 iterations after the burn-in; `NaN` for constant-`eta` columns |
| `t_convergence` | int | iteration of the stop, or `max_iterations` if not converged |

The `info` dictionary has:

| Key | Content |
|---|---|
| `t_burnin` | iteration at which the burn-in actually ended |
| `n_resets` | total outlier chains reset |
| `swap_rates` | acceptance of the swaps at each production ladder level |
| `betas` | the production inverse temperatures |
| `ESS` | dictionary `{column index: bulk ESS}` of the molar fractions, on the returned samples |
| `ESS_min_at_stop` | minimum of the previous values (diagnostic) |

If the run does **not** converge and the burn-in did not finish, `t_burn = t_convergence - 1` and the
post-burn-in sample contains only the last iteration of each cold chain. Check `info["t_burnin"]` before using
`chain_reshaped`.

## MAP surrogate

`find_map(posterior_reshaped, chain_reshaped, Nc)` returns the **sample with the highest log-posterior**
(`argmax`) among the post-burn-in samples:

- `x_MAP`: molar fractions (the last one is `1 - sum`),
- `nc_MAP`: carbon numbers,
- `eta_B_star_MAP`: normalised indices.

It is the posterior mode as visited by the sampler. It is not a marginal mode and it is not a
posterior mean.

`write_map(families, Nc, x_MAP, nc_MAP, eta_B_star_MAP, fuel_name, output_file="MAP.txt", data_dir=None)`
writes a text file and returns the species names. It finds each species with the same lookup used in the
likelihood and computes the blend properties of the MAP mixture. `MAP.txt` contains:

- number of components, families and species names;
- molecular weight, molecular formula (`C nC_mix H nH_mix`), hydrogen-to-carbon ratio;
- derived cetane number, flash point (NaN with a warning if a component has `T_f` below 0 degF), freezing point, lower heating value;
- molar, mass and volume fractions;
- carbon numbers and topochemical indices;
- pseudo-critical temperature and pressure (the latter in bar), from the Lee-Kesler mixing rules
  (`Vc_ij = (Vc_i^(1/3) + Vc_j^(1/3))^3 / 8`, `Tc_ij = sqrt(Tc_i Tc_j)`, mole-averaged acentric factor,
  `Pc_m = (0.2905 - 0.085 omega_m) R Tc_m / Vc_m`).

`data_dir` can point to another copy of the database for the name lookup.

## Pushforward distributions

```python
pf = props_pushforward(full_data, families, classes, chain_reshaped, Nc, variable_names,
                       x_MAP, nc_MAP, eta_B_star_MAP, min_T_array, max_T_array,
                       pressure_distillation, confidence_width=0.9)
```

For each property in `variable_names`, the posterior samples are propagated through the model:

- temperature-dependent properties are evaluated on a grid of 25 temperatures between `min_T_array` and `max_T_array`;
- the distillation curve is evaluated on the volume grid `[0.1, 1, 2, ..., 99, 99.9]` % and interpolated;
- `deltaT_dist` returns `T50 - T10` and `T90 - T10`;
- lumped properties (`molWeight`, `HC`, `DCN`, `flash`, `freezing`, `LHV`) are evaluated once per sample.

To save time, identical samples (rounded to 6 decimals) are evaluated only once. The credible interval of width `w` uses percentile levels from `(1-w)/2` to `1-(1-w)/2`.

The returned dictionary has `temperature_ranges`, `property_MAP`, `volumeFraction_MAP`,
`percentiles_list`, `mean_model`, `percentiles_model`, `model_output_samples`, `vol_frac_interp`
and `sobol_idx`. The Sobol' indices (next section) are computed inside this function.

## Global sensitivity analysis

`gsa/` computes **grouped first-order Sobol' indices** with the estimator of Sobol' and Myshetskaya (2008),
using the posterior samples as the input distribution. The three groups are:

1. molar fractions,
2. carbon numbers,
3. topochemical indices (isomers).

| Function | Property | Output |
|---|---|---|
| `sobol_indices` | temperature-dependent properties other than distillation | `(3, n_T)` indices along the temperature grid |
| `sobol_distillation` | distillation curve | indices along the 101 points of the curve (0-100 % recovered volume) |
| `sobol_lumped` | `molWeight`, `HC`, `DCN`, `flash`, `freezing`, `LHV`, `deltaT_dist` | one value per group (`deltaT_dist`: two columns, one for each range) |

The indices are normalised column by column so that the three groups sum to one (columns that sum
to zero are left unchanged). `sobol_indices`, `sobol_distillation` take a `seed` (default 0) for the resampling.

## Figures

`visualize_results(...)` writes PNG files (600 dpi) into `<output_dir>/Figures/`:

| Subfolder | Files |
|---|---|
| `MCMC_check/` | `Trace_X<k>`, `Trace_nc<k>`, `Trace_eta<k>`, `AcceptanceRate`, `Rhat` |
| `PosteriorPDF/MarginalProbability/` | `Marginal_X<k>`, `Marginal_nc<k>`, `Marginal_eta<k>` |
| `PosteriorPDF/Joint2D/` | `Joint2D_X<m>_X<p>`, `Joint2D_nc<m>_nc<p>` |
| `PropertyCheck_MAP/` | one figure per property: MAP and credible band against the measurements |
| `GSA/` | `Si_<property>`, `Si_lumped` |

The `AcceptanceRate` figure marks the burn-in end and has its y axis limited to 0-100 %: the initial state counts as an accepted sample, so the very first iterations can exceed 100 %. The `Rhat` figure shows the threshold `RHAT_STOP`; its y axis starts at 0.9 and goes up to the largest R-hat (at most 5).
