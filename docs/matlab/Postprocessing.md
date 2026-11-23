# Post-processing and outputs (MATLAB)

## Outputs of `DifferentialEvolutionMarkovChain`

```matlab
[chain, posterior_pdf, chain_reshaped, posterior_pdf_reshaped, AR, R_hat, t_convergence, info] = ...
    DifferentialEvolutionMarkovChain(...);
```

| Output | Size | Content |
|---|---|---|
| `chain` (`x`) | `t_convergence x (3*Nc-1) x n_cold` | **physical** parameters of the cold chains (molar fractions, carbon numbers, `eta_B_star_norm`) |
| `posterior_pdf` (`p_x`) | `t_convergence x n_cold` | log-posterior (with the isomer correction) |
| `chain_reshaped` (`x_nb`) | `((t_convergence - t_burn) * n_cold) x (3*Nc-1)` | post-burn-in samples of all the cold chains stacked |
| `posterior_pdf_reshaped` (`p_x_nb`) | column vector | corresponding log-posterior |
| `AR` | `t_convergence x n_cold` | cumulative acceptance rate in percent |
| `R_hat` | `n_checks x (3*Nc-1)` | R-hat every 2 iterations after the burn-in; `NaN` for the constant `eta` columns |
| `t_convergence` | scalar | iteration of the stop, or `maxIterations` |

The structure `info` has:

| Field | Content |
|---|---|
| `t_burnin` | iteration at which the burn-in ended |
| `n_resets` | total outlier chains reset |
| `swap_rates` | swap acceptance at each level of the production ladder |
| `betas` | inverse temperatures of the production ladder |
| `ESS` | bulk effective sample size of each molar fraction, on the returned samples (`NaN` if the burn-in did not end) |
| `ESS_min` | minimum of `ESS` (diagnostic) |

The log `MCMC.txt` records the outlier resets, the end of the burn-in, the stopping message or "Convergence not reached!", and the effective sample size with a warning if
it is below `ess_min`.

## MAP surrogate

```matlab
[x_MAP, nc_MAP, eta_B_star_MAP] = MAPfinder(posterior_pdf_reshaped, chain_reshaped, numComponents);
species_MAP = MAPwriter(families, numComponents, x_MAP, nc_MAP, eta_B_star_MAP, fuel_name);
```

`MAPfinder` returns the post-burn-in sample with the highest log-posterior (`maxk` of the array it receives, taken as a column, in the same order as `chain_reshaped`). It is the highest-posterior sample visited by the sampler: it is not a marginal mode and not a mean.

`MAPwriter` writes `MAP.txt` in the current folder and returns the species names. The file contains the components, families and species,
molecular weight, molecular formula, H/C ratio, derived cetane number, flash point, freezing point, lower heating value, the molar, mass and volume fractions, the carbon numbers and topochemical indices, and
the pseudo-critical temperature and pressure (bar) from the Lee-Kesler mixing rules. The numbers are written with `num2str`. If a component has a flash point at or below 255.37 K (0 degF), the mixture flash point is `NaN` and a warning is issued.
`MAPwriter` opens the database files by the relative names `<family>/<family>_MatlabData.csv` and `<family>/<family>_Names.csv`; this works with the path set by the examples.

## Pushforward distributions, figures and sensitivity analysis

```matlab
[temperature_range, property_MAP, volumeFraction_MAP, percentiles_list, mean_model, percentiles_model, ...
 model_output_samples, volFrac_interp, sobol_idx] = props_pushforward(fullData, families, classes, chain, AR, R_hat, ...
    chain_reshaped, t_convergence, info.t_burnin, n_ranges, numComponents, variable_names, x_MAP, nc_MAP, eta_B_star_MAP, ...
    minT_array, maxT_array, pressure_distillation, fuel_name, band_percentiles, confidence_width);
```

In MATLAB a single function does the pushforward, the sensitivity analysis **and all the plots** (the equivalent of `props_pushforward` plus `visualize_results` in Python).

- **Grids.** Temperature-dependent properties: 25 temperatures between `minT_array` and `maxT_array` (if they are equal, between 233.15 K and 313.15 K), returned in degrees Celsius.
  Distillation curve: recovered volumes `[0.1, 1, 2, ..., 99, 99.9]` %.
- **Credible band.** Percentile levels from `(1-w)/2` to `1-(1-w)/2` in steps of 0.01.
- **Unique samples.** The samples are rounded to 6 decimals and each distinct one is evaluated once (the distillation curves in `parfor`).
- **Lumped properties** (`molWeight`, `HC`, `DCN`, `flash`, `freezing`, `LHV`) are evaluated once per sample.
- **Sensitivity analysis** is run inside the function ([below](#sensitivity-analysis)) and `sobol_idx` is returned.

### Figures

Written with `print(..., '-dpng', '-r600')` in `./Figures/` (relative to the current folder):

| Folder | Files |
|---|---|
| `MCMC_check/` | `Trace_X<k>`, `Trace_nc<k>`, `Trace_eta<k>` |
| `PosteriorPDF/MarginalProbability/` | `Marginal_X<k>`, `Marginal_nc<k>`, `Marginal_eta<k>` |
| `PosteriorPDF/Joint2D/` | `Joint2D_X<m>_X<p>`, `Joint2D_nc<m>_nc<p>`, `Joint2D_eta<m>_eta<p>`, `Joint2D_X<m>_nc<p>`, `Joint2D_X<m>_eta<p>`, `Joint2D_eta<m>_nc<p>` |
| `PropertyCheck_MAP/` | one figure per property (`rho`, `distillation`, `molWeight`, ...), plus `T5010` and `T9010` if `deltaT_dist` is used |
| `GSA/` | `Si_<property>` for each non-lumped property, `Si_lumped` |

The MATLAB version does **not** plot the acceptance rate and the R-hat trace; `AR` and `R_hat` are returned for you to plot them. The Python version does.
The figures use the utilities `daviolinplot.m` and `slanCM.m` included in `matlab/utilities/`.

### Sensitivity analysis

`gsa/` computes **grouped first-order Sobol' indices** with the estimator of Sobol' and Myshetskaya (2008), on the posterior samples. The groups are: molar fractions,
carbon numbers, topochemical indices.

| Function | Property | Output |
|---|---|---|
| `sobol(classes, numComponents, temperature_range, chain_reshaped, model_output_samples, variable_name, pressure_distillation)` | temperature-dependent properties | `3 x 25` |
| `sobol_distillation(classes, numComponents, chain_reshaped, model_output_samples, pressure_distillation)` | distillation curve | indices along the volume grid |
| `sobol_lumped(families, classes, numComponents, chain_reshaped, model_output_samples, variable_name, pressure_distillation)` | lumped properties and `deltaT_dist` | one value per group |

Procedure: the distinct samples (rounded to 6 decimals) are split at random into two halves A and B; for each group, the matrix AB is B with the columns of the group taken from A, and the index is
`|mean((y_A - c) .* (y_AB - y_B))| / Var(y)`, where `c` and `Var` are the mean and the variance of the outputs. The three indices of each column are then
divided by their sum, so that they add up to 1. For `mu` and `nu` the output is scaled by `1e6`, for `specificHeat` by `1e-3`.

Each function draws its random permutation from a local `RandStream('twister', 'Seed', 0)`, so the global random number generator is left untouched and the indices are reproducible.
