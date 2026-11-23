# API reference (MATLAB)

Source in `matlab/` (about 7000 lines). One function per file, with the name of the file. For the meaning of the parameters see the linked pages.

## Data flow

```
exp/<FUEL>/*.csv ──> p_distillation, data_import ──────────┐
database/<family>/ ──> Hydrocarbons ───────────────────────┤
                                                           v
                                  pdf_build_speed / pdf_build_speed_cheap
                                                           v
                                         DifferentialEvolutionMarkovChain
                                                           v
                    MAPfinder ── MAPwriter ── props_pushforward (pushforward, figures, sobol*)
```

## `thermo-transport/`

| Function | Purpose |
|---|---|
| `species = Hydrocarbons(class, nc_range)` | struct array of the species of a family within a range of carbon numbers |
| `phi = ThermophysicalProperties_SingleLiquid(property, T, W, Tcr, Arho, ..., Bsigma, pressure)` | Yaws correlations of one species and `boilingTemperature` (28 positional arguments) |
| `phi = ThermophysicalProperties_LiquidMixture_speed(property, T, X_all, classes, index_n_eta_all, pressure)` | mixture property, `N_samples x N_T` |

## `distillation/`

| Function | Purpose |
|---|---|
| `F = Flash(T, beta, classes, index_n_eta, X, pressure)` | Rachford-Rice residual |
| `[volumeFraction, T_distillation] = DistillationCurve(X0, classes, index_n_eta, pressure)` | curve of one mixture (101 points) |
| `... = DistillationCurve_cheap(X0, classes, index_n_eta, pressure, vol1)` | same, stopped when the volume reaches `vol1` |
| `[volumeFraction_all, T_distillation_all] = DistillationCurve_speed(X_all, classes, index_n_eta_all, pressure)` | `parfor` over samples |
| `... = DistillationCurve_speed_cheap(X_all, classes, index_n_eta_all, pressure, vol1)` | cheap version |

## `pdfs/`

| Function | Purpose |
|---|---|
| `prob = dirichlet_pdf(x, alpha, LowerBound, UpperBound)` | log prior of the molar fractions |
| `prob = uniform_pdf(x, LowerBound, UpperBound)` | log prior of `eta` |
| `prob = discrete_uniform(n, n_ranges)` | log prior of the carbon numbers |
| `L = LikelihoodFunction_speed(families, classes, data, stdData, X_all, n_all, eta_B_star_all, variable, pressure)` | log-likelihood of one property |
| `L = LikelihoodFunction_speed_cheap(...)` | same, one distillation point |
| `[prior, likelihood, posterior] = pdf_build_speed(fullData, families, classes, variable_names, pressure_distillation, LowerBound_x, UpperBound_x, alpha, n_ranges, LowerBound_eta_B_star, UpperBound_eta_B_star)` | handles `@(X_all, N_all, Eta_all)` |
| `[prior, likelihood, posterior] = pdf_build_speed_cheap(...)` | same arguments, cheap likelihood |

## `mcmc/`

| Function | Purpose |
|---|---|
| `[x, p_x, x_nb, p_x_nb, AR, R_hat, t_convergence, info] = DifferentialEvolutionMarkovChain(...)` | the sampler; 29 arguments in [Configuration reference](Configuration-reference.md) |

Local functions of the same file: `isomer_tables`, `isomer_prior`, `decode_component`, `decode_rows`, `fold_fractions`, `log_sum_exp`, `rhat_classic`, `ess_bulk`.

## `postprocessing/`

| Function | Purpose |
|---|---|
| `[x_MAP, nc_MAP, eta_B_star_MAP] = MAPfinder(posterior, chain, numComponents)` | highest-posterior sample |
| `species_MAP = MAPwriter(families, numComponents, x_MAP, nc_MAP, eta_B_star_MAP, fuel_name)` | writes `MAP.txt` |
| `[temperature_range, property_MAP, volumeFraction_MAP, percentiles_list, mean_model, percentiles_model, model_output_samples, volFrac_interp, sobol_idx] = props_pushforward(fullData, families, classes, chain, AR, R_hat, chain_reshaped, t_convergence, t_burnin, n_ranges, numComponents, variable_names, x_MAP, nc_MAP, eta_B_star_MAP, minT_array, maxT_array, pressure_distillation, fuel_name, band_percentiles, confidence_width)` | pushforward, sensitivity analysis and figures |

## `gsa/`

| Function | Purpose |
|---|---|
| `sobol_idx = sobol(classes, numComponents, temperature_range, chain_reshaped, model_output_samples, variable_name, pressure_distillation)` | temperature-dependent properties |
| `sobol_idx = sobol_distillation(classes, numComponents, chain_reshaped, model_output_samples, pressure_distillation)` | distillation curve |
| `sobol_idx = sobol_lumped(families, classes, numComponents, chain_reshaped, model_output_samples, variable_name, pressure_distillation)` | lumped properties |

## `utilities/`

| Function | Purpose |
|---|---|
| `[fullData, minT_array, maxT_array, variable_names] = data_import(Dataset, families, n_ranges, Tbubble_realFuel, numComponents, pressure_distillation)` | reads the targets ([Input data](../Input-data.md)) |
| `pressure_distillation = p_distillation(Dataset)` | distillation pressure |
| `T_threshold = threshold_temperature(families, n_ranges, Tbubble_realFuel, numComponents, pressure_distillation)` | temperature filter |
| `row_idx = findIndex_eta(class, nC, eta_B_star)` | species index from `(n_C, eta)` |
| `Yi = mol2mass(Xi, Wi)`, `Wavg = molWeight(Xi, Wi)` | composition arithmetic |
| `posterior_val = posterior_handle(x, n, eta_B_star)` | loads a saved posterior handle from `posterior_handle.mat` (which must be on the path) and evaluates it |
| `h = daviolinplot(...)`, `colorList = slanCM(type, num)` | plotting utilities (violin plots, colour maps) used by the figures |
