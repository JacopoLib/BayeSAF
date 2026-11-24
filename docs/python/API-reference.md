# API reference

Source in `python/src/bayesaf/` (about 7000 lines). Signatures are those of the code; for the meaning
of the parameters see the linked pages.

## Data flow

```
exp/<FUEL>/*.csv ──> p_distillation, data_import ──────────┐
database/<family>/ ──> load_hydrocarbons ──────────────────┤
                                                           v
                                       build_pdfs / build_pdfs_cheap
                                                           v
                                                       run_demc
                                                           v
        find_map ── write_map ── props_pushforward ── visualize_results
                                       └─ gsa.sobol_*
```

## `bayesaf.thermo_transport`

| Function / class | Purpose |
|---|---|
| `Species` | dataclass with the data of one species ([Notation](../Notation-and-parametrisation.md#the-species-list-of-a-component)) |
| `load_hydrocarbons(family, nc_range) -> list[Species]` | species of a family within a range of carbon numbers |
| `liquid_property(prop, T, sp, pressure=101325.0)` | pure-species Yaws correlations and `boilingTemperature` |
| `mixture_property(prop, T, X_all, classes, index_n_eta_all, pressure=101325.0)` | mixture property, shape `(N_samples, N_T)` |

Details: [Thermophysical models](../Thermophysical-models.md).

## `bayesaf.distillation`

| Function | Purpose |
|---|---|
| `flash_residual(T, beta, classes, index_n_eta, X, pressure)` | Rachford-Rice residual |
| `distillation_curve(X0, classes, index_n_eta, pressure=101325.0)` | `(volume_fraction[101], T[101])` for one mixture |
| `distillation_curve_cheap(..., vol1)` | same, stopped when the recovered volume reaches `vol1`; the rest is `NaN` |
| `distillation_curve_batch(X_all, classes, index_n_eta, pressure)` | parallel version over samples |
| `distillation_curve_batch_cheap(X_all, classes, index_n_eta, pressure, vol1)` | parallel cheap version |
| `init_pool(n_workers=None, classes=None, batch_size=None, nc=None)` | start the persistent worker pool |
| `shutdown_pool()` | release the pool and the shared memory |

## `bayesaf.pdfs`

| Function | Purpose |
|---|---|
| `log_dirichlet(x, alpha, lower_bound, upper_bound)` | prior of the molar fractions |
| `log_uniform(x, lower_bound, upper_bound)` | prior of `η*_B` |
| `log_discrete_uniform(n, n_ranges)` | prior of the carbon numbers |
| `log_likelihood(families, classes, data, std_data, X_all, n_all, eta_B_star_all, variable, pressure)` | Gaussian log-likelihood of one property |
| `log_likelihood_cheap(...)` | same, with the first-point distillation |
| `build_pdfs(full_data, families, classes, variable_names, pressure_distillation, lower_bound_x, upper_bound_x, alpha, n_ranges, lower_bound_eta, upper_bound_eta)` | returns `(prior, likelihood, posterior)` callables |
| `build_pdfs_cheap(...)` | same arguments, cheap likelihood |

Details: [Bayesian model](../Bayesian-model.md).

## `bayesaf.mcmc`

| Function / class | Purpose |
|---|---|
| `run_demc(families, posterior_fn, posterior_cheap_fn, classes, lower_bound_x, upper_bound_x, n_ranges, lower_bound_eta, upper_bound_eta, max_iterations, t_burnin, scaling_factor_x, scaling_factor_nc, scaling_factor_eta, noise_x, noise_nc, noise_eta, N_chains, ...)` | the sampler; full argument list in [Configuration reference](Configuration-reference.md) |
| `ParameterLayout`, `SurrogateLayout`, `CompositionUnit`, `Block` (`mcmc/layout.py`) | description of the state vector and of the proposal blocks |

`ParameterLayout` is the extension point: a different model with a different state vector can provide its
own subclass (`units`, `blocks`, `n_params`, `posterior_args`, `init_population`) and pass it as `layout`.

## `bayesaf.postprocessing`

| Function | Purpose |
|---|---|
| `find_map(posterior_pdf, chain_reshaped, num_components)` | sample with the highest log-posterior |
| `write_map(families, num_components, x_MAP, nc_MAP, eta_B_star_MAP, fuel_name, output_file="MAP.txt", data_dir=None)` | `MAP.txt` |
| `props_pushforward(full_data, families, classes, chain_reshaped, num_components, variable_names, x_MAP, nc_MAP, eta_B_star_MAP, min_T_array, max_T_array, pressure_distillation, confidence_width=0.9)` | pushforward statistics and Sobol' indices |
| `visualize_results(full_data, families, classes, chain, AR, R_hat, chain_reshaped, t_burnin, R_hat_threshold, n_ranges, num_components, variable_names, pushforward_results, fuel_name, confidence_width=0.9, band_percentiles=False, output_dir=".")` | all figures |

Details: [Post-processing and outputs](Postprocessing.md).

## `bayesaf.gsa`

| Function | Purpose |
|---|---|
| `sobol_indices(classes, num_components, temperature_range, chain_reshaped, model_output_samples, variable_name, pressure_distillation, seed=0)` | temperature-dependent properties, `(3, n_T)` |
| `sobol_distillation(classes, num_components, chain_reshaped, model_output_samples, pressure_distillation, seed=0)` | distillation curve |
| `sobol_lumped(families, classes, num_components, chain_reshaped, model_output_samples, variable_name, pressure_distillation)` | lumped properties |

## `bayesaf.utilities`

| Function | Purpose |
|---|---|
| `data_import(dataset, families, n_ranges, Tbubble_real_fuel, num_components, pressure_distillation, data_dir=None)` | read the targets ([Input data](../Input-data.md)) |
| `p_distillation(dataset, data_dir=None)` | distillation pressure from the data |
| `threshold_temperature(families, n_ranges, Tbubble_real_fuel, num_components, pressure_distillation)` | temperature filter |
| `find_index_eta(species_list, nC, eta_B_star)` | species index from `(n_C, η*_B)` |
| `mol_weight(xi, wi)`, `mol_to_mass(xi, wi)` | composition arithmetic |
