# Bayesian model

> **Implementation map**
>
> | Concept | Python (`python/src/bayesaf/pdfs/`) | MATLAB (`matlab/pdfs/`) |
> |---|---|---|
> | Dirichlet prior | `log_dirichlet` | `dirichlet_pdf.m` |
> | Isomer prior (`eta`) | `log_uniform` | `uniform_pdf.m` |
> | Carbon-number prior | `log_discrete_uniform` | `discrete_uniform.m` |
> | Isomer correction | `_uniform_isomer_prior` in `mcmc/demc.py` | `isomer_prior`, `isomer_tables` in `DifferentialEvolutionMarkovChain.m` |
> | Likelihood of one property | `log_likelihood` | `LikelihoodFunction_speed.m` |
> | Cheap likelihood | `log_likelihood_cheap` | `LikelihoodFunction_speed_cheap.m` |
> | Posterior handles | `build_pdfs`, `build_pdfs_cheap` | `pdf_build_speed.m`, `pdf_build_speed_cheap.m` |
>
> The model is the same in both implementations. The differences are listed in
> [Differences between the Python and MATLAB implementations](Differences-Python-MATLAB.md).

Prior, likelihood and posterior as implemented in `pdfs/distributions.py`, `pdfs/pdf_build.py`
and `pdfs/likelihood.py`. All densities are handled in log scale and vectorised over samples:
every callable has the signature `fn(X_all, N_all, Eta_all) -> ndarray(N_samples,)`, with
`X_all` of shape `(N, Nc-1)`, and `N_all`, `Eta_all` of shape `(N, Nc)`.

## Posterior

```math
\log p(\lambda \mid \mathcal{D}) = \log p(\lambda) + \sum_{k \in \text{properties}} \log L_k(\lambda) + \text{const}
```

`build_pdfs(...)` returns `(prior_fn, likelihood_fn, posterior_fn)`. The likelihood is the sum of
one handle per property in `DATASET`.

## Prior

`prior_fn` is the sum of three terms.

**Molar fractions: Dirichlet.** `log_dirichlet(x, alpha, lower, upper)`:

```math
\log p(x) = \log\Gamma\Big(\sum_j \alpha_j\Big) - \sum_j \log\Gamma(\alpha_j) + \sum_j (\alpha_j - 1)\log x_j
```

evaluated on the full vector `[x_1..x_{Nc-1}, 1 - sum]`. It returns `-inf` if any fraction is
outside `[lower_j, upper_j]`. The density is **not** renormalised for the truncation, so with
tight bounds its value is a constant offset that does not affect the sampling.

**Carbon number: discrete uniform.** `log_discrete_uniform(n, n_ranges)` is `-log(prod_j |range_j|)`
if every `n_j` is in the corresponding `N_RANGES` entry, `-inf` otherwise.

**Isomer index: continuous uniform.** `log_uniform(eta, lower, upper)` is `-sum_j log(upper_j - lower_j)`
inside the bounds, `-inf` outside. With the default bounds `[0, 1]` it is a constant. Narrower
bounds exclude species whose `eta_B_star_norm` falls outside them.

### Isomer prior correction in the sampler

The sampler works in the `u` coordinates of [Notation and parametrisation](Notation-and-parametrisation.md#decoder-cells),
where a uniform proposal gives each isomer a probability equal to the measure `m` of its decoder
cell. `run_demc` wraps both posteriors (full and cheap) with a correction that makes the prior
**uniform over the reachable isomers of each carbon number**:

```math
\text{correction}_j = -\log(\#\text{reachable isomers at } n_{C,j}) - \log m_j
```

summed over components. It is `-inf` for species that are unreachable. This correction is added
inside `run_demc`: `build_pdfs` alone does not include it. When using `posterior_fn` outside the
sampler (for example to evaluate a given species), remember that it is the *uncorrected*
posterior.

## Likelihood

`log_likelihood(...)` evaluates, for one property, a Gaussian likelihood of the measurements
`d_i` with standard deviations `s_i` given the model value `phi_i(lambda)`:

```math
\log L = -\tfrac{1}{2} \sum_i \Big[\log(2\pi s_i^2) + \frac{(\phi_i - d_i)^2}{s_i^2}\Big]
```

(`_gaussian_log_likelihood`). The measurements are treated as independent: there is no model-error
term and no correlation between points. The standard deviations are taken from the data files
(see [Input data](Input-data.md)); they are not inferred.

The per-property model values are:

| `variable` | Model value `phi` |
|---|---|
| `molWeight` | `sum_j x_j W_j` |
| `HC` | atomic H/C ratio from the family formulas ([below](#hc-ratio)) |
| `DCN` | volume-weighted average of the component DCN |
| `LHV` | mass-weighted average of `Hc` |
| `flash` | blending index method of ASTM D-7215 on component flash points |
| `freezing` | volume-weighted blend of `Tfz^(1/0.05)`, raised to `0.05` |
| `distillation` | temperature of the computed distillation curve at the measured recovered volume |
| `deltaT_dist` | `T50 - T10` (code `5010`) or `T90 - T10` (code `9010`) of the computed curve |
| `rho`, `mu`, `nu`, `kappa`, `specificHeat`, `latentHeat`, `vaporPressure`, `sigma` | mixture property at the measured temperatures ([Thermophysical models](Thermophysical-models.md)) |

Details:

- **Index lookup.** For every sample and component, the species is found with `find_index_eta(classes[j], nC, eta)`.
- **HC ratio.** Hydrogen atoms per molecule are `2n+2` (n- and isoparaffins), `2n` (cycloparaffins),
  `2n-2` (dicycloparaffins), `2n-6` (alkylbenzenes), `2n-12` (alkylnaphthalenes), `2n-8` (cycloaromatics),
  averaged with the molar fractions. The mean carbon number is rounded to 2 decimals and the ratio to
  3 decimals before the Gaussian is evaluated.
- **DCN.** Volume fractions are `V_i = rho_mix * Y_i / rho_i`, using the liquid densities at 300 K and
  the distillation pressure and mass fractions `Y_i`; `rho_mix` is the ideal mixture density.
- **LHV.** `sum_j Y_j Hc_j`, with `Y_j` the mass fractions and `Hc` in MJ/kg.
- **Flash point.** For each component, `T_f` is converted to degrees Fahrenheit,
  `BI_j = 51708 * exp((ln F_j - 2.6287)^2 / (-0.91725))`. The blending index is the mass-weighted
  sum, converted back with `F = exp(sqrt(-0.91725 ln(BI/51708)) + 2.6287)` and then to kelvin.
  If any component has `T_f <= 0 degF` the method is undefined and the likelihood is `-inf`.
- **Freezing point.** Volume-weighted (same volume fractions as DCN) blend of `Tfz^(1/0.05)`; the
  result is raised to `0.05`.
- **Distillation.** The full curve is computed for each sample ([Thermophysical models](Thermophysical-models.md#distillation-curve)),
  then linearly interpolated at the measured recovered volumes with linear extrapolation outside the
  curve. Any non-finite interpolated value gives `-inf`.
- **`deltaT_dist`.** Computes its own distillation curves in the Python implementation: the cache that would
  share them with the `distillation` handle is created anew at every call and is never reused. The MATLAB implementation
  does share them when the batch is the same (see the box below).
- **Temperature-dependent properties.** Evaluated at the measured temperatures with `mixture_property`.

## Cheap likelihood and cheap posterior

`build_pdfs_cheap` builds a second posterior with `log_likelihood_cheap`. It is identical to the
full one except for `distillation`: only the **first** distillation data point is used, with
a distillation curve that is stopped as soon as the recovered volume reaches that point.
This reduces the cost of one evaluation substantially. The sampler uses it for the first stage of
delayed acceptance ([Sampler](Sampler.md#delayed-acceptance)).

The CHRJ example script removes `deltaT_dist` from the data passed to `build_pdfs_cheap` (it has no
effect there because that dataset is not in its `DATASET`). If you pass `deltaT_dist` to
`build_pdfs_cheap`, it is evaluated in full, with the complete distillation curve. If your `DATASET` includes
`deltaT_dist`, remove it from the cheap data set as the script does, otherwise the first stage loses
its cost advantage.

The cheap posterior does not need to be accurate. It needs to be a reasonable proposal filter:
the second stage corrects the difference exactly, so the chain samples the full posterior.

## What the model does not include

- No inference of the measurement uncertainties: they are fixed by the data files.
- No term for model error of the property correlations or of the mixing rules.
- No correlation between measurements of the same property.
- No prior information on the composition beyond the bounds and the Dirichlet concentration. Informative
  priors from GCxGC can be introduced through `ALPHA_DIRICHLET` and the bounds.

## MATLAB implementation: notes

These points apply to `matlab/pdfs/` only.

- **Property variables.** `LikelihoodFunction_speed` handles `molWeight`, `HC`, `DCN`, `LHV`, `flash`, `freezing`,
  `distillation`, `deltaT_dist` and, in the `otherwise` branch, any property accepted by
  `ThermophysicalProperties_LiquidMixture_speed` (`rho`, `mu`, `nu`, `kappa`, `specificHeat`, `latentHeat`,
  `vaporPressure`, `sigma`).
- **Standard deviation column.** As in Python: column 2 for `molWeight`, `HC`, `DCN`, `flash`, `freezing`, `LHV`;
  column 3 for the others (`pdf_build_speed.m`).
- **H/C ratio.** The hydrogen count uses a `switch` on the family; a family without a matching case stops with an error
  (`Unknown hydrocarbon family for the H/C ratio`).
- **Flash point.** As in Python, a sample with a component whose flash point is at or below 255.37 K (0 degF) gets a likelihood of `-Inf`
  (ASTM D-7215 is undefined there).
- **`deltaT_dist` and the cache.** The function keeps a `persistent` cache that `distillation` fills and `deltaT_dist` reuses, provided that the pressure, the
  molar fractions and the species indices of the current batch are the same as those stored. Otherwise `deltaT_dist` computes its own curves. The order of the two properties in `Dataset`
  therefore does not matter.
- **Cheap likelihood.** The cheap function keeps its own cache. In the example `deltaT_dist` is removed from the cheap data
  set, so the cheap posterior contains one distillation point and no `deltaT_dist`.
- **Parallelism.** The interpolation of the distillation curves and the curves themselves run in `parfor`.
