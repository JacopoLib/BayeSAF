# Known limitations and pitfalls

Observations from reading the code.

## 1. Data files are positional

`data_import` reads columns by position. A column added or reordered in a CSV changes the meaning
silently. The property is also detected by substring from the file name ([Input data](../Input-data.md#mapping-from-file-name-to-property)).

## 2. Prior bounds on `eta` do not constrain the sampler

`LOWER_ETA` and `UPPER_ETA` only appear in `log_uniform` inside the prior. The sampler explores the
full range of `u` regardless; values outside the bounds are given `-inf` and rejected.

## 3. Fixed measurement uncertainty and independent errors

The Gaussian likelihood uses the `std` of the data files as fixed standard deviations and no model error.
Underestimated `std` values make the posterior too narrow; the correlations in the distillation curve are not modelled.

## 4. Mixing rules

`vaporPressure` and the other `kappa`, `specificHeat`, `latentHeat` mixture values are
mass-weighted averages, not mechanistic mixing rules ([Thermophysical models](../Thermophysical-models.md#mixing-rules-mixture_property)).

## 5. No automated tests

The repository has no `tests/` directory. Changes to the sampler should be validated by running a short
configuration of the example with a fixed seed and comparing the chains.
