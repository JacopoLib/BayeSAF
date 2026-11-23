# Known limitations and pitfalls (MATLAB)

Observations from reading the code.

## 1. Settings fixed in the code

The seed (1234), the production ladder, the outlier check, the burn-in cap and the R-hat check spacing are variables inside
`DifferentialEvolutionMarkovChain.m` and cannot be passed from the example. The log is always `MCMC.txt` in the current folder.
If you change `rhat_spacing`, keep it a multiple of 100: R-hat is stored only at multiples of 100, so any other value
means the run never stops on convergence. (The Python version checks this and raises an error.)

## 2. The name of the data files selects the property

`data_import` detects the property with `contains` on the name, with the first match in the order given in
[Input data](../Input-data.md#mapping-from-file-name-to-property); a name that matches none stops the run with an error. A file named
`mu_nu`, for example, is read as `mu`. Use exactly the names in that list.

## 3. Prior bounds on `eta` do not constrain the sampler

`LowerBound_eta_B_star` and `UpperBound_eta_B_star` act only in the prior (`uniform_pdf`). The sampler explores the full range of its coordinates; values outside the bounds
get `-Inf` and are rejected.

## 4. Fixed measurement uncertainty, independent errors

The Gaussian likelihood uses the `std` of the data files and no model error; the points of a curve are treated as independent.

## 5. Cost of the setup

`data_import` computes the temperature threshold (with `fsolve`, on every candidate species) once for **each entry** of `Dataset`, and `Hydrocarbons` reads the
database files at every call. It adds time before the sampling starts; it does not affect the result.

## 6. Parallel pool

Several steps use `parfor`; if a pool is not running, MATLAB starts one with the default profile. If the pool is closed during the run, the next `parfor` starts a new one.

## 7. No automated tests

The repository has no tests. After modifying the sampler, run a short configuration of the example and compare the chains and `MCMC.txt`.
