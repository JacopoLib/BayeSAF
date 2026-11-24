# Quick start (MATLAB): three-component surrogate of CHRJ (POSF-6152)

The example `matlab/examples/CHRJ_3comp/example_HEFA.m` builds a three-component surrogate of the HEFA-SPK fuel POSF-6152 (CHRJ).

1. Open MATLAB and set the **current folder** to `matlab/examples/CHRJ_3comp/`.
2. Open `example_HEFA.m` and run it.

On the first run MATLAB starts a parallel pool. A waitbar shows the iterations of the sampler, with a *Cancel* button.

## What the script does

| Section | Content |
|---|---|
| Paths | `addpath` of the code, of the database and of the data of the fuel |
| Real fuel data | `Dataset`, `fuel_name`, `Tbubble_realFuel` |
| Inputs of the inference | `numComponents`, `families`, `n_ranges`, bounds on `η*_B` and on the molar fractions; `classes{i} = Hydrocarbons(families{i}, n_ranges{i})` |
| Import of the data | `p_distillation`, `data_import` |
| Posterior | `pdf_build_speed` (full) and `pdf_build_speed_cheap` (first distillation point only; `deltaT_dist` is removed from its data) |
| DE-MC | parameters and the call to `DifferentialEvolutionMarkovChain` |
| Report | converged or not, maximum R-hat, minimum effective sample size (with a warning if it is below `ess_min`) |
| Post-processing | `MAPfinder`, `MAPwriter`, `props_pushforward` (figures and sensitivity analysis) |
| Saving | `save('CHRJ_3comp.mat')` |

## The configuration of the example

| Item | Value |
|---|---|
| Properties | `rho`, `distillation`, `molWeight`, `HC`, `DCN`, `LHV`, `flash` |
| Components | 3: `nparaffins` (C8-C15), `isoparaffins` (C9-C12), `isoparaffins` (C13-C16) |
| Molar-fraction bounds | 0.1 to 0.75 for each component |
| Chains | `2*(3*numComponents-1) = 16`, `n_cold = 8` |
| Maximum iterations | 20000 |
| `t_burnin` | 500 |
| Bubble temperature | 445 K |
| Not used | `deltaT_dist.csv` exists in `exp/CHRJ_POSF6152/` but is not in `Dataset` |

## Files written

In the current folder: `MCMC.txt` (log of the sampler), `MAP.txt`, `Figures/` and `CHRJ_3comp.mat`. See [Post-processing and outputs](Postprocessing.md).

## Reading the end of the run

The script prints one of

```
DE-MC converged at t = ... (burn-in: ... iterations).
DE-MC did NOT converge at t = ... (burn-in: ... iterations): increase maxIterations.
```

then the maximum R-hat of the last check and the minimum effective sample size. A run that reaches `maxIterations` has not met the R-hat criterion ([Sampler](../Sampler.md#convergence-and-stopping)).

## Trying it quickly

For a first test, reduce the cost of the run in the example: for instance `maxIterations = 1500` and `t_burnin = 100`. The run will not converge, but
you will see the full chain of operations (sampler, MAP, figures, sensitivity analysis). With a short run, the effective sample size is small and a warning is printed.

## Adapting it to another fuel

1. Create `exp/<FUEL>/` with one CSV per property ([Input data](../Input-data.md)) and replace the `addpath` of the data folder.
2. Change `Dataset`, `fuel_name` and `Tbubble_realFuel`.
3. Choose `numComponents`, `families` and `n_ranges`, then the bounds, and `N_chains` (the example uses `2*(3*numComponents-1)`).
