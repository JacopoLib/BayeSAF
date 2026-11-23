# Quick start: three-component surrogate of CHRJ (POSF-6152)

The example `python/examples/CHRJ_3comp/example_hefa.py` builds a three-component surrogate of the
HEFA-SPK fuel POSF-6152 (CHRJ).

```bash
cd python/examples/CHRJ_3comp
python example_hefa.py
```

It writes `MAP.txt`, `demc_log.txt` and a `Figures/` directory in the current folder.

## What the script does

1. **Species databases**: `load_hydrocarbons(family, nC_range)` loads, for each component, the
   candidate species of its family within a carbon-number range.
2. **Experimental data**: `p_distillation` reads the pressure of the distillation curve and
   `data_import` assembles one array per target property.
3. **Posterior**: `build_pdfs` returns prior, likelihood and posterior callables.
   `build_pdfs_cheap` builds the cheap posterior used for delayed acceptance: it uses only the first
   distillation point ([Bayesian model](../Bayesian-model.md#cheap-likelihood-and-cheap-posterior)). The script also
   removes `deltaT_dist` from the cheap data set, which has no effect here because `DATASET` does not contain it.
4. **Sampling**: `run_demc(...)`. See [Configuration reference](Configuration-reference.md).
5. **MAP**: `find_map` and `write_map` (writes `MAP.txt`).
6. **Pushforward and figures**: `props_pushforward` and `visualize_results`.

## The configuration of the example

| Item | Value |
|---|---|
| Target properties | `rho`, `distillation`, `molWeight`, `HC`, `DCN`, `LHV`, `flash` |
| Components | 3 |
| Families | `nparaffins`, `isoparaffins`, `isoparaffins` |
| Carbon ranges | C8-C15, C9-C12, C13-C16 |
| Molar-fraction bounds | 0.1 to 0.75 for each component |
| Chains | `2*(3*Nc-1) = 16` |
| Maximum iterations | 20000 |
| Burn-in | 500 |
| Bubble temperature | 445 K |
| Not used | `exp/CHRJ_POSF6152/deltaT_dist.csv` exists but `deltaT_dist` is not in `DATASET` |

## Pointers

- The meaning of each parameter: [Configuration reference](Configuration-reference.md).
- What the sampler does in each phase: [Sampler](../Sampler.md).
- How to read `MAP.txt` and the figures: [Post-processing and outputs](Postprocessing.md).
- The script must be run from its own directory (relative path to `exp/`) and contains the
  `if __name__ == "__main__":` guard needed by the worker pool.

## Reading the end of the run

The script prints whether DE-MC converged (`t_convergence < MAX_ITERATIONS`), the maximum R-hat
at the last check and the minimum effective sample size (with a warning if it is below `ESS_MIN`). If the stopping rule
(R-hat) is not satisfied, it suggests increasing `MAX_ITERATIONS`; see [Sampler](../Sampler.md#convergence-and-stopping).

## Adapting it to another fuel

1. Create `exp/<FUEL>/` with one CSV per property ([Input data](../Input-data.md)).
2. Change `DATASET`, `DATA_DIR`, `FUEL_NAME` and `T_BUBBLE`.
3. Choose `NUM_COMPONENTS`, `FAMILIES` and `N_RANGES`, then adjust bounds and `N_CHAINS`
   (the example uses `2*(3*Nc-1)`).
