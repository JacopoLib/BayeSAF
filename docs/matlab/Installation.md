# Installation (MATLAB)

## Requirements

| Item | Needed for |
|---|---|
| MATLAB R2022b or later | `spmdIndex` (seeding of the parallel workers). The code also uses `readtable(..., 'VariableNamingRule', 'preserve')` (R2019b), `maxk` (R2017b) and `all(x, [2 3])` (R2018b) |
| Parallel Computing Toolbox | `parfor` (distillation curves, likelihood, sensitivity analysis, pushforward), `gcp`, `spmd` |
| Statistics and Machine Learning Toolbox | `pdist2`, `randsample`, `unifrnd`, `norminv`, `tiedrank`, `ksdensity`, `quantile`, `prctile` |
| Optimization Toolbox | `fsolve`, used for the boiling temperatures in `threshold_temperature` |

You can check the toolboxes in the Command Window:

```matlab
license('test', 'distrib_computing_toolbox')
license('test', 'statistics_toolbox')
license('test', 'optimization_toolbox')
```

Each call returns 1 if the toolbox is available. The first two were checked by the author; the requirement list is derived from the functions used in the code.

A parallel pool is started automatically at the first `parfor` (the default profile of your installation is used). The first
call is slower for this reason.

## Get the code

```bash
git clone https://github.com/JacopoLib/BayeSAF.git
```

There is nothing to install: the example adds the folders to the MATLAB path itself.

## Folder structure the code expects

The examples use relative paths from `matlab/examples/CHRJ_3comp/`:

```matlab
addpath('../../thermo-transport') ... addpath('../../utilities')      % code
addpath('../../../database/nparaffins') ...                            % species database, one folder per family
addpath('../../../exp/CHRJ_POSF6152')                                  % experimental data of the fuel
```

So:

- do not move `matlab/`, `database/` and `exp/` relative to one another;
- make the folder of the example the **current folder** before running it;
- the data files and the database files are found through the path, not through a file name with a folder.

## Where the outputs go

The files are written in the **current folder**: `MCMC.txt`, `MAP.txt`, the `Figures/` folder, and the workspace file saved by the example (`CHRJ_3comp.mat`).
Do not run from a folder that is tracked by git without checking what is created there: the `.mat` file can be large (tens of MB).
