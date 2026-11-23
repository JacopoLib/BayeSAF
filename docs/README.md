# BayeSAF: documentation

BayeSAF formulates **surrogate mixtures** that emulate the physicochemical properties of a real (or to-be-designed) sustainable aviation fuel, using Bayesian inference on a database of
hydrocarbons.

> Liberatori J., Cavalieri D., Blandino M., Valorani M., Ciottoli P.P., *BayeSAF: Emulation and
> Design of Sustainable Alternative Fuels via Bayesian Inference and Descriptors-Based Machine
> Learning*, Fuel 419, 138835 (2026), https://doi.org/10.1016/j.fuel.2026.138835

The repository contains two implementations of the same method, **MATLAB** (`matlab/`) and **Python** (`python/`), which share the database (`database/`) and the experimental data (`exp/`).
This documentation describes both as they are coded, including the points where they differ and the pitfalls found while documenting them.

## In one paragraph

The user chooses the number of surrogate components and the hydrocarbon family and range of carbon numbers of each one. BayeSAF infers, for every component, a molar fraction, a number of carbon atoms
and a normalised topochemical branching index. Together they identify one species of the database. The properties of the mixture are computed with Yaws' correlations and mixing rules, and compared
with the measurements through a Gaussian likelihood. The posterior is sampled with a Differential Evolution Markov Chain with parallel tempering, delayed acceptance and an exact Gibbs-type move.
The outputs are the maximum-a-posteriori (MAP) surrogate, the posterior samples, the pushforward distributions of the properties and sensitivity indices.

## Choose your implementation

| | Python | MATLAB |
|---|---|---|
| Install | [Installation](python/Installation.md) | [Installation](matlab/Installation.md) |
| Run the example | [Quick start](python/Quick-start.md) | [Quick start](matlab/Quick-start.md) |
| Parameters | [Configuration reference](python/Configuration-reference.md) | [Configuration reference](matlab/Configuration-reference.md) |
| Results | [Post-processing and outputs](python/Postprocessing.md) | [Post-processing and outputs](matlab/Postprocessing.md) |
| Functions | [API reference](python/API-reference.md) | [API reference](matlab/API-reference.md) |
| Pitfalls | [Known limitations](python/Known-limitations.md) | [Known limitations](matlab/Known-limitations.md) |
| Contributing | [Development](python/Development.md) | |

## Common to both: data and methodology

| Page | Content |
|---|---|
| [Input data](Input-data.md) | species database and experimental files: formats, units, filters |
| [Notation and parametrisation](Notation-and-parametrisation.md) | what is inferred; normalised coordinates, decoding, decoder cells |
| [Bayesian model](Bayesian-model.md) | prior, isomer correction, likelihood of each property, cheap posterior |
| [Thermophysical models](Thermophysical-models.md) | Yaws formulas, mixing rules, flash and distillation curve |
| [Sampler](Sampler.md) | DE-MC moves, delayed acceptance, tempering, burn-in, convergence |
| [Differences between the Python and MATLAB implementations](Differences-Python-MATLAB.md) | what differs, and what is the same |

Each methodology page starts with a table that maps the concepts to the Python and MATLAB code, and ends with a section on the MATLAB implementation. The text in between
quotes Python names (`run_demc`, `log_likelihood`, ...).

## Status

Python `0.1.0` (Python 3.9 or later; numpy, scipy, pandas, matplotlib, numba). MATLAB R2022b or later with the Parallel Computing, Statistics and Machine Learning and Optimization Toolboxes.
