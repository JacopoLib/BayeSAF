# Installation

## Requirements

- Python `>=3.9`
- numpy `>=1.24`, scipy `>=1.10`, pandas `>=2.0`, matplotlib `>=3.7`
- optional (development): pytest, pytest-cov
- numba `>=0.58`: installed automatically. It compiles the root finder of the distillation curve, which
  is the most expensive part of the likelihood. If the import fails the code still runs with SciPy, but more slowly.

## Install from the repository

```bash
git clone git@<host-gitlab>:<group>/saf-uq/bayesaf.git
cd bayesaf/python
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
```

The package is installed in editable mode, so changes to `python/src/bayesaf/` take effect
immediately.

## Check the installation

```bash
python -c "import bayesaf; print(bayesaf.__version__)"
```

## Data location

The species database is located by `thermo_transport/hydrocarbons.py` relative to the source file
(`_DB_ROOT` is the `database/` directory at the repository root). For this reason:

- install in **editable mode** (`pip install -e`) from a clone of the repository; a regular
  install that copies the package into `site-packages` will not find `database/`;
- do not move `database/` away from the repository root.

The examples also read `exp/` from the repository root (see [Input data](../Input-data.md)) and use
relative paths such as `../../../exp/CHRJ_POSF6152`, so run them from their own directory.

## Common issues

- **`FileNotFoundError` on a CSV**: you are not running from the example directory, or the
  `DATA_DIR` path is wrong.
- **Slow runs**: the distillation curve is the expensive term; see
  [Sampler](../Sampler.md#delayed-acceptance) and [Configuration reference](Configuration-reference.md).
