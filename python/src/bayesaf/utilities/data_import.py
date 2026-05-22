"""
    ____                  _____ ___    ______
   / __ )____ ___  _____ / ___//   |  / ____/
  / __  / __ `/ / / / _ \\__ \\/ /| | / /_
 / /_/ / /_/ / /_/ /  __/__/ / ___ |/ __/
/_____/\\__,_\\__, /\\___/____/_/  |_/_/
            /____/

BayeSAF: Emulation and Design of Sustainable Alternative Fuels
via Bayesian Inference and Descriptors-Based Machine Learning

Contributors / Copyright Notice
© 2026 Jacopo Liberatori — jacopo.liberatori@centralesupelec.fr
Postdoctoral Researcher @ Laboratoire EM2C, CentraleSupélec (CNRS)

© 2026 Davide Cavalieri — davide.cavalieri@uniroma1.it
Postdoctoral Researcher @ Sapienza University of Rome,
Department of Mechanical and Aerospace Engineering (DIMA)

© 2026 Matteo Blandino, Ph.D.

Reference:
J. Liberatori, D. Cavalieri, M. Blandino, M. Valorani, and P.P. Ciottoli.
BayeSAF: Emulation and Design of Sustainable Alternative Fuels via Bayesian
Inference and Descriptors-Based Machine Learning. Fuel 419, 138835 (2026).
Available at: https://doi.org/10.1016/j.fuel.2026.138835.

------------------------------------------------------------------------

Description:
The data_import module provides tools to import information about the
experimental measurements concerning the real fuel. For each thermophysical
property dataset it reads the CSV file, infers the variable name from the
filename, applies temperature threshold filtering where appropriate, and
returns the structured data arrays together with the temperature ranges.

Auxiliary parameters:
Np: number of thermophysical properties targeted during the formulation
    of the surrogate mixture  [-]

Inputs:
1) Dataset                : (1 x Np) list of strings denoting the name of
                            the experimental dataset for each property
2) families               : (1 x numComponents) list of hydrocarbon family
                            strings
3) n_ranges               : (1 x numComponents) carbon-atom ranges per
                            component
4) Tbubble_realFuel       : bubble temperature of the real fuel  [K]
5) numComponents          : number of surrogate mixture components  [-]
6) pressure_distillation  : pressure for distillation curve evaluation  [Pa]

Outputs:
1) fullData      : (1 x Np) list of (Nd x 3) or (Nd x 4) arrays of
                   experimental measurements (independent variable,
                   thermophysical property value, standard deviation)
2) minT_array    : (1 x Np) array of minimum temperatures per property
3) maxT_array    : (1 x Np) array of maximum temperatures per property
4) variable_names: (1 x Np) list of property name strings
------------------------------------------------------------------------
"""

from __future__ import annotations

from typing import Sequence

import numpy as np
import pandas as pd

from .misc import threshold_temperature

# Properties whose data are NOT filtered by temperature threshold
_NO_TEMP_FILTER = {
    "distillation", "deltaT_dist",
    "molWeight", "HC", "DCN", "flash", "freezing", "LHV",
}

# Mapping from dataset name substring → variable name
_VAR_MAP = [
    ("deltaT_dist", "deltaT_dist"),   # must come before "distillation"
    ("distillation", "distillation"),
    ("molWeight", "molWeight"),
    ("specificHeat", "specificHeat"),
    ("latentHeat", "latentHeat"),
    ("vaporPressure", "vaporPressure"),
    ("kappa", "kappa"),
    ("sigma", "sigma"),
    ("rho", "rho"),
    ("mu", "mu"),
    ("nu", "nu"),
    ("HC", "HC"),
    ("DCN", "DCN"),
    ("flash", "flash"),
    ("freezing", "freezing"),
    ("LHV", "LHV"),
]


def _detect_variable(name: str) -> str:
    for substr, var in _VAR_MAP:
        if substr in name:
            return var
    raise ValueError(f"Cannot detect variable name from dataset entry: '{name}'")


def data_import(
    dataset: list[str],
    families: list[str],
    n_ranges: list[Sequence[int]],
    Tbubble_real_fuel: float,
    num_components: int,
    pressure_distillation: float,
    data_dir: str | None = None,
) -> tuple[list[np.ndarray], list[float], list[float], list[str]]:
    """
    Import experimental thermophysical-property data.

    Parameters
    ----------
    dataset : list[str]
        Property names (e.g. ``["rho", "nu"]``) when *data_dir* is given,
        or explicit CSV file paths otherwise.
    families : list[str]
        Hydrocarbon family for each surrogate component.
    n_ranges : list of sequences
        Carbon-atom ranges for each component.
    Tbubble_real_fuel : float
        Bubble temperature of the real fuel [K].
    num_components : int
        Number of surrogate components.
    pressure_distillation : float
        Distillation pressure [Pa].
    data_dir : str or None
        Directory containing ``<name>.csv`` files.  When supplied, each entry
        in *dataset* is treated as a bare property name and the path is built
        as ``data_dir/<name>.csv``.  When *None*, entries are used as-is
        (must be complete file paths).

    Returns
    -------
    full_data : list of ndarray
        Each element is an (N_d, 3) or (N_d, 4) array
        [indep_var, property, std_dev (, pressure)].
    min_T_array : list of float
        Minimum temperature (or independent variable) per dataset.
    max_T_array : list of float
        Maximum temperature per dataset.
    variable_names : list of str
        Detected property name for each dataset entry.
    """
    import pathlib

    T_threshold = threshold_temperature(
        families, n_ranges, Tbubble_real_fuel, num_components, pressure_distillation
    )

    full_data: list[np.ndarray] = []
    min_T_array: list[float] = []
    max_T_array: list[float] = []
    variable_names: list[str] = []

    for name in dataset:
        var = _detect_variable(name)
        if data_dir is not None:
            csv_path = str(pathlib.Path(data_dir) / f"{name}.csv")
        else:
            csv_path = name
        raw = pd.read_csv(csv_path, header=0).to_numpy(dtype=float)

        if var in _NO_TEMP_FILTER:
            data = raw
        else:
            data = raw[0.999 * T_threshold > raw[:, 0], :]

        if var not in ("distillation", "deltaT_dist"):
            min_T = float(data[:, 0].min())
            max_T = float(data[:, 0].max())
        else:
            min_T = float(data[:, 1].min())
            max_T = float(data[:, 1].max())

        full_data.append(data)
        min_T_array.append(min_T)
        max_T_array.append(max_T)
        variable_names.append(var)

    return full_data, min_T_array, max_T_array, variable_names
