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
The misc module provides miscellaneous utility functions used across the
BayeSAF pipeline: threshold_temperature computes the temperature below which
the lightest component in the surrogate mixture remains in liquid phase, used
to filter temperature-dependent experimental data; p_distillation extracts
the distillation pressure from the dataset list.

Inputs (threshold_temperature):
1) families              : (1 x numComponents) list of hydrocarbon family
                           strings
2) n_ranges              : (1 x numComponents) carbon-atom ranges per
                           component
3) Tbubble_realFuel      : bubble temperature of the real fuel  [K]
4) pressure_distillation : pressure for distillation curve evaluation  [Pa]

Outputs (threshold_temperature):
1) T_threshold: temperature threshold below which the lightest surrogate
                component remains in liquid phase  [K]

Inputs (p_distillation):
1) dataset  : list of CSV file names / identifiers for each property
2) data_dir : optional path to the data directory

Outputs (p_distillation):
1) pressure: distillation pressure  [Pa]
------------------------------------------------------------------------
"""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

from bayesaf.thermo_transport.hydrocarbons import load_hydrocarbons
from bayesaf.thermo_transport.properties import liquid_property


def p_distillation(dataset: list[str], data_dir: str | None = None) -> float:
    """
    Extract the distillation pressure from the dataset list.

    The pressure is stored in column 4 (index 3) of the first row of the
    ``'distillation'`` or ``'deltaT_dist'`` CSV file.  Defaults to 101325 Pa.

    Parameters
    ----------
    dataset : list[str]
        CSV file names / identifiers for each property.
    data_dir : str or None
        Directory containing ``<name>.csv`` files.  When supplied, each entry
        in *dataset* is treated as a bare property name and the path is built
        as ``data_dir/<name>.csv``.  When *None*, entries are used as-is.

    Returns
    -------
    float
        Pressure [Pa].
    """
    def _csv_path(name: str) -> str:
        if data_dir is not None:
            return str(Path(data_dir) / f"{name}.csv")
        return name

    for name in dataset:
        if "distillation" in name and "deltaT" not in name:
            df = pd.read_csv(_csv_path(name), header=0)
            return float(df.iloc[0, 3])
    for name in dataset:
        if "deltaT_dist" in name:
            df = pd.read_csv(_csv_path(name), header=0)
            return float(df.iloc[0, 3])
    return 101325.0


def threshold_temperature(
    families: list[str],
    n_ranges: list[Sequence[int]],
    Tbubble_real_fuel: float,
    num_components: int,
    pressure_distillation: float,
) -> float:
    """
    Threshold temperature for filtering experimental data.

    Defined as the minimum between *Tbubble_real_fuel* and the lowest
    boiling temperature of all candidate surrogate components.

    Parameters
    ----------
    families : list[str]
        Hydrocarbon family names for each surrogate component.
    n_ranges : list of sequences of int
        Carbon-atom ranges for each component.
    Tbubble_real_fuel : float
        Bubble temperature of the real fuel [K].
    num_components : int
        Number of surrogate mixture components.
    pressure_distillation : float
        Pressure for boiling-temperature computation [Pa].

    Returns
    -------
    float
        Threshold temperature [K].
    """
    Tb_min = np.inf
    for j in range(num_components):
        species_list = load_hydrocarbons(families[j], n_ranges[j])
        for sp in species_list:
            Tb = liquid_property("boilingTemperature", 0.0, sp, pressure_distillation)
            if Tb < Tb_min:
                Tb_min = Tb
    return min(Tb_min, Tbubble_real_fuel)
