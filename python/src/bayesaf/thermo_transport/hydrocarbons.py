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
The hydrocarbons module provides a structure array (here a list of Species
dataclass instances) for a given hydrocarbon family and range of number of
carbon atoms. It reads species data from CSV files in the database directory
and populates each Species instance with the number of carbon atoms,
normalised topochemical atom index, molecular weight [g/mol], critical
temperature [K], and Yaws' polynomial coefficients for liquid-phase dynamic
viscosity, liquid-phase density, vapour pressure, liquid-phase thermal
conductivity, liquid-phase specific heat capacity, heat of vaporization, and
liquid-phase surface tension.

Auxiliary parameters:
Nd: number of experimental measurements for the thermophysical property
    under consideration  [-]
Np: number of thermophysical properties targeted during the formulation
    of the surrogate mixture  [-]

Inputs:
1) class   : string denoting the hydrocarbon family under consideration
2) nc_range: range of the number of carbon atoms

Outputs:
1) species: list of Species instances for the given hydrocarbon family and
            range of number of carbon atoms, containing detailed information
            about the whole set of chemical species (number of carbon atoms,
            normalised topochemical atom index, molecular weight [g/mol],
            critical temperature [K], and Yaws' polynomial coefficients for
            liquid-phase dynamic viscosity, liquid-phase density, vapour
            pressure, liquid-phase thermal conductivity, liquid-phase
            specific heat capacity, heat of vaporization, and liquid-phase
            surface tension)
------------------------------------------------------------------------
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

# Path to the CSV database relative to this file:
# src/bayesaf/thermo_transport/ → ../../../database/
_DB_ROOT = Path(__file__).parent.parent.parent.parent.parent / "database"

_FAMILY_CSV: dict[str, str] = {
    "nparaffins": "nparaffins/nparaffins_MatlabData.csv",
    "isoparaffins": "isoparaffins/isoparaffins_MatlabData.csv",
    "isoparaffins_mono_bis": "isoparaffins_mono_bis/isoparaffins_mono_bis_MatlabData.csv",
    "cycloparaffins": "cycloparaffins/cycloparaffins_MatlabData.csv",
    "dicycloparaffins": "dicycloparaffins/dicycloparaffins_MatlabData.csv",
    "alkylbenzenes": "alkylbenzenes/alkylbenzenes_MatlabData.csv",
    "alkylnaphthalenes": "alkylnaphtalenes/alkylnaphtalenes_MatlabData.csv",
    "cycloaromatics": "cycloaromatics/cycloaromatics_MatlabData.csv",
}


@dataclass
class Species:
    """
    Thermophysical data for a single hydrocarbon species.

    All Yaws' polynomial coefficients follow the naming conventions in
    ThermophysicalProperties_SingleLiquid.m.
    """

    nC: int
    eta_B_star: float
    eta_B_star_norm: float
    mol_weight: float          # [g/mol]
    Tc: float                  # critical temperature [K]
    coeff_mu: np.ndarray       # (4,) — dynamic viscosity
    coeff_rho: np.ndarray      # (3,) — liquid density
    coeff_psat: np.ndarray     # (5,) — vapour pressure
    coeff_k: np.ndarray        # (3,) — thermal conductivity
    coeff_cl: np.ndarray       # (4,) — specific heat
    coeff_hv: np.ndarray       # (2,) — latent heat
    coeff_sigma: np.ndarray    # (2,) — surface tension
    DCN: float
    Tf: float                  # flash-point temperature [K]
    Tfz: float                 # freezing temperature [K]
    Hc: float                  # lower heating value [MJ/kg]


# Cache loaded DataFrames to avoid re-reading CSV on every call
_df_cache: dict[str, pd.DataFrame] = {}


def _load_df(family: str) -> pd.DataFrame:
    if family not in _df_cache:
        csv_path = _DB_ROOT / _FAMILY_CSV[family]
        _df_cache[family] = pd.read_csv(csv_path, sep=";")
    return _df_cache[family]


def load_hydrocarbons(family: str, nc_range: Sequence[int]) -> list[Species]:
    """
    Build a list of Species for *family* filtered to *nc_range*.

    Parameters
    ----------
    family : str
        One of the keys in ``_FAMILY_CSV``
        (e.g. ``'nparaffins'``, ``'isoparaffins'``, …).
    nc_range : sequence of int
        Carbon-atom counts to include (e.g. ``range(8, 16)``).

    Returns
    -------
    list[Species]
        One entry per row in the database that matches ``nc_range``,
        ordered by *nc_range* then by row order within each nC group.
    """
    df = _load_df(family)
    nc_set = set(nc_range)
    mask = df["nC"].isin(nc_set)
    # preserve the order given by nc_range
    rows = pd.concat(
        [df[df["nC"] == nc] for nc in nc_range if nc in df["nC"].values],
        ignore_index=True,
    )

    species_list: list[Species] = []
    for _, row in rows.iterrows():
        sp = Species(
            nC=int(row["nC"]),
            eta_B_star=float(row["eta_B_star"]),
            eta_B_star_norm=float(row["eta_B_star_norm"]),
            mol_weight=float(row["W"]),
            Tc=float(row["Tc"]),
            coeff_mu=np.array([row["Amu"], row["Bmu"], row["Cmu"], row["Dmu"]], dtype=float),
            coeff_rho=np.array([row["Arho"], row["Brho"], row["Crho"]], dtype=float),
            coeff_psat=np.array([row["Asat"], row["Bsat"], row["Csat"], row["Dsat"], row["Esat"]], dtype=float),
            coeff_k=np.array([row["Ak"], row["Bk"], row["Ck"]], dtype=float),
            coeff_cl=np.array([row["Ac"], row["Bc"], row["Cc"], row["Dc"]], dtype=float),
            coeff_hv=np.array([row["Avap"], row["Bvap"]], dtype=float),
            coeff_sigma=np.array([row["Asigma"], row["Bsigma"]], dtype=float),
            DCN=float(row["DCN"]),
            Tf=float(row["Tf"]),
            Tfz=float(row["Tfz"]),
            Hc=float(row["Hc"]),
        )
        species_list.append(sp)

    return species_list
