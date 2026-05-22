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
The find_index module returns the list index in the given hydrocarbon family
species list corresponding to the chemical species under consideration. It
first filters species that match the target number of carbon atoms nC, then
selects among them the species whose normalised topochemical atom index
eta_B_star_norm is closest to the target value. A module-level cache avoids
repeated array extractions for the same species list.

Auxiliary parameters:
Np: number of thermophysical properties targeted during surrogate
    formulation  [-]

Inputs:
1) class      : list of Species instances for the hydrocarbon family the
                chemical species belongs to, containing detailed information
                about the whole set of chemical species (number of carbon
                atoms, normalised topochemical atom index, molecular weight
                [g/mol], critical temperature [K], and Yaws' polynomial
                coefficients for all liquid-phase thermophysical properties)
2) nC         : number of carbon atoms for the chemical species  [-]
3) eta_B_star : normalised topochemical atom index for the species  [-]

Outputs:
1) row_idx: list index in the species list corresponding to the species
            under consideration
------------------------------------------------------------------------
"""

from __future__ import annotations

import numpy as np

from bayesaf.thermo_transport.hydrocarbons import Species


# Module-level cache: id(species_list) → (nC_arr, eta_arr)
# Species lists are created once per run and never mutated, so id() is stable.
_SPECIES_CACHE: dict[int, tuple[np.ndarray, np.ndarray]] = {}


def find_index_eta(species_list: list[Species], nC: int, eta_B_star: float) -> int:
    """
    Return the list index of the species matching *nC* and closest to *eta_B_star*.

    Parameters
    ----------
    species_list : list[Species]
        All candidate species for one surrogate component.
    nC : int
        Target number of carbon atoms.
    eta_B_star : float
        Target normalised topochemical atom index (eta_B_star_norm).

    Returns
    -------
    int
        0-based index into *species_list*.
    """
    key = id(species_list)
    if key not in _SPECIES_CACHE:
        _SPECIES_CACHE[key] = (
            np.array([sp.nC for sp in species_list], dtype=int),
            np.array([sp.eta_B_star_norm for sp in species_list], dtype=float),
        )
    nC_arr, eta_arr = _SPECIES_CACHE[key]

    filtered = np.where(nC_arr == nC)[0]
    diffs = np.abs(eta_arr[filtered] - eta_B_star)
    return int(filtered[np.argmin(diffs)])
