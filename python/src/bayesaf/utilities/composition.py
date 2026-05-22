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
The composition module provides utility functions for mixture composition
arithmetic: computing the mean molecular weight of a surrogate mixture from
molar fractions and component molecular weights (mol_weight), and converting
molar fractions to mass fractions (mol_to_mass).

Inputs:
1) xi : (Nc,) array of molar fractions of the surrogate components  [-]
2) wi : (Nc,) array of molecular weights of the surrogate components  [g/mol]

Outputs:
1) mol_weight  : mean molecular weight of the mixture  [g/mol]
2) mol_to_mass : (Nc,) array of mass fractions of the surrogate components [-]
------------------------------------------------------------------------
"""

import numpy as np


def mol_weight(xi: np.ndarray, wi: np.ndarray) -> float:
    """
    Mean molecular weight of a mixture.

    Parameters
    ----------
    xi : array_like, shape (Nc,)
        Molar fractions of the surrogate components.
    wi : array_like, shape (Nc,)
        Molecular weights of the surrogate components [g/mol].

    Returns
    -------
    float
        Mean molecular weight [g/mol].
    """
    return float(np.dot(xi, wi))


def mol_to_mass(xi: np.ndarray, wi: np.ndarray) -> np.ndarray:
    """
    Convert molar fractions to mass fractions.

    Parameters
    ----------
    xi : array_like, shape (Nc,)
        Molar fractions of the surrogate components.
    wi : array_like, shape (Nc,)
        Molecular weights of the surrogate components [g/mol].

    Returns
    -------
    np.ndarray, shape (Nc,)
        Mass fractions of the surrogate components.
    """
    xi = np.asarray(xi, dtype=float)
    wi = np.asarray(wi, dtype=float)
    w_avg = np.dot(xi, wi)
    return (xi * wi) / w_avg
