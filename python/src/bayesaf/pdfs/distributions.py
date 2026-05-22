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
The distributions module provides log-probability evaluation for the prior
distributions used in the Bayesian inference framework: a Dirichlet
distribution for molar fractions (log_dirichlet), a continuous uniform
distribution for normalised topochemical atom indices (log_uniform), and a
discrete uniform distribution for numbers of carbon atoms
(log_discrete_uniform). All functions return the log-probability to avoid
numerical underflow.

Inputs:
1) x         : molar fractions or parameter values to evaluate  [-]
2) alpha     : concentration parameters of the Dirichlet prior
3) lower_bound: lower bound of the support  [-]
4) upper_bound: upper bound of the support  [-]
5) n_ranges  : (1 x Nc) list of discrete carbon-atom ranges per component

Outputs:
1) log_prob: log-prior probability value  [-]
------------------------------------------------------------------------
"""

from __future__ import annotations

from typing import Sequence

import numpy as np
from scipy.special import gammaln


def log_dirichlet(
    x: np.ndarray,
    alpha: np.ndarray,
    lower_bound: np.ndarray,
    upper_bound: np.ndarray,
) -> float:
    """
    Log-PDF of the Dirichlet distribution for molar fractions.

    Parameters
    ----------
    x : ndarray, shape (Nc-1,)
        Molar fractions of the first Nc-1 components.
    alpha : ndarray, shape (Nc,)
        Concentration parameters (all > 0).
    lower_bound : ndarray, shape (Nc,)
        Lower bound for each molar fraction.
    upper_bound : ndarray, shape (Nc,)
        Upper bound for each molar fraction.

    Returns
    -------
    float
        Log-probability density (−∞ if out-of-bounds).
    """
    x = np.asarray(x, dtype=float)
    alpha = np.asarray(alpha, dtype=float)
    lower_bound = np.asarray(lower_bound, dtype=float)
    upper_bound = np.asarray(upper_bound, dtype=float)

    x_full = np.append(x, 1.0 - x.sum())

    if np.any(x_full < lower_bound) or np.any(x_full > upper_bound):
        return -np.inf

    log_const = gammaln(alpha.sum()) - gammaln(alpha).sum()
    return float(log_const + ((alpha - 1.0) * np.log(x_full)).sum())


def log_uniform(
    x: np.ndarray,
    lower_bound: np.ndarray,
    upper_bound: np.ndarray,
) -> float:
    """
    Log-PDF of independent continuous uniform distributions for eta_B_star.

    Parameters
    ----------
    x : ndarray, shape (Nc,)
        Normalised topochemical atom indices.
    lower_bound, upper_bound : ndarray, shape (Nc,)
        Bounds for each component.

    Returns
    -------
    float
        Sum of log-uniform densities (−∞ if any value is out of bounds).
    """
    x = np.asarray(x, dtype=float)
    lower_bound = np.asarray(lower_bound, dtype=float)
    upper_bound = np.asarray(upper_bound, dtype=float)

    prob = 0.0
    for j in range(len(x)):
        if x[j] < lower_bound[j] or x[j] > upper_bound[j]:
            return -np.inf
        prob -= np.log(upper_bound[j] - lower_bound[j])
    return prob


def log_discrete_uniform(
    n: np.ndarray,
    n_ranges: list[Sequence[int]],
) -> float:
    """
    Log-PDF of a discrete uniform distribution for carbon-atom counts.

    Parameters
    ----------
    n : ndarray, shape (Nc,)
        Carbon-atom counts (integers).
    n_ranges : list of sequences
        Allowed values for each component.

    Returns
    -------
    float
        Log-probability (−∞ if any value is outside its range).
    """
    n = np.asarray(n, dtype=int)
    num_outcomes = 1
    for j, rng in enumerate(n_ranges):
        if n[j] not in rng:
            return -np.inf
        num_outcomes *= len(rng)
    return -np.log(float(num_outcomes))
