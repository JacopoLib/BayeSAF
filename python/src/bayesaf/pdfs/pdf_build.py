r"""
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
The pdf_build module assembles vectorised prior, likelihood, and posterior
callables for use with the DE-MC sampler. Informative priors are constructed
for molar fractions (Dirichlet), numbers of carbon atoms (discrete uniform),
and normalised topochemical atom indices (continuous uniform). Each
likelihood handle calls the vectorised log-likelihood function for one
experimental dataset; all likelihoods are summed to form the total log-
likelihood. The posterior is the sum of log-prior and log-likelihood.

Inputs:
1) fullData              : (1 x Np) list of (Nd x 3) or (Nd x 4) arrays of
                           experimental measurements
2) families              : (1 x Nc) list of hydrocarbon family strings
3) classes               : (1 x Nc) list of lists of Species
4) variable_names        : (1 x Np) list of property name strings
5) pressure_distillation : pressure for distillation curve evaluation  [Pa]
6) LowerBound_x          : (1 x Nc) lower bounds for molar fractions
7) UpperBound_x          : (1 x Nc) upper bounds for molar fractions
8) alpha                 : Dirichlet concentration parameters
9) n_ranges              : (1 x Nc) carbon-atom ranges per component
10) LowerBound_eta_B_star: (1 x Nc) lower bounds for topochemical indices
11) UpperBound_eta_B_star: (1 x Nc) upper bounds for topochemical indices

Outputs:
1) prior      : vectorised log-prior callable
2) likelihood : vectorised log-likelihood callable
3) posterior  : vectorised log-posterior callable
------------------------------------------------------------------------
"""

from __future__ import annotations

from typing import Callable, Sequence

import numpy as np

from bayesaf.thermo_transport.hydrocarbons import Species
from .distributions import log_dirichlet, log_discrete_uniform, log_uniform
from .likelihood import log_likelihood, log_likelihood_cheap

# Type alias for the vectorised PDF callables:
#   fn(X_all, N_all, Eta_all) → ndarray(N_samples,)
PDFFn = Callable[[np.ndarray, np.ndarray, np.ndarray], np.ndarray]

# Standard deviations: for scalar properties use col 1 (*10), else col 2 (*10)
_SCALAR_VARS = {"molWeight", "HC", "DCN", "flash", "freezing", "LHV"}


def _std_data(data: np.ndarray, variable: str) -> np.ndarray:
    if variable in _SCALAR_VARS:
        return data[:, 1]
    return data[:, 2]


def build_pdfs(
    full_data: list[np.ndarray],
    families: list[str],
    classes: list[list[Species]],
    variable_names: list[str],
    pressure_distillation: float,
    lower_bound_x: np.ndarray,
    upper_bound_x: np.ndarray,
    alpha: np.ndarray,
    n_ranges: list[Sequence[int]],
    lower_bound_eta: np.ndarray,
    upper_bound_eta: np.ndarray,
) -> tuple[PDFFn, PDFFn, PDFFn]:
    """
    Build vectorised prior, likelihood, and posterior callables.

    Each callable has the signature::

        fn(X_all, N_all, Eta_all) → ndarray(N_samples,)

    Parameters
    ----------
    full_data : list of ndarray
        Experimental data for each property.
    families : list[str]
    classes : list of list[Species]
    variable_names : list[str]
    pressure_distillation : float
    lower_bound_x, upper_bound_x : ndarray, shape (Nc,)
        Molar-fraction bounds.
    alpha : ndarray, shape (Nc,)
        Dirichlet concentration parameters.
    n_ranges : list of sequences
        Carbon-atom ranges per component.
    lower_bound_eta, upper_bound_eta : ndarray, shape (Nc,)
        eta_B_star bounds.

    Returns
    -------
    prior, likelihood, posterior : callable
    """
    alpha = np.asarray(alpha, dtype=float)
    lower_bound_x = np.asarray(lower_bound_x, dtype=float)
    upper_bound_x = np.asarray(upper_bound_x, dtype=float)
    lower_bound_eta = np.asarray(lower_bound_eta, dtype=float)
    upper_bound_eta = np.asarray(upper_bound_eta, dtype=float)

    def prior_fn(X_all: np.ndarray, N_all: np.ndarray, Eta_all: np.ndarray) -> np.ndarray:
        N_samples = X_all.shape[0]
        log_px = np.array([
            log_dirichlet(X_all[k], alpha, lower_bound_x, upper_bound_x)
            for k in range(N_samples)
        ])
        log_pn = np.array([
            log_discrete_uniform(N_all[k].astype(int), n_ranges)
            for k in range(N_samples)
        ])
        log_pe = np.array([
            log_uniform(Eta_all[k], lower_bound_eta, upper_bound_eta)
            for k in range(N_samples)
        ])
        return log_px + log_pn + log_pe

    # Build per-property likelihood handles.
    # X_all has shape (N_samples, Nc-1): DEMC stores only the first Nc-1
    # molar fractions; log_likelihood reconstructs the last from 1-sum.
    lik_handles: list[Callable] = []
    for i, data in enumerate(full_data):
        std = _std_data(data, variable_names[i])
        var = variable_names[i]

        def _lik(
            X_all, N_all, Eta_all,
            _data=data, _std=std, _var=var,
        ) -> np.ndarray:
            return log_likelihood(
                families, classes, _data, _std,
                X_all, N_all, Eta_all, _var, pressure_distillation,
            )

        lik_handles.append(_lik)

    def likelihood_fn(X_all: np.ndarray, N_all: np.ndarray, Eta_all: np.ndarray) -> np.ndarray:
        total = np.zeros(X_all.shape[0])
        for h in lik_handles:
            total += h(X_all, N_all, Eta_all)
        return total

    def posterior_fn(X_all: np.ndarray, N_all: np.ndarray, Eta_all: np.ndarray) -> np.ndarray:
        return prior_fn(X_all, N_all, Eta_all) + likelihood_fn(X_all, N_all, Eta_all)

    return prior_fn, likelihood_fn, posterior_fn


def build_pdfs_cheap(
    full_data: list[np.ndarray],
    families: list[str],
    classes: list[list[Species]],
    variable_names: list[str],
    pressure_distillation: float,
    lower_bound_x: np.ndarray,
    upper_bound_x: np.ndarray,
    alpha: np.ndarray,
    n_ranges: list[Sequence[int]],
    lower_bound_eta: np.ndarray,
    upper_bound_eta: np.ndarray,
) -> tuple[PDFFn, PDFFn, PDFFn]:
    """
    Like :func:`build_pdfs` but uses :func:`log_likelihood_cheap` for
    ``'distillation'``: evaluates only the first distillation data point
    with an early-terminated distillation curve.
    """
    alpha = np.asarray(alpha, dtype=float)
    lower_bound_x = np.asarray(lower_bound_x, dtype=float)
    upper_bound_x = np.asarray(upper_bound_x, dtype=float)
    lower_bound_eta = np.asarray(lower_bound_eta, dtype=float)
    upper_bound_eta = np.asarray(upper_bound_eta, dtype=float)

    def prior_fn(X_all, N_all, Eta_all):
        N_samples = X_all.shape[0]
        log_px = np.array([
            log_dirichlet(X_all[k], alpha, lower_bound_x, upper_bound_x)
            for k in range(N_samples)
        ])
        log_pn = np.array([
            log_discrete_uniform(N_all[k].astype(int), n_ranges)
            for k in range(N_samples)
        ])
        log_pe = np.array([
            log_uniform(Eta_all[k], lower_bound_eta, upper_bound_eta)
            for k in range(N_samples)
        ])
        return log_px + log_pn + log_pe

    lik_handles: list[Callable] = []
    for i, data in enumerate(full_data):
        std = _std_data(data, variable_names[i])
        var = variable_names[i]

        def _lik(
            X_all, N_all, Eta_all,
            _data=data, _std=std, _var=var,
        ) -> np.ndarray:
            return log_likelihood_cheap(
                families, classes, _data, _std,
                X_all, N_all, Eta_all, _var, pressure_distillation,
            )

        lik_handles.append(_lik)

    def likelihood_fn(X_all, N_all, Eta_all):
        total = np.zeros(X_all.shape[0])
        for h in lik_handles:
            total += h(X_all, N_all, Eta_all)
        return total

    def posterior_fn(X_all, N_all, Eta_all):
        return prior_fn(X_all, N_all, Eta_all) + likelihood_fn(X_all, N_all, Eta_all)

    return prior_fn, likelihood_fn, posterior_fn
