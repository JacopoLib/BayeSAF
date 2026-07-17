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
The likelihood module provides a vectorised Gaussian log-likelihood function
for multiple MCMC samples. For each sample (defined by molar fractions X_all,
numbers of carbon atoms n_all, and normalised topochemical atom indices
eta_B_star_all) it evaluates the thermophysical property model and computes
the sum of squared standardised residuals with respect to the experimental
data. Supported properties include: molecular weight, H/C ratio, derived
cetane number (DCN), flash point, freezing point, lower heating value (LHV),
liquid density, dynamic viscosity, kinematic viscosity, thermal conductivity,
specific heat capacity, latent heat of vaporization, vapour pressure, surface
tension, and distillation curve.

Inputs:
1) families        : (1 x Nc) list of strings denoting the hydrocarbon family
                     of each surrogate component
2) classes         : (1 x Nc) list of lists of Species — candidate species
                     per surrogate component
3) data            : (Nd x 3) or (Nd x 4) array of experimental measurements
                     (independent variable, property value, standard deviation)
4) stdData         : (Nd,) array of standard deviations for each measurement
5) X_all           : (numSamples x Nc-1) array of molar fractions
6) n_all           : (numSamples x Nc) array of numbers of carbon atoms
7) eta_B_star_all  : (numSamples x Nc) array of normalised topochemical
                     atom indices
8) variable        : string denoting the thermophysical property
9) pressure        : pressure at which the distillation curve is computed  [Pa]

Outputs:
1) L: (numSamples x 1) array of log-likelihood values
------------------------------------------------------------------------
"""

from __future__ import annotations

import numpy as np

from bayesaf.thermo_transport.hydrocarbons import Species
from bayesaf.thermo_transport.properties import liquid_property, mixture_property
from bayesaf.utilities.composition import mol_to_mass
from bayesaf.utilities.find_index import find_index_eta
from bayesaf.distillation.distillation_curve import (
    distillation_curve_batch,
    distillation_curve_batch_cheap,
)


def _interp_extrap(xq: float | np.ndarray, xp: np.ndarray, fp: np.ndarray) -> np.ndarray:
    """
    1-D linear interpolation with linear extrapolation outside [xp[0], xp[-1]].

    Matches MATLAB's ``interp1(xp, fp, xq, 'linear', 'extrap')``.
    Only finite (xp, fp) pairs are used; returns NaN when fewer than 2
    valid points are available.

    Parameters
    ----------
    xq : float or array-like
        Query points.
    xp : array-like, shape (N,)
        Data x-values (must be monotonically increasing after filtering).
    fp : array-like, shape (N,)
        Data y-values.

    Returns
    -------
    ndarray, shape matching xq
    """
    xq = np.atleast_1d(np.asarray(xq, dtype=float))
    xp = np.asarray(xp, dtype=float)
    fp = np.asarray(fp, dtype=float)

    valid = np.isfinite(xp) & np.isfinite(fp)
    if valid.sum() < 2:
        return np.full(xq.shape, np.nan)

    xp_v = xp[valid]
    fp_v = fp[valid]

    # Interior: standard linear interpolation (no extrapolation)
    result = np.interp(xq, xp_v, fp_v)

    # Left extrapolation
    left = xq < xp_v[0]
    if left.any():
        slope = (fp_v[1] - fp_v[0]) / (xp_v[1] - xp_v[0])
        result[left] = fp_v[0] + slope * (xq[left] - xp_v[0])

    # Right extrapolation
    right = xq > xp_v[-1]
    if right.any():
        slope = (fp_v[-1] - fp_v[-2]) / (xp_v[-1] - xp_v[-2])
        result[right] = fp_v[-1] + slope * (xq[right] - xp_v[-1])

    return result


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _gaussian_log_likelihood(
    phi_model: np.ndarray,   # (N_samples,) or (N_samples, N_exp)
    data_vals: np.ndarray,   # (N_exp,)
    std_vals: np.ndarray,    # (N_exp,)
) -> np.ndarray:
    """Vectorised Gaussian log-likelihood summed over experimental points."""
    # Ensure phi_model is 2-D: (N_samples, N_exp)
    if phi_model.ndim == 1:
        phi_model = phi_model[:, np.newaxis]
    data_mat = data_vals[np.newaxis, :]
    std_mat = std_vals[np.newaxis, :]
    arr = np.log(2 * np.pi * std_mat**2) + ((phi_model - data_mat)**2) / std_mat**2
    return -0.5 * arr.sum(axis=1)


def _build_index_matrix(
    classes: list[list[Species]],
    n_all: np.ndarray,
    eta_B_star_all: np.ndarray,
) -> np.ndarray:
    """Compute index_n_eta[k, j] for all samples k and components j."""
    N_samples, Nc = n_all.shape
    idx = np.empty((N_samples, Nc), dtype=int)
    for j in range(Nc):
        for k in range(N_samples):
            idx[k, j] = find_index_eta(classes[j], int(n_all[k, j]), float(eta_B_star_all[k, j]))
    return idx


def _h_to_carbon_ratio(
    families: list[str],
    n_all: np.ndarray,
    mol_frac_mat: np.ndarray,
) -> np.ndarray:
    """H/C atom ratio for each sample."""
    N_samples, Nc = mol_frac_mat.shape
    n_hydrogen = np.zeros(N_samples)
    for j in range(Nc):
        fam = families[j]
        if fam in ("nparaffins", "isoparaffins"):
            n_hydrogen += mol_frac_mat[:, j] * (2 * n_all[:, j] + 2)
        elif fam == "cycloparaffins":
            n_hydrogen += mol_frac_mat[:, j] * (2 * n_all[:, j])
        elif fam == "dicycloparaffins":
            n_hydrogen += mol_frac_mat[:, j] * (2 * n_all[:, j] - 2)
        elif fam == "alkylbenzenes":
            n_hydrogen += mol_frac_mat[:, j] * (2 * n_all[:, j] - 6)
        elif fam == "alkylnaphtalenes":
            n_hydrogen += mol_frac_mat[:, j] * (2 * n_all[:, j] - 12)
        elif fam == "cycloaromatics":
            n_hydrogen += mol_frac_mat[:, j] * (2 * n_all[:, j] - 8)
    n_carbon = np.round((mol_frac_mat * n_all).sum(axis=1), 2)
    return np.round(n_hydrogen / n_carbon, 3)


# ---------------------------------------------------------------------------
# Main likelihood
# ---------------------------------------------------------------------------

def log_likelihood(
    families: list[str],
    classes: list[list[Species]],
    data: np.ndarray,
    std_data: np.ndarray,
    X_all: np.ndarray,
    n_all: np.ndarray,
    eta_B_star_all: np.ndarray,
    variable: str,
    pressure: float,
    *,
    _cache: dict | None = None,
) -> np.ndarray:
    """
    Vectorised Gaussian log-likelihood for multiple MCMC samples.

    Parameters
    ----------
    families : list[str]
        Hydrocarbon family names.
    classes : list of list[Species]
    data : ndarray, shape (N_exp, ≥2)
        Experimental data columns: [indep_var, property, std, …].
    std_data : ndarray, shape (N_exp,)
        Standard deviations.
    X_all : ndarray, shape (N_samples, Nc-1)
    n_all : ndarray, shape (N_samples, Nc)
    eta_B_star_all : ndarray, shape (N_samples, Nc)
    variable : str
        Property name.
    pressure : float
    _cache : dict, optional
        Mutable dict used to cache distillation curves between calls.

    Returns
    -------
    ndarray, shape (N_samples,)
        Log-likelihood for each sample.
    """
    if _cache is None:
        _cache = {}

    N_samples, Nc_m1 = X_all.shape
    Nc = Nc_m1 + 1
    std_data = std_data.ravel()

    mol_frac_mat = np.empty((N_samples, Nc), dtype=float)
    mol_frac_mat[:, :Nc_m1] = X_all
    mol_frac_mat[:, Nc_m1] = 1.0 - X_all.sum(axis=1)

    index_n_eta = _build_index_matrix(classes, n_all, eta_B_star_all)

    # ── scalar / non-temperature properties ──────────────────────────────────
    if variable == "molWeight":
        W_mat = np.array(
            [[classes[j][index_n_eta[k, j]].mol_weight for j in range(Nc)]
             for k in range(N_samples)]
        )
        phi = (mol_frac_mat * W_mat).sum(axis=1)
        return _gaussian_log_likelihood(phi, data[:, 0], std_data)

    elif variable == "HC":
        phi = _h_to_carbon_ratio(families, n_all, mol_frac_mat)
        return _gaussian_log_likelihood(phi, data[:, 0], std_data)

    elif variable == "DCN":
        W_mat = np.zeros((N_samples, Nc))
        dcn_mat = np.zeros((N_samples, Nc))
        rho_mat = np.zeros((N_samples, Nc))
        for j in range(Nc):
            for k in range(N_samples):
                sp = classes[j][index_n_eta[k, j]]
                W_mat[k, j] = sp.mol_weight
                dcn_mat[k, j] = sp.DCN
                rho_mat[k, j] = liquid_property("rho", 300.0, sp, pressure)
        rho_mix = (mol_frac_mat * W_mat).sum(axis=1) / ((mol_frac_mat * W_mat) / rho_mat).sum(axis=1)
        Y_i = np.array([mol_to_mass(mol_frac_mat[k], W_mat[k]) for k in range(N_samples)])
        V_i = rho_mix[:, np.newaxis] * Y_i / rho_mat
        phi = (V_i * dcn_mat).sum(axis=1)
        return _gaussian_log_likelihood(phi, data[:, 0], std_data)

    elif variable == "LHV":
        W_mat = np.zeros((N_samples, Nc))
        lhv_mat = np.zeros((N_samples, Nc))
        for j in range(Nc):
            for k in range(N_samples):
                sp = classes[j][index_n_eta[k, j]]
                W_mat[k, j] = sp.mol_weight
                lhv_mat[k, j] = sp.Hc
        Y_i = np.array([mol_to_mass(mol_frac_mat[k], W_mat[k]) for k in range(N_samples)])
        phi = (Y_i * lhv_mat).sum(axis=1)
        return _gaussian_log_likelihood(phi, data[:, 0], std_data)

    elif variable == "flash":
        W_mat = np.zeros((N_samples, Nc))
        bi_flash = np.zeros((N_samples, Nc))
        has_invalid_flash = np.zeros(N_samples, dtype=bool)
        for j in range(Nc):
            for k in range(N_samples):
                sp = classes[j][index_n_eta[k, j]]
                W_mat[k, j] = sp.mol_weight
                flash_F = (sp.Tf - 273.15) * 9.0 / 5.0 + 32.0
                if flash_F > 0.0:
                    bi_flash[k, j] = 51708.0 * np.exp((np.log(flash_F) - 2.6287)**2 / (-0.91725))
                else:
                    has_invalid_flash[k] = True  # ASTM D-7215 undefined for Tf ≤ 0 °F
        Y_i = np.array([mol_to_mass(mol_frac_mat[k], W_mat[k]) for k in range(N_samples)])
        bi_blend = (Y_i * bi_flash).sum(axis=1)
        ratio = np.maximum(bi_blend / 51708.0, 1e-300)
        phi_F = np.exp(np.sqrt((-0.91725) * np.log(ratio)) + 2.6287)
        phi = (phi_F - 32.0) * 5.0 / 9.0 + 273.15
        L = _gaussian_log_likelihood(phi, data[:, 0], std_data)
        L[has_invalid_flash] = -np.inf
        return L

    elif variable == "freezing":
        W_mat = np.zeros((N_samples, Nc))
        rho_mat = np.zeros((N_samples, Nc))
        bi_freeze = np.zeros((N_samples, Nc))
        for j in range(Nc):
            for k in range(N_samples):
                sp = classes[j][index_n_eta[k, j]]
                W_mat[k, j] = sp.mol_weight
                rho_mat[k, j] = liquid_property("rho", 300.0, sp, pressure)
                bi_freeze[k, j] = sp.Tfz ** (1.0 / 0.05)
        rho_mix = (mol_frac_mat * W_mat).sum(axis=1) / ((mol_frac_mat * W_mat) / rho_mat).sum(axis=1)
        Y_i = np.array([mol_to_mass(mol_frac_mat[k], W_mat[k]) for k in range(N_samples)])
        V_i = rho_mix[:, np.newaxis] * Y_i / rho_mat
        bi_blend = np.maximum(0.0, (V_i * bi_freeze).sum(axis=1))
        phi = bi_blend ** 0.05
        return _gaussian_log_likelihood(phi, data[:, 0], std_data)

    elif variable == "distillation":
        vol_all, T_all = distillation_curve_batch(X_all, classes, index_n_eta, pressure)
        _cache["vol_all"] = vol_all
        _cache["T_all"] = T_all
        _cache["pressure"] = pressure
        _cache["index_n_eta"] = index_n_eta
        _cache["X_all"] = X_all

        T_target = data[:, 0]
        T_interp = np.array([
            _interp_extrap(T_target, vol_all[s], T_all[s]) for s in range(N_samples)
        ])
        L = _gaussian_log_likelihood(T_interp, data[:, 1], std_data)
        # Set -inf for samples where any interpolated value is non-finite
        # (matches MATLAB: L(any(~isfinite(T_interpolated), 2)) = -Inf)
        L[~np.isfinite(T_interp).all(axis=1)] = -np.inf
        return L

    elif variable == "deltaT_dist":
        if "vol_all" in _cache and _cache.get("pressure") == pressure:
            vol_all = _cache["vol_all"]
            T_all = _cache["T_all"]
        else:
            vol_all, T_all = distillation_curve_batch(X_all, classes, index_n_eta, pressure)

        T10 = np.array([_interp_extrap(10.0,  vol_all[s], T_all[s]) for s in range(N_samples)]).ravel()
        T50 = np.array([_interp_extrap(50.0,  vol_all[s], T_all[s]) for s in range(N_samples)]).ravel()
        T90 = np.array([_interp_extrap(90.0,  vol_all[s], T_all[s]) for s in range(N_samples)]).ravel()

        N_exp = len(data[:, 0])
        phi = np.zeros((N_samples, N_exp))
        for i in range(N_exp):
            if data[i, 0] == 5010:
                phi[:, i] = T50 - T10
            elif data[i, 0] == 9010:
                phi[:, i] = T90 - T10

        L = _gaussian_log_likelihood(phi, data[:, 1], std_data)
        L[~np.isfinite(phi).all(axis=1)] = -np.inf
        return L

    else:
        # General temperature-dependent property
        T_exp = data[:, 0]
        phi = mixture_property(variable, T_exp, X_all, classes, index_n_eta, pressure)
        return _gaussian_log_likelihood(phi, data[:, 1], std_data)


# ---------------------------------------------------------------------------
# Cheap likelihood
# ---------------------------------------------------------------------------

def log_likelihood_cheap(
    families: list[str],
    classes: list[list[Species]],
    data: np.ndarray,
    std_data: np.ndarray,
    X_all: np.ndarray,
    n_all: np.ndarray,
    eta_B_star_all: np.ndarray,
    variable: str,
    pressure: float,
    *,
    _cache: dict | None = None,
) -> np.ndarray:
    """
    Lightweight variant of :func:`log_likelihood`.

    For ``'distillation'``, only the first data point is evaluated using an
    early-terminated distillation curve (stops at ``vol1 = data[0, 0]``).
    All other properties are identical to :func:`log_likelihood`.
    """
    if _cache is None:
        _cache = {}

    if variable == "distillation":
        N_samples = X_all.shape[0]
        Nc = X_all.shape[1] + 1
        mol_frac_mat = np.empty((N_samples, Nc), dtype=float)
        mol_frac_mat[:, :Nc - 1] = X_all
        mol_frac_mat[:, Nc - 1] = 1.0 - X_all.sum(axis=1)

        index_n_eta = _build_index_matrix(classes, n_all, eta_B_star_all)

        vol1 = float(data[0, 0])
        vol_all, T_all = distillation_curve_batch_cheap(X_all, classes, index_n_eta, pressure, vol1)

        _cache["vol_all"] = vol_all
        _cache["T_all"] = T_all
        _cache["pressure"] = pressure

        std_data = std_data.ravel()
        T_interp = np.array([
            _interp_extrap(vol1, vol_all[s], T_all[s])[0]
            for s in range(N_samples)
        ])
        L = _gaussian_log_likelihood(T_interp, data[0:1, 1], std_data[0:1])
        L[~np.isfinite(T_interp)] = -np.inf
        return L

    else:
        return log_likelihood(
            families, classes, data, std_data,
            X_all, n_all, eta_B_star_all, variable, pressure,
            _cache=_cache,
        )
