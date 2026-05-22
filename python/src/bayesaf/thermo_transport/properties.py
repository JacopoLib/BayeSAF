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
The properties module calculates the value of a liquid-phase thermophysical
property of a pure chemical species or of a mixture at given temperature and
pressure. All formulae are identical to the MATLAB originals; only the calling
convention has been modernised (NumPy arrays instead of separate scalar
coefficients). Supported properties include: liquid density (rho), dynamic
viscosity (mu), kinematic viscosity (nu), thermal conductivity (kappa),
isobaric specific heat capacity (specificHeat), latent heat of vaporization
(latentHeat), vapour pressure (vaporPressure), and surface tension (sigma).

Inputs:
1) property : string denoting the thermophysical property under consideration
2) T        : temperature  [K]
3) W        : molecular weight  [g/mol]
4) Tcr      : critical temperature  [K]
5) Arho     : A-coefficient for liquid density (Yaws')
6) Brho     : B-coefficient for liquid density (Yaws')
7) Crho     : C-coefficient for liquid density (Yaws')
8) Amu      : A-coefficient for liquid dynamic viscosity (Yaws')
9) Bmu      : B-coefficient for liquid dynamic viscosity (Yaws')
10) Cmu     : C-coefficient for liquid dynamic viscosity (Yaws')
11) Dmu     : D-coefficient for liquid dynamic viscosity (Yaws')
12) Ak      : A-coefficient for liquid thermal conductivity (Yaws')
13) Bk      : B-coefficient for liquid thermal conductivity (Yaws')
14) Ck      : C-coefficient for liquid thermal conductivity (Yaws')
15) Ac      : A-coefficient for liquid specific heat capacity (Yaws')
16) Bc      : B-coefficient for liquid specific heat capacity (Yaws')
17) Cc      : C-coefficient for liquid specific heat capacity (Yaws')
18) Dc      : D-coefficient for liquid specific heat capacity (Yaws')
19) Ah      : A-coefficient for latent heat of vaporization (Yaws')
20) Bh      : B-coefficient for latent heat of vaporization (Yaws')
21) Ap      : A-coefficient for vapour pressure (Yaws')
22) Bp      : B-coefficient for vapour pressure (Yaws')
23) Cp      : C-coefficient for vapour pressure (Yaws')
24) Dp      : D-coefficient for vapour pressure (Yaws')
25) Ep      : E-coefficient for vapour pressure (Yaws')
26) Asigma  : A-coefficient for liquid surface tension (Yaws')
27) Bsigma  : B-coefficient for liquid surface tension (Yaws')
28) pressure: pressure  [Pa]

Outputs:
1) phi: value of the thermophysical property under consideration
------------------------------------------------------------------------
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

from .hydrocarbons import Species


# ---------------------------------------------------------------------------
# Single-species properties
# ---------------------------------------------------------------------------

def liquid_property(
    prop: str,
    T: float | np.ndarray,
    sp: Species,
    pressure: float = 101325.0,
) -> float | np.ndarray:
    """
    Liquid-phase thermophysical property of a pure species at temperature *T*.

    Parameters
    ----------
    prop : str
        One of: ``'rho'``, ``'mu'``, ``'nu'``, ``'kappa'``,
        ``'specificHeat'``, ``'latentHeat'``, ``'vaporPressure'``,
        ``'boilingTemperature'``, ``'sigma'``.
    T : float or ndarray
        Temperature [K].  Not used for ``'boilingTemperature'``.
    sp : Species
        Species data (Yaws' coefficients + critical temperature).
    pressure : float
        Pressure [Pa].  Relevant for ``'boilingTemperature'`` and
        ``'vaporPressure'``.

    Returns
    -------
    float or ndarray
        Property value in SI units (see below).

    Units
    -----
    rho            kg/m³
    mu             Pa·s
    nu             m²/s
    kappa          W/(m·K)
    specificHeat   J/(kg·K)
    latentHeat     J/kg
    vaporPressure  Pa
    boilingTemperature  K
    sigma          N/m
    """
    T = np.asarray(T, dtype=float)

    Arho, Brho, Crho = sp.coeff_rho
    Amu, Bmu, Cmu, Dmu = sp.coeff_mu
    Ak, Bk, Ck = sp.coeff_k
    Ac, Bc, Cc, Dc = sp.coeff_cl
    Ah, Bh = sp.coeff_hv
    Ap, Bp, Cp, Dp, Ep = sp.coeff_psat
    Asigma, Bsigma = sp.coeff_sigma
    W = sp.mol_weight
    Tc = sp.Tc

    if prop == "rho":
        return 1000.0 * Arho * Brho ** (-(1.0 - T / Tc) ** Crho)

    elif prop == "mu":
        return 0.001 * 10.0 ** (Amu + Bmu / T + Cmu * T + Dmu * T**2)

    elif prop == "nu":
        rho = liquid_property("rho", T, sp, pressure)
        mu = liquid_property("mu", T, sp, pressure)
        return mu / rho

    elif prop == "kappa":
        return Ak + Bk * T + Ck * T**2

    elif prop == "specificHeat":
        return 1000.0 / W * (Ac + Bc * T + Cc * T**2 + Dc * T**3)

    elif prop == "latentHeat":
        return (Ah * (1.0 - T / Tc) ** Bh) / W * 1e6

    elif prop == "vaporPressure":
        return 10.0 ** (Ap + Bp / T + Cp * np.log10(T) + Dp * T + Ep * T**2) * 133.322387415

    elif prop == "boilingTemperature":
        def psat_minus_p(T_val: float) -> float:
            psat = 10.0 ** (Ap + Bp / T_val + Cp * np.log10(T_val) + Dp * T_val + Ep * T_val**2) * 133.322387415
            return psat - pressure

        return float(brentq(psat_minus_p, 200.0, 1200.0, xtol=1e-6))

    elif prop == "sigma":
        return (Asigma * (1.0 - T / Tc) ** Bsigma) / 1e3

    else:
        raise ValueError(f"Unknown thermophysical property: '{prop}'")


# ---------------------------------------------------------------------------
# Mixture properties (vectorised over N_samples × N_exp)
# ---------------------------------------------------------------------------

def mixture_property(
    prop: str,
    T: np.ndarray,
    X_all: np.ndarray,
    classes: list[list[Species]],
    index_n_eta_all: np.ndarray,
    pressure: float = 101325.0,
) -> np.ndarray:
    """
    Vectorised liquid-mixture property for multiple MCMC samples.

    Parameters
    ----------
    prop : str
        Property name (same options as :func:`liquid_property`).
    T : ndarray, shape (N_exp,)
        Experimental temperatures [K].
    X_all : ndarray, shape (N_samples, Nc-1)
        Molar fractions (last component computed automatically).
    classes : list of list[Species]
        ``classes[j]`` is the species list for component *j*.
    index_n_eta_all : ndarray of int, shape (N_samples, Nc)
        Species index in ``classes[j]`` for each sample/component.
    pressure : float
        Pressure [Pa].

    Returns
    -------
    ndarray, shape (N_samples, N_exp)
        Property values.
    """
    N_samples, Nc_m1 = X_all.shape
    Nc = Nc_m1 + 1
    N_exp = len(T)
    T = np.asarray(T, dtype=float)

    # Full molar fraction matrix [N_samples × Nc]
    X_full = np.empty((N_samples, Nc), dtype=float)
    X_full[:, :Nc_m1] = X_all
    X_full[:, Nc_m1] = 1.0 - X_all.sum(axis=1)

    if prop == "rho":
        # mixture density: sum(Xi*Wi) / sum(Xi*Wi/rho_i)
        numerator = np.zeros((N_samples, N_exp, Nc))
        rho_mat = np.zeros((N_samples, N_exp, Nc))
        for i in range(Nc):
            for j in range(N_samples):
                sp = classes[i][index_n_eta_all[j, i]]
                rho_mat[j, :, i] = liquid_property("rho", T, sp, pressure)
                numerator[j, :, i] = X_full[j, i] * sp.mol_weight
        return numerator.sum(axis=2) / (numerator / rho_mat).sum(axis=2)

    elif prop == "mu":
        # log mixing rule: ln(mu_mix) = sum(Xi * ln(mu_i))
        mu_log = np.zeros((N_samples, N_exp, Nc))
        for i in range(Nc):
            for j in range(N_samples):
                sp = classes[i][index_n_eta_all[j, i]]
                mu_log[j, :, i] = X_full[j, i] * np.log(
                    liquid_property("mu", T, sp, pressure)
                )
        return np.exp(mu_log.sum(axis=2))

    elif prop == "nu":
        rho = mixture_property("rho", T, X_all, classes, index_n_eta_all, pressure)
        mu = mixture_property("mu", T, X_all, classes, index_n_eta_all, pressure)
        return mu / rho

    elif prop in ("kappa", "specificHeat", "latentHeat", "vaporPressure"):
        # mass-weighted average
        numerator = np.zeros((N_samples, N_exp, Nc))
        denominator = np.zeros((N_samples, Nc))
        for i in range(Nc):
            for j in range(N_samples):
                sp = classes[i][index_n_eta_all[j, i]]
                numerator[j, :, i] = X_full[j, i] * sp.mol_weight * liquid_property(prop, T, sp, pressure)
                denominator[j, i] = X_full[j, i] * sp.mol_weight
        return numerator.sum(axis=2) / denominator.sum(axis=1, keepdims=True)

    elif prop == "sigma":
        sigma_mat = np.zeros((N_samples, N_exp, Nc))
        Xs = np.zeros((N_samples, N_exp, Nc))
        for i in range(Nc):
            for j in range(N_samples):
                sp = classes[i][index_n_eta_all[j, i]]
                sigma_mat[j, :, i] = liquid_property("sigma", T, sp, pressure)
                Xs[j, :, i] = X_full[j, i]
        return (sigma_mat * Xs).sum(axis=2)

    else:
        raise ValueError(f"Unknown thermophysical property: '{prop}'")
