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
The sobol module calculates grouped first-order Sobol' indices based on the
posterior PDF about three sets of input parameters and their impact on
non-lumped thermophysical properties other than the distillation curve:
i) molar fractions, ii) numbers of carbon atoms, and iii) topochemical atom
indices. First-order Sobol' indices are calculated according to the algorithm
proposed in Sobol' and Myshetskaya (2008).

References:
Sobol', I.M. and Myshetskaya, E.E., 2008.
Monte Carlo estimators for small sensitivity indices.
Monte Carlo Methods Appl 13 (5-6), 455-65.

Auxiliary parameters:
N_chains     : number of DE-MC chains
t_convergence: number of iterations at which convergence of the DE-MC run
               is reached according to the R-hat statistic (= maxIterations
               if convergence is not reached)
t_burnin     : number of iterations to be discarded as burn-in
maxIterations: maximum number of iterations for each chain in DE-MC

Inputs:
1) classes              : (1 x numComponents) list of lists of Species —
                          candidate species per surrogate component
2) numComponents        : number of surrogate mixture components  [-]
3) temperature_range    : (1 x 25) array of temperature values ranging from
                          the minimum to the maximum temperature in the
                          dataset of the thermophysical property under
                          consideration
4) chain_reshaped       : ((t_convergence - t_burnin)*N_chains x
                          3*numComponents-1) array of post-burnin posterior
                          samples (molar fractions, nC, eta_B_star_norm)
5) model_output_samples : ((t_convergence - t_burnin)*N_chains x 25) array
                          of model evaluations along temperature_range for
                          each sample in chain_reshaped
6) variable_name        : string denoting the thermophysical property
7) pressure_distillation: pressure at which the distillation curve is
                          computed  [Pa]

Outputs:
1) sobol_idx: (3 x 25) array of grouped first-order Sobol' indices along
              temperature_range
------------------------------------------------------------------------
"""

from __future__ import annotations

import numpy as np

from bayesaf.thermo_transport.hydrocarbons import Species
from bayesaf.thermo_transport.properties import mixture_property
from bayesaf.utilities.find_index import find_index_eta
from bayesaf.distillation.distillation_curve import distillation_curve, distillation_curve_batch


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _build_index_matrix(
    classes: list[list[Species]],
    chain_block: np.ndarray,   # (N, Nc)  — nC values
    eta_block: np.ndarray,      # (N, Nc)  — eta_B_star_norm values
) -> np.ndarray:
    N, Nc = chain_block.shape
    idx = np.empty((N, Nc), dtype=int)
    for j in range(Nc):
        for k in range(N):
            idx[k, j] = find_index_eta(classes[j], int(chain_block[k, j]), float(eta_block[k, j]))
    return idx


def _eval_property(
    variable: str,
    temperature_range: np.ndarray,
    chain_ab: np.ndarray,       # (N, 3*Nc-1)
    classes: list[list[Species]],
    Nc: int,
    pressure: float,
    unit_scale: float = 1.0,
) -> np.ndarray:
    """Evaluate thermophysical property for each sample × temperature."""
    N = chain_ab.shape[0]
    n_T = len(temperature_range)
    X_ab = chain_ab[:, :Nc - 1]
    nC_ab = chain_ab[:, Nc - 1 : 2 * Nc - 1]
    eta_ab = chain_ab[:, 2 * Nc - 1 :]

    idx = _build_index_matrix(classes, nC_ab, eta_ab)
    phi = mixture_property(variable, temperature_range, X_ab, classes, idx, pressure)
    return phi * unit_scale  # (N, n_T)


def _sobol_myshetskaya(
    y_A: np.ndarray,  # (N//2, n_out)
    y_B: np.ndarray,  # (N//2, n_out)
    y_AB: np.ndarray, # (N//2, n_out)
    c: np.ndarray,    # (n_out,)
    Vy: np.ndarray,   # (n_out,)
) -> np.ndarray:
    """First-order Sobol' index per output using Sobol'-Myshetskaya estimator."""
    V = np.abs(np.mean((y_A - c) * (y_AB - y_B), axis=0))
    with np.errstate(divide="ignore", invalid="ignore"):
        S = np.where(Vy > 0, V / Vy, 0.0)
    return S


# ---------------------------------------------------------------------------
# Sobol' indices for temperature-dependent properties
# ---------------------------------------------------------------------------

def sobol_indices(
    classes: list[list[Species]],
    num_components: int,
    temperature_range: np.ndarray,
    chain_reshaped: np.ndarray,
    model_output_samples: np.ndarray,
    variable_name: str,
    pressure_distillation: float,
    seed: int = 0,
) -> np.ndarray:
    """
    Grouped first-order Sobol' indices for a temperature-dependent property.

    Parameters
    ----------
    classes : list of list[Species]
    num_components : int
    temperature_range : ndarray, shape (n_T,)
    chain_reshaped : ndarray, shape (N, 3*Nc-1)
        Post-burnin samples (molar fractions, nC, eta_B_star_norm).
    model_output_samples : ndarray, shape (N, n_T)
        Model evaluations at each temperature for each sample.
    variable_name : str
    pressure_distillation : float
    seed : int

    Returns
    -------
    ndarray, shape (3, n_T)
        Normalised first-order Sobol' indices for the three parameter groups:
        [0] molar fractions, [1] nC, [2] eta_B_star.
    """
    rng = np.random.default_rng(seed)
    Nc = num_components

    # Remove duplicates
    rounded = np.round(chain_reshaped, 6)
    _, uniq_idx = np.unique(rounded, axis=0, return_index=True)
    samples = chain_reshaped[uniq_idx]
    output = model_output_samples[uniq_idx]

    if samples.shape[0] % 2 != 0:
        samples = samples[:-1]
        output = output[:-1]

    N = samples.shape[0]
    perm = rng.permutation(N)
    A = samples[perm[: N // 2]]
    B = samples[perm[N // 2 :]]
    y_A = output[perm[: N // 2]]
    y_B = output[perm[N // 2 :]]

    all_out = np.vstack([y_A, y_B])
    Vy = all_out.var(axis=0)
    c = all_out.mean(axis=0)

    n_T = len(temperature_range)
    sobol = np.zeros((3, n_T))

    # Unit scaling
    scale = 1e6 if variable_name in ("mu", "nu") else (1e-3 if variable_name == "specificHeat" else 1.0)

    for group in range(3):
        A_B = B.copy()
        if group == 0:
            A_B[:, :Nc - 1] = A[:, :Nc - 1]
        elif group == 1:
            A_B[:, Nc - 1 : 2 * Nc - 1] = A[:, Nc - 1 : 2 * Nc - 1]
        else:
            A_B[:, 2 * Nc - 1 :] = A[:, 2 * Nc - 1 :]

        y_AB = _eval_property(variable_name, temperature_range, A_B, classes, Nc, pressure_distillation, scale)
        sobol[group] = _sobol_myshetskaya(y_A, y_B, y_AB, c, Vy)

    # Normalise rows
    col_sum = sobol.sum(axis=0, keepdims=True)
    col_sum[col_sum == 0] = 1.0
    return sobol / col_sum


# ---------------------------------------------------------------------------
# Sobol' indices for the distillation curve
# ---------------------------------------------------------------------------

def sobol_distillation(
    classes: list[list[Species]],
    num_components: int,
    chain_reshaped: np.ndarray,
    model_output_samples: np.ndarray,
    pressure_distillation: float,
    seed: int = 0,
) -> np.ndarray:
    """
    Grouped first-order Sobol' indices for the distillation curve.

    Parameters
    ----------
    model_output_samples : ndarray, shape (N, 101)
        Distillation temperatures at vol-fractions 0…100 for each sample.

    Returns
    -------
    ndarray, shape (3, 101)
        Normalised Sobol' indices at 101 volume-fraction points.
    """
    rng = np.random.default_rng(seed)
    Nc = num_components

    vol_frac_interp = np.concatenate([[0.1], np.linspace(1, 99, 99), [99.9]])

    rounded = np.round(chain_reshaped, 6)
    _, uniq_idx = np.unique(rounded, axis=0, return_index=True)
    samples = chain_reshaped[uniq_idx]
    output = model_output_samples[uniq_idx]   # (N, 101)

    if samples.shape[0] % 2 != 0:
        samples = samples[:-1]
        output = output[:-1]

    N = samples.shape[0]
    perm = rng.permutation(N)
    A = samples[perm[: N // 2]]
    B = samples[perm[N // 2 :]]

    # Interpolate model outputs to vol_frac_interp
    vol_grid = np.arange(101, dtype=float)

    def interp_output(out_block):
        N_b = out_block.shape[0]
        res = np.zeros((N_b, len(vol_frac_interp)))
        for s in range(N_b):
            valid = ~np.isnan(out_block[s])
            res[s] = np.interp(vol_frac_interp, vol_grid[valid], out_block[s][valid])
        return res

    # Output for A and B samples (indices already applied via perm)
    y_A = interp_output(output[perm[: N // 2]])
    y_B = interp_output(output[perm[N // 2 :]])

    all_out = np.vstack([y_A, y_B])
    Vy = all_out.var(axis=0)
    c = all_out.mean(axis=0)

    sobol = np.zeros((3, len(vol_frac_interp)))

    def eval_dist_batch(ab_block):
        N_b = ab_block.shape[0]
        X_ab = ab_block[:, :Nc - 1]
        nC_ab = ab_block[:, Nc - 1 : 2 * Nc - 1]
        eta_ab = ab_block[:, 2 * Nc - 1 :]

        # Pre-compute species indices (fast lookup, not the bottleneck)
        idx_all = np.array([
            [find_index_eta(classes[j], int(nC_ab[s, j]), float(eta_ab[s, j]))
             for j in range(Nc)]
            for s in range(N_b)
        ], dtype=int)

        # Compute all distillation curves in one parallel batch call
        vf_batch, T_batch = distillation_curve_batch(X_ab, classes, idx_all, pressure_distillation)

        res = np.zeros((N_b, len(vol_frac_interp)))
        for s in range(N_b):
            valid = ~np.isnan(T_batch[s])
            res[s] = np.interp(vol_frac_interp, vf_batch[s, valid], T_batch[s, valid])
        return res

    for group in range(3):
        A_B = B.copy()
        if group == 0:
            A_B[:, :Nc - 1] = A[:, :Nc - 1]
        elif group == 1:
            A_B[:, Nc - 1 : 2 * Nc - 1] = A[:, Nc - 1 : 2 * Nc - 1]
        else:
            A_B[:, 2 * Nc - 1 :] = A[:, 2 * Nc - 1 :]

        y_AB = eval_dist_batch(A_B)
        sobol[group] = _sobol_myshetskaya(y_A, y_B, y_AB, c, Vy)

    col_sum = sobol.sum(axis=0, keepdims=True)
    col_sum[col_sum == 0] = 1.0
    return sobol / col_sum
