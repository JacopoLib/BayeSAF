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
The sobol_lumped module calculates grouped first-order Sobol' indices based
on the posterior PDF about three sets of input parameters and their impact on
lumped mixture properties (molecular weight, hydrogen-to-carbon ratio, and
derived cetane number): i) molar fractions, ii) numbers of carbon atoms, and
iii) topochemical atom indices. First-order Sobol' indices are calculated
according to the algorithm proposed in Sobol' and Myshetskaya (2008).

References:
Sobol', I.M. and Myshetskaya, E.E., 2008.
Monte Carlo estimators for small sensitivity indices.
Monte Carlo Methods Appl 13 (5-6), 455-65.

Hu, J. and Burns, A., 1970.
Index predicts cloud, pour and flash points of distillates fuel blends.
Oil Gas J. 68(45):66.

Riazi, R., 2005.
Characterization and Properties of Petroleum Fractions.
ASTM International, Philadelphia.

Auxiliary parameters:
N_chains     : number of DE-MC chains
t_convergence: number of iterations at which convergence of the DE-MC run
               is reached according to the R-hat statistic (= maxIterations
               if convergence is not reached)
t_burnin     : number of iterations to be discarded as burn-in
maxIterations: maximum number of iterations for each chain in DE-MC

Inputs:
1) families             : (1 x numComponents) list of strings denoting the
                          hydrocarbon family of each surrogate component
2) classes              : (1 x numComponents) list of lists of Species —
                          candidate species per surrogate component
3) numComponents        : number of surrogate mixture components  [-]
4) chain_reshaped       : ((t_convergence - t_burnin)*N_chains x
                          3*numComponents-1) array of post-burnin posterior
                          samples (molar fractions, nC, eta_B_star_norm)
5) model_output_samples : ((t_convergence - t_burnin)*N_chains x 1) array
                          of lumped property values for each sample
6) variable_name        : string denoting the lumped property under
                          consideration
7) pressure_distillation: pressure at which the distillation curve is
                          computed  [Pa]

Outputs:
1) sobol_idx: (3 x 1) array of grouped first-order Sobol' indices
------------------------------------------------------------------------
"""

from __future__ import annotations

import numpy as np

from bayesaf.thermo_transport.hydrocarbons import Species
from bayesaf.thermo_transport.properties import liquid_property
from bayesaf.utilities.find_index import find_index_eta
from bayesaf.utilities.composition import mol_to_mass, mol_weight
from bayesaf.distillation.distillation_curve import distillation_curve, distillation_curve_batch

# H-count formulas per family
_H_FORMULA: dict[str, object] = {
    "nparaffins": lambda n: 2 * n + 2,
    "isoparaffins": lambda n: 2 * n + 2,
    "isoparaffins_mono_bis": lambda n: 2 * n + 2,
    "cycloparaffins": lambda n: 2 * n,
    "dicycloparaffins": lambda n: 2 * n - 2,
    "alkylbenzenes": lambda n: 2 * n - 6,
    "alkylnaphtalenes": lambda n: 2 * n - 12,
    "cycloaromatics": lambda n: 2 * n - 8,
}


def _eval_lumped(
    variable: str,
    families: list[str],
    classes: list[list[Species]],
    mol_frac: np.ndarray,        # (Nc,)
    nC: np.ndarray,              # (Nc,)  integer values
    idx: np.ndarray,             # (Nc,)  species indices
    pressure: float,
) -> float | np.ndarray:
    """Evaluate a lumped property for one sample."""
    Nc = len(mol_frac)
    Wl = np.array([classes[l][idx[l]].mol_weight for l in range(Nc)])

    if variable == "molWeight":
        return mol_weight(mol_frac, Wl)

    elif variable == "HC":
        nC_arr = np.asarray(nC, dtype=float)
        n_carbon = round(float(np.dot(mol_frac, nC_arr)), 2)
        n_hydrogen = sum(
            float(mol_frac[l]) * _H_FORMULA[families[l]](int(nC_arr[l]))
            for l in range(Nc)
        )
        n_hydrogen = round(n_hydrogen, 2)
        return round(n_hydrogen / n_carbon, 3) if n_carbon > 0 else 0.0

    elif variable == "DCN":
        rho_l = np.array([
            liquid_property("rho", 300.0, classes[l][idx[l]], 101325.0)
            for l in range(Nc)
        ], dtype=float)
        DCN_l = np.array([classes[l][idx[l]].DCN for l in range(Nc)])
        numerator = mol_frac * Wl
        y_l = mol_to_mass(mol_frac, Wl)
        rho_mix = numerator.sum() / (numerator / rho_l).sum()
        V_l = rho_mix * y_l / rho_l
        return float(np.dot(V_l, DCN_l))

    elif variable == "flash":
        Tf_l = np.array([classes[l][idx[l]].Tf for l in range(Nc)])
        Tf_F = (Tf_l - 273.15) * 9.0 / 5.0 + 32.0
        Tf_F = np.clip(Tf_F, 1e-6, None)
        BI_flash = 51708.0 * np.exp((np.log(Tf_F) - 2.6287) ** 2 / (-0.91725))
        y_l = mol_to_mass(mol_frac, Wl)
        BI_blend = float(np.dot(y_l, BI_flash))
        if BI_blend <= 0.0:
            return float("nan")
        ratio = BI_blend / 51708.0
        inner = (-0.91725) * np.log(ratio)
        if inner < 0.0:
            return float("nan")
        y_AB_F = float(np.exp(inner ** 0.5 + 2.6287))
        return (y_AB_F - 32.0) * 5.0 / 9.0 + 273.15

    elif variable == "freezing":
        rho_l = np.array([
            liquid_property("rho", 300.0, classes[l][idx[l]], 101325.0)
            for l in range(Nc)
        ], dtype=float)
        Tfz_l = np.array([classes[l][idx[l]].Tfz for l in range(Nc)])
        # Integer exponent 20 = 1/0.05 keeps result real for negative Tfz values
        BI_fz = Tfz_l ** 20
        numerator = mol_frac * Wl
        y_l = mol_to_mass(mol_frac, Wl)
        rho_mix = numerator.sum() / (numerator / rho_l).sum()
        V_l = rho_mix * y_l / rho_l
        BI_blend = float(np.dot(V_l, BI_fz))
        sign = 1.0 if BI_blend >= 0.0 else -1.0
        return sign * abs(BI_blend) ** 0.05

    elif variable == "LHV":
        LHV_l = np.array([classes[l][idx[l]].Hc for l in range(Nc)])
        y_l = mol_to_mass(mol_frac, Wl)
        return float(np.dot(y_l, LHV_l))

    elif variable == "deltaT_dist":
        # distillation_curve expects (Nc-1,); mol_frac is (Nc,) so strip last
        vf, T = distillation_curve(mol_frac[:-1], classes, idx, pressure)
        valid = ~np.isnan(T)
        T10 = float(np.interp(10.0, vf[valid], T[valid]))
        T50 = float(np.interp(50.0, vf[valid], T[valid]))
        T90 = float(np.interp(90.0, vf[valid], T[valid]))
        return np.array([T50 - T10, T90 - T10])

    raise ValueError(f"Unknown variable: {variable!r}")


def sobol_lumped(
    families: list[str],
    classes: list[list[Species]],
    num_components: int,
    chain_reshaped: np.ndarray,
    model_output_samples: np.ndarray,
    variable_name: str,
    pressure_distillation: float,
    seed: int = 0,
) -> np.ndarray:
    """
    Grouped first-order Sobol' indices for lumped mixture properties.

    Parameters
    ----------
    families : list[str]
    classes : list of list[Species]
    num_components : int
    chain_reshaped : ndarray, shape (N, 3*Nc-1)
        Post-burnin samples (molar fractions, nC, eta_B_star_norm).
    model_output_samples : ndarray, shape (N,) or (N, 2) for deltaT_dist
        Pre-computed model evaluations.
    variable_name : str
        One of 'molWeight', 'HC', 'DCN', 'flash', 'freezing', 'LHV', 'deltaT_dist'.
    pressure_distillation : float
    seed : int

    Returns
    -------
    ndarray, shape (3,) or (3, 2) for deltaT_dist
        Normalised first-order Sobol' indices for the three parameter groups:
        [0] molar fractions, [1] nC, [2] eta_B_star.
    """
    rng = np.random.default_rng(seed)
    Nc = num_components
    is_delta = variable_name == "deltaT_dist"

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

    all_out = np.concatenate([y_A, y_B], axis=0)
    Vy = np.var(all_out, axis=0)
    c = np.mean(all_out, axis=0)

    out_shape = (3, 2) if is_delta else (3,)
    sobol = np.zeros(out_shape)

    def _eval_batch(ab_block: np.ndarray) -> np.ndarray:
        N_b = ab_block.shape[0]
        X_ab  = ab_block[:, : Nc - 1]
        nC_ab = ab_block[:, Nc - 1 : 2 * Nc - 1]
        eta_ab = ab_block[:, 2 * Nc - 1 :]

        # Fast path for deltaT_dist: use the parallelised batch distillation
        # call instead of N sequential distillation_curve() calls — same
        # approach as sobol_distillation, giving an N-fold speed-up.
        if is_delta:
            idx_all = np.array([
                [find_index_eta(classes[j], int(nC_ab[s, j]), float(eta_ab[s, j]))
                 for j in range(Nc)]
                for s in range(N_b)
            ], dtype=int)
            vf_batch, T_batch = distillation_curve_batch(
                X_ab, classes, idx_all, pressure_distillation
            )
            out = np.full((N_b, 2), np.nan)
            for s in range(N_b):
                valid = ~np.isnan(T_batch[s])
                if valid.sum() >= 2:
                    vf_v, T_v = vf_batch[s][valid], T_batch[s][valid]
                    T10 = float(np.interp(10.0, vf_v, T_v))
                    T50 = float(np.interp(50.0, vf_v, T_v))
                    T90 = float(np.interp(90.0, vf_v, T_v))
                    out[s] = [T50 - T10, T90 - T10]
            return out

        # General path for all other lumped properties (sequential, cheap).
        out = np.full(N_b, np.nan)
        for s in range(N_b):
            idx = np.array([
                find_index_eta(classes[j], int(nC_ab[s, j]), float(eta_ab[s, j]))
                for j in range(Nc)
            ], dtype=int)
            mol_frac = np.append(X_ab[s], 1.0 - X_ab[s].sum())
            try:
                out[s] = float(_eval_lumped(
                    variable_name, families, classes, mol_frac,
                    nC_ab[s], idx, pressure_distillation
                ))
            except Exception:
                pass   # leave row as NaN; handled by Sobol' estimator via Vy
        return out

    def _sobol_myshetskaya(y_ab: np.ndarray) -> np.ndarray:
        diff = (y_A - c) * (y_ab - y_B)
        V = np.abs(np.mean(diff, axis=0))
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.where(Vy > 0, V / Vy, 0.0)

    for group in range(3):
        A_B = B.copy()
        if group == 0:
            A_B[:, : Nc - 1] = A[:, : Nc - 1]
        elif group == 1:
            A_B[:, Nc - 1 : 2 * Nc - 1] = A[:, Nc - 1 : 2 * Nc - 1]
        else:
            A_B[:, 2 * Nc - 1 :] = A[:, 2 * Nc - 1 :]

        y_AB = _eval_batch(A_B)
        sobol[group] = _sobol_myshetskaya(y_AB)

    # Normalise
    col_sum = sobol.sum(axis=0, keepdims=True)
    col_sum[col_sum == 0] = 1.0
    return sobol / col_sum
