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
The props_pushforward module provides tools to analyze and visualize results
from the Bayesian inference analysis. It evaluates MAP surrogate properties
and computes pushforward statistics (mean, confidence bands, percentiles) from
the post-burnin chain samples. The module is equivalent to the computation
part of the MATLAB props_pushforward.m function.

Auxiliary parameters:
Nd           : number of experimental measurements for each property  [-]
Np           : number of thermophysical properties targeted during surrogate
               formulation  [-]
N_chains     : number of DE-MC chains
maxIterations: maximum number of iterations for each chain in DE-MC

Inputs:
1) fullData              : (1 x Np) list of (Nd x 3) or (Nd x 4) arrays of
                           experimental measurements
2) families              : (1 x numComponents) list of hydrocarbon family
                           strings
3) classes               : (1 x numComponents) list of lists of Species
4) chain                 : (t_convergence x 3*numComponents-1 x N_chains)
                           full posterior sample array
5) AR                    : (t_convergence x N_chains) acceptance rate array
6) R_hat                 : R-hat statistic array,
                           shape (floor((t_convergence-t_burnin)/2) x
                           3*numComponents-1)
7) chain_reshaped        : ((t_convergence - t_burnin)*N_chains x
                           3*numComponents-1) post-burnin sample array
8) t_convergence         : iteration at which convergence was declared
9) t_burnin              : number of burn-in iterations discarded
10) n_ranges             : (1 x numComponents) carbon-atom ranges per
                           component
11) numComponents        : number of surrogate mixture components  [-]
12) variable_names       : (1 x Np) list of property name strings
13) x_MAP                : (1 x numComponents) MAP molar fractions
14) nc_MAP               : (1 x numComponents) MAP numbers of carbon atoms
15) eta_B_star_MAP       : (1 x numComponents) MAP topochemical indices
16) minT_array           : (1 x Np) minimum temperatures per property
17) maxT_array           : (1 x Np) maximum temperatures per property
18) pressure_distillation: pressure for distillation curve evaluation  [Pa]
19) fuel_name            : legend string for real fuel experimental data
20) band_percentiles     : whether to colour the confidence band with
                           percentiles ('True' or 'False')
21) confidence_width     : width of the confidence interval

Outputs:
1) temperature_range     : temperature grid per property
2) property_MAP          : MAP surrogate property values
3) volumeFraction_MAP    : MAP surrogate distillation volume fractions
4) percentiles_list      : list of percentile arrays
5) mean_model            : posterior mean model evaluations
6) percentiles_model     : model percentile evaluations
7) model_output_samples  : all post-burnin model evaluations
8) volFrac_interp        : interpolated distillation volume fractions
9) sobol_idx             : grouped first-order Sobol' indices
------------------------------------------------------------------------
"""

from __future__ import annotations

from typing import Sequence

import numpy as np

from bayesaf.thermo_transport.hydrocarbons import Species
from bayesaf.thermo_transport.properties import liquid_property, mixture_property
from bayesaf.utilities.find_index import find_index_eta
from bayesaf.utilities.composition import mol_to_mass, mol_weight
from bayesaf.distillation.distillation_curve import (
    distillation_curve,
    distillation_curve_batch,
    init_pool,
    shutdown_pool,
)

# Families whose H count formula is 2n+2
_PARAFFIN_FAMILIES = {"nparaffins", "isoparaffins", "isoparaffins_mono_bis"}
_H_COUNT = {
    "nparaffins": lambda n: 2 * n + 2,
    "isoparaffins": lambda n: 2 * n + 2,
    "isoparaffins_mono_bis": lambda n: 2 * n + 2,
    "cycloparaffins": lambda n: 2 * n,
    "dicycloparaffins": lambda n: 2 * n - 2,
    "alkylbenzenes": lambda n: 2 * n - 6,
    "alkylnaphthalenes": lambda n: 2 * n - 12,
    "cycloaromatics": lambda n: 2 * n - 8,
}

_LUMPED_VARS = {"molWeight", "HC", "DCN", "flash", "freezing", "LHV"}
_TEMP_VARS = {"rho", "mu", "nu", "kappa", "specificHeat", "latentHeat",
              "vaporPressure", "sigma", "distillation"}

_SCALE = {
    "mu": 1e6, "nu": 1e6,
    "specificHeat": 1e-3,
}


# ---------------------------------------------------------------------------
# Helpers for lumped scalar properties
# ---------------------------------------------------------------------------

def _hc_ratio(families: list[str], mol_frac: np.ndarray, nC: np.ndarray) -> float:
    Nc = len(families)
    # Reconstruct full mol_frac if only Nc-1 fractions are provided (last component implicit)
    if len(mol_frac) == Nc - 1:
        mol_frac = np.append(mol_frac, 1.0 - mol_frac.sum())
    n_carbon = round(float(np.dot(mol_frac, nC)), 2)
    n_hydrogen = sum(
        float(mol_frac[l]) * _H_COUNT[families[l]](int(nC[l]))
        for l in range(Nc)
    )
    n_hydrogen = round(n_hydrogen, 2)
    return round(n_hydrogen / n_carbon, 3) if n_carbon > 0 else 0.0


def _dcn(classes: list[list[Species]], mol_frac: np.ndarray, idx: np.ndarray, pressure: float) -> float:
    Nc = len(mol_frac)
    Wl = np.array([classes[l][idx[l]].mol_weight for l in range(Nc)])
    rho_l = np.array([liquid_property("rho", 300.0, classes[l][idx[l]], pressure) for l in range(Nc)])
    DCN_l = np.array([classes[l][idx[l]].DCN for l in range(Nc)])
    numerator = mol_frac * Wl
    y_l = mol_to_mass(mol_frac, Wl)
    rho_mix = numerator.sum() / (numerator / rho_l).sum()
    V_l = rho_mix * y_l / rho_l
    return float(np.dot(V_l, DCN_l))


def _flash_point(classes: list[list[Species]], mol_frac: np.ndarray, idx: np.ndarray) -> float:
    Nc = len(mol_frac)
    Wl = np.array([classes[l][idx[l]].mol_weight for l in range(Nc)])
    Tf_l = np.array([classes[l][idx[l]].Tf for l in range(Nc)])
    Tf_F = (Tf_l - 273.15) * 9.0 / 5.0 + 32.0
    # guard: Tf_F must be > 0 for log to be real
    Tf_F = np.clip(Tf_F, 1e-6, None)
    BI = 51708.0 * np.exp((np.log(Tf_F) - 2.6287) ** 2 / (-0.91725))
    y_l = mol_to_mass(mol_frac, Wl)
    BI_blend = float(np.dot(y_l, BI))
    # guard: argument of outer log must be > 0
    if BI_blend <= 0.0:
        return float("nan")
    ratio = BI_blend / 51708.0
    inner = (-0.91725) * np.log(ratio)
    # guard: inner must be >= 0 for sqrt to be real
    if inner < 0.0:
        return float("nan")
    Tf_F_blend = float(np.exp(inner ** 0.5 + 2.6287))
    return (Tf_F_blend - 32.0) * 5.0 / 9.0 + 273.15


def _freezing_point(classes: list[list[Species]], mol_frac: np.ndarray, idx: np.ndarray, pressure: float) -> float:
    Nc = len(mol_frac)
    Wl = np.array([classes[l][idx[l]].mol_weight for l in range(Nc)])
    rho_l = np.array([liquid_property("rho", 300.0, classes[l][idx[l]], pressure) for l in range(Nc)])
    Tfz_l = np.array([classes[l][idx[l]].Tfz for l in range(Nc)])
    # Use integer exponent 20 = 1/0.05 so that negative Tfz values (extreme
    # MCMC proposals) yield a real result, matching MATLAB's behaviour.
    BI = Tfz_l ** 20
    numerator = mol_frac * Wl
    y_l = mol_to_mass(mol_frac, Wl)
    rho_mix = numerator.sum() / (numerator / rho_l).sum()
    V_l = rho_mix * y_l / rho_l
    BI_blend = float(np.dot(V_l, BI))
    # BI_blend^0.05 = BI_blend^(1/20); take abs to stay real, restore sign
    sign = 1.0 if BI_blend >= 0.0 else -1.0
    return sign * abs(BI_blend) ** 0.05


def _lhv(classes: list[list[Species]], mol_frac: np.ndarray, idx: np.ndarray) -> float:
    Nc = len(mol_frac)
    Wl = np.array([classes[l][idx[l]].mol_weight for l in range(Nc)])
    LHV_l = np.array([classes[l][idx[l]].Hc for l in range(Nc)])
    y_l = mol_to_mass(mol_frac, Wl)
    return float(np.dot(y_l, LHV_l))


def _delta_t_dist(classes: list[list[Species]], mol_frac: np.ndarray, idx: np.ndarray, pressure: float) -> np.ndarray:
    # distillation_curve expects (Nc-1,); mol_frac may be (Nc,) so strip last
    x_in = mol_frac[:-1] if len(mol_frac) == len(classes) else mol_frac
    vf, T = distillation_curve(x_in, classes, idx, pressure)
    valid = ~np.isnan(T)
    T10 = float(np.interp(10.0, vf[valid], T[valid]))
    T50 = float(np.interp(50.0, vf[valid], T[valid]))
    T90 = float(np.interp(90.0, vf[valid], T[valid]))
    return np.array([T50 - T10, T90 - T10])


# ---------------------------------------------------------------------------
# Main function
# ---------------------------------------------------------------------------

def props_pushforward(
    full_data: list[np.ndarray],
    families: list[str],
    classes: list[list[Species]],
    chain_reshaped: np.ndarray,
    num_components: int,
    variable_names: list[str],
    x_MAP: np.ndarray,
    nc_MAP: np.ndarray,
    eta_B_star_MAP: np.ndarray,
    min_T_array: list[float],
    max_T_array: list[float],
    pressure_distillation: float,
    confidence_width: float = 0.9,
) -> dict:
    """
    Compute pushforward statistics for all target properties.

    Parameters
    ----------
    full_data : list of ndarray
        Experimental data per property (columns: T/x, value, std[, pressure]).
    families : list[str]
    classes : list of list[Species]
    chain_reshaped : ndarray, shape (N, 3*Nc-1)
        Post-burnin samples.
    num_components : int
    variable_names : list[str]
    x_MAP, nc_MAP, eta_B_star_MAP : ndarray, shape (Nc,)
    min_T_array, max_T_array : list[float]  (K)
    pressure_distillation : float
    confidence_width : float
        Width of the confidence interval (e.g. 0.9 → 5th–95th percentile band).

    Returns
    -------
    dict with keys:
        temperature_ranges : list[ndarray]   — °C, or vol-frac for distillation
        property_MAP       : list             — MAP value(s)
        volumeFraction_MAP : ndarray | None   — vol-frac grid for distillation
        percentiles_list   : ndarray          — quantile levels
        mean_model         : list[ndarray]
        percentiles_model  : list[ndarray]
        model_output_samples : list[ndarray]
        vol_frac_interp    : ndarray          — interpolation grid for distillation
        sobol_idx          : list             — Sobol' indices per variable
    """
    from bayesaf.gsa.sobol import sobol_indices, sobol_distillation
    from bayesaf.gsa.sobol_lumped import sobol_lumped

    Nc = num_components
    N = chain_reshaped.shape[0]

    # MAP species indices
    idx_MAP = np.array([
        find_index_eta(classes[i], int(nc_MAP[i]), float(eta_B_star_MAP[i]))
        for i in range(Nc)
    ], dtype=int)

    # Confidence interval percentile levels
    lo = (1.0 - confidence_width) / 2.0
    hi = 1.0 - lo
    percentiles_list = np.arange(lo, hi + 0.005, 0.01)

    # Unique-sample cache for expensive evaluations
    rounded = np.round(chain_reshaped, 6)
    _, uniq_rows_idx, uniq_inv = np.unique(rounded, axis=0, return_index=True, return_inverse=True)
    unique_rows = chain_reshaped[uniq_rows_idx]
    N_uniq = unique_rows.shape[0]

    # Locate the MAP sample inside the unique-row table so that all MAP
    # outputs are drawn from the same cache as the CI bands (guaranteeing
    # numerical identity and avoiding any discrepancy between
    # distillation_curve and distillation_curve_batch for the same input).
    _map_row_rounded = np.round(
        np.concatenate([x_MAP[:-1], nc_MAP, eta_B_star_MAP]), 6
    )
    _map_matches = np.where(np.all(rounded == _map_row_rounded, axis=1))[0]
    _uniq_map_idx: int | None = int(uniq_inv[_map_matches[0]]) if len(_map_matches) > 0 else None

    # Species indices for unique rows
    idx_uniq = np.array([
        [find_index_eta(classes[k], int(unique_rows[j, Nc - 1 + k]), float(unique_rows[j, 2 * Nc - 1 + k]))
         for k in range(Nc)]
        for j in range(N_uniq)
    ], dtype=int)

    vol_frac_interp = np.concatenate([[0.1], np.linspace(1, 99, 99), [99.9]])

    # Output containers
    temperature_ranges: list = []
    property_MAP: list = []
    volumeFraction_MAP = None
    mean_model: list = []
    percentiles_model: list = []
    model_output_samples: list = []
    sobol_idx: list = []

    # ---------- Cache distillation outputs (shared between 'distillation' and 'deltaT_dist') ----------
    _dist_cache: dict = {}  # {j_uniq: T_interp_array_degC}
    _delta_cache: dict = {}  # {j_uniq: [T50-T10, T90-T10]}

    # ── Pre-compute all unique distillation curves in parallel ────────────────
    # If the run includes either 'distillation' or 'deltaT_dist', compute all
    # N_uniq curves in a single parallelised batch call before entering the
    # property loop.  The results are stored in _dist_cache / _delta_cache so
    # the per-property sections below simply read from memory.
    needs_dist = any(v in ("distillation", "deltaT_dist") for v in variable_names)
    if needs_dist:
        init_pool(classes=classes)
        print(f"  Pre-computing distillation curves for {N_uniq} unique samples "
              f"(parallel)...")
        X_uniq_part = unique_rows[:, :Nc - 1]
        vf_batch, T_batch = distillation_curve_batch(
            X_uniq_part, classes, idx_uniq, pressure_distillation
        )
        for j in range(N_uniq):
            vf_j, T_j = vf_batch[j], T_batch[j]
            valid = ~np.isnan(T_j)
            T_interp = np.interp(vol_frac_interp, vf_j[valid], T_j[valid]) - 273.15
            _dist_cache[j] = T_interp
            T10 = float(np.interp(10.0, vf_j[valid], T_j[valid]))
            T50 = float(np.interp(50.0, vf_j[valid], T_j[valid]))
            T90 = float(np.interp(90.0, vf_j[valid], T_j[valid]))
            _delta_cache[j] = np.array([T50 - T10, T90 - T10])

    try:
        for h, var in enumerate(variable_names):
            print(f"  [{h+1}/{len(variable_names)}] Computing pushforward for '{var}'...")

            if var not in _LUMPED_VARS:
                # ---------- Temperature-dependent property ----------

                minT = min_T_array[h]
                maxT = max_T_array[h]
                if minT != maxT:
                    T_range_K = np.linspace(minT, maxT, 25)
                else:
                    T_range_K = np.linspace(233.15, 313.15, 25)
                T_range_C = T_range_K - 273.15

                if var not in ("distillation", "deltaT_dist"):
                    # MAP property
                    scale = _SCALE.get(var, 1.0)
                    prop_map_arr = np.array([
                        scale * float(mixture_property(var, np.array([T]), x_MAP[: Nc - 1].reshape(1, -1),
                                                       classes, idx_MAP.reshape(1, -1), pressure_distillation)[0, 0])
                        for T in T_range_K
                    ])
                    property_MAP.append(prop_map_arr)
                    temperature_ranges.append(T_range_C)

                    # Unique-sample evaluations: shape (N_uniq, n_T)
                    scale = _SCALE.get(var, 1.0)
                    X_uniq = unique_rows[:, : Nc - 1]
                    uniq_out = mixture_property(var, T_range_K, X_uniq, classes, idx_uniq, pressure_distillation) * scale

                    # Map back to all samples
                    out_all = uniq_out[uniq_inv]  # (N, n_T)
                    model_output_samples.append(out_all)
                    mean_model.append(out_all.mean(axis=0))
                    pcts = np.quantile(out_all, percentiles_list, axis=0).T  # (n_T, n_pct)
                    percentiles_model.append(pcts)

                    # Sobol indices
                    si = sobol_indices(classes, Nc, T_range_K, chain_reshaped, out_all,
                                       var, pressure_distillation)
                    sobol_idx.append(si)

                elif var == "distillation":
                    # Cache already populated in the pre-computation block above.
                    # Map to all samples: (N, 101)
                    out_all = np.array([_dist_cache[uniq_inv[s]] for s in range(N)])
                    model_output_samples.append(out_all)
                    mean_model.append(out_all.mean(axis=0))
                    pcts = np.quantile(out_all, percentiles_list, axis=0).T
                    percentiles_model.append(pcts)

                    # MAP distillation curve — read from the batch cache so it is
                    # numerically identical to the CI bands (same code path).
                    if _uniq_map_idx is not None:
                        property_MAP.append(_dist_cache[_uniq_map_idx])   # already °C
                        volumeFraction_MAP = vol_frac_interp
                    else:
                        # Fallback (MAP not found in chain — should not happen)
                        vf_map, T_map = distillation_curve(x_MAP[:-1], classes, idx_MAP, pressure_distillation)
                        volumeFraction_MAP = vf_map
                        property_MAP.append(T_map - 273.15)
                    temperature_ranges.append(vol_frac_interp)

                    si = sobol_distillation(classes, Nc, chain_reshaped, out_all + 273.15,
                                            pressure_distillation)
                    sobol_idx.append(si)

                elif var == "deltaT_dist":
                    # Cache already populated in the pre-computation block above.
                    # Map to all samples: (N, 2)
                    out_all = np.array([_delta_cache[uniq_inv[s]] for s in range(N)])
                    model_output_samples.append(out_all)
                    mean_model.append(out_all.mean(axis=0))
                    pcts = np.quantile(out_all, percentiles_list, axis=0).T
                    percentiles_model.append(pcts)

                    # MAP delta T values — read from the batch cache for consistency.
                    if _uniq_map_idx is not None:
                        property_MAP.append(_delta_cache[_uniq_map_idx])  # [T50-T10, T90-T10] in K
                    else:
                        vf_map, T_map = distillation_curve(x_MAP[:-1], classes, idx_MAP, pressure_distillation)
                        valid = ~np.isnan(T_map)
                        T10_map = float(np.interp(10.0, vf_map[valid], T_map[valid]))
                        T50_map = float(np.interp(50.0, vf_map[valid], T_map[valid]))
                        T90_map = float(np.interp(90.0, vf_map[valid], T_map[valid]))
                        property_MAP.append(np.array([T50_map - T10_map, T90_map - T10_map]))
                    temperature_ranges.append(np.array([0.0, 1.0]))  # placeholder

                    si = sobol_lumped(families, classes, Nc, chain_reshaped, out_all,
                                      "deltaT_dist", pressure_distillation)
                    sobol_idx.append(si)

            else:
                # ---------- Lumped property ----------
                temperature_ranges.append(None)

                # Evaluate MAP property
                if var == "molWeight":
                    Wl_map = np.array([classes[l][idx_MAP[l]].mol_weight for l in range(Nc)])
                    prop_val = mol_weight(x_MAP, Wl_map)
                elif var == "HC":
                    prop_val = _hc_ratio(families, x_MAP, nc_MAP)
                elif var == "DCN":
                    prop_val = _dcn(classes, x_MAP, idx_MAP, 101325.0)
                elif var == "flash":
                    prop_val = _flash_point(classes, x_MAP, idx_MAP)
                elif var == "freezing":
                    prop_val = _freezing_point(classes, x_MAP, idx_MAP, 101325.0)
                elif var == "LHV":
                    prop_val = _lhv(classes, x_MAP, idx_MAP)
                else:
                    prop_val = 0.0
                property_MAP.append(prop_val)

                # Evaluate for all chain samples (unique-row cache via uniq_inv)
                out_all = np.zeros(N)
                prev_row = None
                for j in range(N):
                    row = chain_reshaped[j]
                    if prev_row is not None and np.array_equal(row, prev_row):
                        out_all[j] = out_all[j - 1]
                    else:
                        x_part = row[: Nc - 1]
                        x_full = np.append(x_part, 1.0 - x_part.sum())
                        nC_j = row[Nc - 1 : 2 * Nc - 1]
                        idx_j = np.array([
                            find_index_eta(classes[k], int(row[Nc - 1 + k]), float(row[2 * Nc - 1 + k]))
                            for k in range(Nc)
                        ], dtype=int)
                        if var == "molWeight":
                            Wl = np.array([classes[l][idx_j[l]].mol_weight for l in range(Nc)])
                            out_all[j] = mol_weight(x_full, Wl)
                        elif var == "HC":
                            out_all[j] = _hc_ratio(families, x_full, nC_j)
                        elif var == "DCN":
                            out_all[j] = _dcn(classes, x_full, idx_j, 101325.0)
                        elif var == "flash":
                            out_all[j] = _flash_point(classes, x_full, idx_j)
                        elif var == "freezing":
                            out_all[j] = _freezing_point(classes, x_full, idx_j, 101325.0)
                        elif var == "LHV":
                            out_all[j] = _lhv(classes, x_full, idx_j)
                    prev_row = row

                model_output_samples.append(out_all)
                mean_model.append(float(out_all.mean()))
                pcts = np.quantile(out_all, percentiles_list)
                percentiles_model.append(pcts)

                si = sobol_lumped(families, classes, Nc, chain_reshaped, out_all, var, pressure_distillation)
                sobol_idx.append(si)

        return {
            "temperature_ranges": temperature_ranges,
            "property_MAP": property_MAP,
            "volumeFraction_MAP": volumeFraction_MAP,
            "percentiles_list": percentiles_list,
            "mean_model": mean_model,
            "percentiles_model": percentiles_model,
            "model_output_samples": model_output_samples,
            "vol_frac_interp": vol_frac_interp,
            "sobol_idx": sobol_idx,
        }
    finally:
        if needs_dist:
            shutdown_pool()
