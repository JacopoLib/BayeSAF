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
The demc module adopts the differential evolution Markov chain (DE-MC)
algorithm proposed by ter Braak (2006) to explore and sample from the
posterior probability density function (PDF). The implementation is adapted
from the DiffeRential Evolution Adaptive Metropolis (DREAM) toolbox developed
by Vrugt et al. (2008, 2009, 2016). It supports parallel tempering, delayed
acceptance, snooker jumps, blocked updates, and outlier-chain detection.
Convergence is monitored via the R-hat Gelman-Rubin statistic.

References:
ter Braak, C.J.F. (2006). A Markov chain Monte Carlo version of the genetic
algorithm differential evolution: easy Bayesian computing for real parameter
spaces. Stat. Comput. 16, 239-249.

Vrugt, J.A. et al. (2008, 2009, 2016). DREAM toolbox.

Auxiliary parameters:
Nc: number of surrogate mixture components  [-]

Inputs:
1) families             : (1 x Nc) list of strings denoting the hydrocarbon
                          family of each surrogate component
2) Posterior_PDF        : callable representing the vectorised posterior PDF
3) Posterior_PDF_cheap  : callable representing a cheap surrogate of the
                          posterior PDF (used for delayed acceptance)
4) classes              : (1 x Nc) list of lists of Species — candidate
                          species per surrogate component
5) LowerBound_x         : (1 x Nc) array with lower bounds for molar fractions
6) UpperBound_x         : (1 x Nc) array with upper bounds for molar fractions
7) n_ranges             : (1 x Nc) list of carbon-atom ranges per component
8) LowerBound_eta_B_star: (1 x Nc) array with lower bounds for normalised
                          topochemical atom indices
9) UpperBound_eta_B_star: (1 x Nc) array with upper bounds for normalised
                          topochemical atom indices
10) maxIterations       : maximum number of iterations for each chain
12) t_burnin            : number of iterations to discard as burn-in
13) t_check             : iteration from which to start the R-hat check
14) scaling_factor_X    : jump-rate scaling for molar fractions
15) scaling_factor_nc   : jump-rate scaling for numbers of carbon atoms
16) scaling_factor_eta  : jump-rate scaling for topochemical atom indices
17) noise_X             : noise parameter for molar-fraction proposals
18) noise_nc            : noise parameter for carbon-atom proposals
19) noise_eta           : noise parameter for topochemical-index proposals
20) N_chains            : number of chains
21) outlier_method      : outlier-chain detection method ('iqr' or 'mad')
22) R_hat_threshold     : R-hat convergence threshold
23) PT_switch           : use parallel tempering ('True') or not ('False')
24) beta_min            : minimum inverse temperature for parallel tempering
25) T_ladder            : temperature ladder form ('linear' or 'geometric')
26) swap_freq           : parallel tempering swap frequency
27) n_pairs             : maximum number of chain pairs used to propose the new sample
28) p_gibbs             : probability of a Gibbs move
29) p_snooker           : probability of a snooker jump

Outputs:
1) x            : (t_convergence x 3*Nc-1 x N_chains) array of posterior
                  samples (molar fractions, nC, eta_B_star_norm)
2) p_x          : (t_convergence x N_chains) array of log-posterior values
3) x_nb         : post-burnin samples, shape ((t_convergence-t_burnin)*N_chains
                  x 3*Nc-1)
4) p_x_nb       : post-burnin log-posterior values
5) AR            : (t_convergence x N_chains) acceptance rate array
6) R_hat         : R-hat statistic array, shape
                   (floor((t_convergence-t_burnin)/2) x 3*Nc-1)
7) t_convergence : iteration at which convergence was declared (or
                   maxIterations if convergence was not reached)
------------------------------------------------------------------------
"""

from __future__ import annotations

import sys
from typing import Callable, Sequence

import numpy as np

from bayesaf.thermo_transport.hydrocarbons import Species
from bayesaf.distillation.distillation_curve import init_pool, shutdown_pool
from bayesaf.mcmc.layout import (  # noqa: F401  (decoders re-exported for compatibility)
    ParameterLayout,
    SurrogateLayout,
    _decode_eta,
    _decode_nc,
    _triangle_fold,
)

# Type alias for vectorised posterior callable
PostFn = Callable[[np.ndarray, np.ndarray, np.ndarray], np.ndarray]


# ---------------------------------------------------------------------------
# Marsaglia-Tsang (2000) normal ziggurat — matches MATLAB's randn(MT19937)
# ---------------------------------------------------------------------------
# Reference: G. Marsaglia & W.W. Tsang, "The Ziggurat Method for Generating
#   Random Variables", J. Stat. Software 5(8), 2000.
# Parameters (N=128 strips, r=3.442619855899, v=9.91256303526217e-3) are
# identical to MATLAB's internal randn implementation.
#
# MT integer consumption per sample:
#   fast path (~97.5 %):  1 raw 32-bit MT integer   ← matches MATLAB randn
#   slow path tail:       1 int + 2×random_sample() ← matches MATLAB randn
#   slow path inner:      1 int + 1×random_sample() ← matches MATLAB randn
# This ensures the MT state stays aligned with MATLAB after every call.

import math as _math

def _build_normal_ziggurat() -> tuple:
    """Compute kn/wn/fn tables for the normal ziggurat (float64, N=128)."""
    N   = 128
    M1  = 2_147_483_648.0   # 2^31
    dn  = 3.442619855899    # starting x-value (= r)
    tn  = dn
    vn  = 9.91256303526217e-3

    kn = np.empty(N, dtype=np.uint64)
    wn = np.empty(N, dtype=np.float64)
    fn = np.empty(N, dtype=np.float64)

    q      = vn / _math.exp(-0.5 * dn * dn)
    kn[0]  = int((dn / q) * M1)
    kn[1]  = 0
    wn[0]  = q / M1
    wn[127] = dn / M1
    fn[0]  = 1.0
    fn[127] = _math.exp(-0.5 * dn * dn)

    for i in range(126, 0, -1):          # i = 126 down to 1
        dn        = _math.sqrt(-2.0 * _math.log(vn / dn + _math.exp(-0.5 * dn * dn)))
        kn[i + 1] = int((dn / tn) * M1)
        tn        = dn
        fn[i]     = _math.exp(-0.5 * dn * dn)
        wn[i]     = dn / M1

    return kn, wn, fn


_ZIGG_KN, _ZIGG_WN, _ZIGG_FN = _build_normal_ziggurat()
_ZIGG_R     = 3.442619855899
_ZIGG_RINV  = 1.0 / _ZIGG_R          # 1/r used in tail sampling


def _mt_randn(rng: np.random.RandomState) -> float:
    """One N(0,1) sample via MATLAB-identical ziggurat.

    Draws exactly one raw 32-bit MT integer in the fast path (same as
    MATLAB's genrand_int32 call inside randn), keeping the MT state
    perfectly aligned with a concurrent MATLAB run.
    """
    while True:
        # One raw 32-bit MT integer — matches MATLAB's genrand_int32()
        u32 = rng.randint(0, 0x1_0000_0000, dtype=np.uint32)
        hz  = int(u32.view(np.int32))  # reinterpret bits as signed 32-bit
        iz  = int(u32) & 127           # strip index (bottom 7 bits)

        if abs(hz) < int(_ZIGG_KN[iz]):
            return hz * _ZIGG_WN[iz]   # fast path (~97.5 %)

        # ── Slow path ─────────────────────────────────────────────────────
        if iz == 0:
            # Tail sampling (x > r): exponential with rate 1/r
            while True:
                x = -_ZIGG_RINV * _math.log(rng.random_sample())  # 2 MT ints
                y = -_math.log(rng.random_sample())                 # 2 MT ints
                if y + y >= x * x:
                    return _ZIGG_R + x if hz > 0 else -(_ZIGG_R + x)
        else:
            x = hz * _ZIGG_WN[iz]
            if (_ZIGG_FN[iz]
                    + rng.random_sample() * (_ZIGG_FN[iz - 1] - _ZIGG_FN[iz])  # 2 MT ints
                    < _math.exp(-0.5 * x * x)):
                return x


def _mt_randn_vec(rng: np.random.RandomState, n: int) -> np.ndarray:
    """Return an array of n i.i.d. N(0,1) samples from the MATLAB ziggurat."""
    return np.array([_mt_randn(rng) for _ in range(n)], dtype=np.float64)


# ---------------------------------------------------------------------------
# Convergence diagnostics
# ---------------------------------------------------------------------------

def _r_hat(
    chains: np.ndarray,   # (T, n_params, N_chains)
    t_burnin: int,
    t: int,
) -> np.ndarray:
    """Compute Gelman-Rubin R-hat for each parameter.

    Uses the same formula as the MATLAB implementation:
      x_j_r_bar = 2/(T_eff-2) * sum(x_sel, axis=time)
    instead of the plain sample mean, with matching W/B/sigma coefficients.
    The slice includes the current sample at index t (matching MATLAB's 1-indexed
    x_rhat(idx_start:t, :, :) which includes row t).
    """
    idx_start = t_burnin + (t - t_burnin) // 2
    x_sel = chains[idx_start : t + 1]       # (T_eff, n_params, N_chains) — includes row t
    T_eff, n_params, N_chains = x_sel.shape

    if T_eff <= 2:
        return np.full(n_params, np.nan)

    # Chain "centre" — matches MATLAB: x_j_r_bar = 2/(T_eff-2) * sum(x_sel, 1)
    x_j_bar = (2.0 / (T_eff - 2)) * x_sel.sum(axis=0)   # (n_params, N_chains)
    x_bar_bar = x_j_bar.mean(axis=1, keepdims=True)       # (n_params, 1)

    # Within-chain variance
    W = ((x_sel - x_j_bar[np.newaxis]) ** 2).sum(axis=(0, 2)) / (N_chains * (T_eff - 2)) * 2

    # Between-chain variance
    B = T_eff / (2 * (N_chains - 1)) * ((x_j_bar - x_bar_bar) ** 2).sum(axis=1)

    sigma = ((T_eff - 2) / T_eff) * W + (2 / T_eff) * B
    with np.errstate(divide="ignore", invalid="ignore"):
        r = np.sqrt(((N_chains + 1) / N_chains) * sigma / W - (T_eff - 2) / (N_chains * T_eff))
    # W=0 → parameter constant across all chains (common for discrete nC early
    # in the chain).  Return NaN so plots show gaps; callers use the sentinel
    # version for convergence testing.
    return r


def _detect_outliers(
    p_x: np.ndarray,   # (t_current, N_chains)
    t: int,
    method: str,
) -> list[int]:
    """Return indices of outlier chains."""
    half = max(1, t // 2)
    omega = p_x[half:t].mean(axis=0)    # mean posterior per chain
    outliers = []
    if method == "mad":
        med = np.median(omega)
        mad = 1.4826 * np.median(np.abs(omega - med))
        for j in range(len(omega)):
            if omega[j] < med - 3 * mad:
                outliers.append(j)
    elif method == "iqr":
        q1, q3 = np.percentile(omega, [25, 75])
        iqr = q3 - q1
        for j in range(len(omega)):
            if omega[j] < q1 - 1.5 * iqr:
                outliers.append(j)
    return outliers


# ---------------------------------------------------------------------------
# Gibbs moves
# ---------------------------------------------------------------------------

def _build_gibbs_candidates(
    *,
    i: int,
    t: int,
    t_burnin: int,
    Nc: int,
    N_chains: int,
    components: np.ndarray,
    classes: list[list[Species]],
    n_ranges_list: list[list[int]],
    mol_temp: np.ndarray,
    nC_temp: np.ndarray,
    eta_temp: np.ndarray,
    nC_vals: np.ndarray,
    eta_vals: np.ndarray,
    max_comb: int | float,
    rng: np.random.RandomState,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Replicate the MATLAB Gibbs candidate construction more closely.

    Returns
    -------
    mol_all, nC_all, eta_all
        Candidate states to score with posterior_cheap_fn.
    """
    n_comp_update = len(components)

    pair_lists = []
    cur_idxs = []

    # MATLAB uses NBH = Inf and maxSamples = Inf*100 in the shared file,
    # so practically it keeps all local pairs.
    for cc in range(n_comp_update):
        j_comp = int(components[cc])
        sp_list = classes[j_comp]
        n_rng = n_ranges_list[j_comp]

        pair_list = []
        for nC_cand in n_rng:
            eta_norm_list = [sp.eta_B_star_norm for sp in sp_list if sp.nC == nC_cand]
            if len(eta_norm_list) == 0:
                continue
            for eta_cand in eta_norm_list:
                pair_list.append([nC_cand, eta_cand])

        pair_list = np.asarray(pair_list, dtype=float)

        cur_pair = np.array([nC_temp[j_comp], eta_temp[j_comp]], dtype=float)
        if pair_list.size == 0:
            pair_list = cur_pair[None, :]
            cur_idx = 0
        else:
            cur_rows = np.where(np.all(np.isclose(pair_list, cur_pair), axis=1))[0]
            if len(cur_rows) == 0:
                pair_list = np.vstack([pair_list, cur_pair])
                cur_idx = len(pair_list) - 1
            else:
                cur_idx = int(cur_rows[0])

        # MATLAB adds donor pairs to the per-component list after burn-in
        if t > t_burnin:
            donor_pairs = np.column_stack([
                nC_vals[np.arange(N_chains) != i, j_comp],
                eta_vals[np.arange(N_chains) != i, j_comp],
            ])
            if donor_pairs.size > 0:
                pair_list = np.vstack([pair_list, donor_pairs])
                pair_list = np.unique(pair_list, axis=0)

            # recompute current index after dedup/reorder
            cur_rows = np.where(np.all(np.isclose(pair_list, cur_pair), axis=1))[0]
            cur_idx = int(cur_rows[0])

        pair_lists.append(pair_list)
        cur_idxs.append(cur_idx)

    # ---------- Cartesian product with MATLAB-like stratified truncation ----------
    num_comb_each = [pl.shape[0] for pl in pair_lists]
    n_comb_full = int(np.prod(num_comb_each)) if len(num_comb_each) > 0 else 1

    if np.isfinite(max_comb) and n_comb_full > max_comb:
        idx_sub = []
        targ = max(1, int(np.ceil(max_comb ** (1.0 / n_comp_update))))

        for cc in range(n_comp_update):
            P = pair_lists[cc]
            n = P.shape[0]
            inc = cur_idxs[cc]

            if n <= targ:
                idx_sub.append(np.arange(n, dtype=int))
                continue

            # Group by nC as in MATLAB
            nC_vals_cc = P[:, 0]
            nC_unique, grp = np.unique(nC_vals_cc, return_inverse=True)
            G = len(nC_unique)
            targ_cc = min(n, max(targ, G))

            rows_by_grp = [np.where(grp == g)[0] for g in range(G)]

            selected = []
            # one random row per nC group
            for rows in rows_by_grp:
                selected.append(int(rows[rng.randint(len(rows))]))

            # ensure current included
            if inc not in selected:
                g_inc = grp[inc]
                repl = next((k for k, idx in enumerate(selected) if grp[idx] == g_inc), None)
                if repl is not None:
                    selected[repl] = inc
                else:
                    selected.append(int(inc))

            # fill remaining quota round-robin
            rem_by_grp = []
            for rows in rows_by_grp:
                rem = np.array([r for r in rows if r not in selected], dtype=int)
                if len(rem) > 0:
                    rem = rem[rng.permutation(len(rem))]
                rem_by_grp.append(rem)

            g_ptr = np.ones(G, dtype=int) * 0
            while len(selected) < targ_cc:
                progressed = False
                for g in range(G):
                    if len(selected) >= targ_cc:
                        break
                    if g_ptr[g] < len(rem_by_grp[g]):
                        selected.append(int(rem_by_grp[g][g_ptr[g]]))
                        g_ptr[g] += 1
                        progressed = True
                if not progressed:
                    break

            idx_sub.append(np.asarray(selected, dtype=int))
    else:
        idx_sub = [np.arange(pl.shape[0], dtype=int) for pl in pair_lists]

    grids = np.meshgrid(*idx_sub, indexing="ij")
    flat = [g.ravel() for g in grids]
    idx_pairs = np.column_stack(flat) if len(flat) > 1 else flat[0][:, None]
    n_comb = idx_pairs.shape[0]

    mol_all = np.tile(mol_temp, (n_comb, 1))
    nC_all = np.tile(nC_temp, (n_comb, 1)).astype(float)
    eta_all = np.tile(eta_temp, (n_comb, 1))

    for cc in range(n_comp_update):
        j_comp = int(components[cc])
        P = pair_lists[cc]
        nC_all[:, j_comp] = P[idx_pairs[:, cc], 0]
        eta_all[:, j_comp] = P[idx_pairs[:, cc], 1]

    # ---------- ensure current joint state exists ----------
    cur_idx_each = np.asarray(cur_idxs, dtype=int)
    if n_comp_update == 1:
        cur_idx_in_grid = int(cur_idx_each[0])
    else:
        matches = np.where(np.all(idx_pairs == cur_idx_each[None, :], axis=1))[0]
        if len(matches) == 0:
            cur_idx_in_grid = n_comb
            mol_all = np.vstack([mol_all, mol_temp[None, :]])
            nC_all = np.vstack([nC_all, nC_temp[None, :].astype(float)])
            eta_all = np.vstack([eta_all, eta_temp[None, :]])
        else:
            cur_idx_in_grid = int(matches[0])

    # ---------- after burn-in append other chains' full joint states ----------
    if t > t_burnin:
        other_idx = np.where(np.arange(N_chains) != i)[0]
        if len(other_idx) > 0:
            donor_mol = np.tile(mol_temp, (len(other_idx), 1))
            donor_nC = nC_vals[other_idx].astype(float)
            donor_eta = eta_vals[other_idx].copy()

            mol_all = np.vstack([mol_all, donor_mol])
            nC_all = np.vstack([nC_all, donor_nC])
            eta_all = np.vstack([eta_all, donor_eta])

            idx_all_states = np.arange(len(mol_all) - len(other_idx), len(mol_all), dtype=int)
        else:
            idx_all_states = np.array([], dtype=int)

        n_comb = len(mol_all)
        if np.isfinite(max_comb) and n_comb > max_comb:
            mandatory = np.unique(np.concatenate([[cur_idx_in_grid], idx_all_states]))
            n_free = max(0, int(max_comb) - len(mandatory))

            pool = np.setdiff1d(np.arange(n_comb, dtype=int), mandatory, assume_unique=False)
            if n_free > 0 and len(pool) > 0:
                take = rng.permutation(pool)[: min(n_free, len(pool))]
                idx_keep = np.unique(np.concatenate([mandatory, take]))
            else:
                idx_keep = mandatory

            mol_all = mol_all[idx_keep]
            nC_all = nC_all[idx_keep]
            eta_all = eta_all[idx_keep]
    else:
        n_comb = len(mol_all)
        if np.isfinite(max_comb) and n_comb > max_comb:
            mandatory = np.array([cur_idx_in_grid], dtype=int)
            n_free = max(0, int(max_comb) - len(mandatory))
            pool = np.setdiff1d(np.arange(n_comb, dtype=int), mandatory, assume_unique=False)
            if n_free > 0 and len(pool) > 0:
                take = rng.permutation(pool)[: min(n_free, len(pool))]
                idx_keep = np.unique(np.concatenate([mandatory, take]))
            else:
                idx_keep = mandatory

            mol_all = mol_all[idx_keep]
            nC_all = nC_all[idx_keep]
            eta_all = eta_all[idx_keep]

    return mol_all, nC_all, eta_all


# ---------------------------------------------------------------------------
# Main DE-MC function
# ---------------------------------------------------------------------------

def run_demc(
    families: list[str],
    posterior_fn: PostFn,
    posterior_cheap_fn: PostFn,
    classes: list[list[Species]],
    lower_bound_x: np.ndarray,
    upper_bound_x: np.ndarray,
    n_ranges: list[Sequence[int]],
    lower_bound_eta: np.ndarray,
    upper_bound_eta: np.ndarray,
    max_iterations: int,
    t_burnin: int,
    t_check: int,
    scaling_factor_x: float,
    scaling_factor_nc: float | np.ndarray,
    scaling_factor_eta: float,
    noise_x: float,
    noise_nc: float,
    noise_eta: float,
    N_chains: int,
    outlier_method: str = "mad",
    R_hat_threshold: float = 1.1,
    PT_switch: str = "True",
    beta_min: float = 0.0001,
    T_ladder: str = "geometric",
    swap_freq: int = 20,
    n_pairs: int = 1,
    p_gibbs: float = 0.0,
    p_snooker: float = 0.0,
    log_file: str = "MCMC.txt",
    seed: int = 1234,
    x_init: np.ndarray | None = None,
    layout: ParameterLayout | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    """
    Run the DE-MC algorithm with delayed acceptance and parallel tempering.

    Parameters mirror the MATLAB function signature exactly, plus ``layout``:
    a ParameterLayout describing the state vector. If None, the standard
    single-surrogate SurrogateLayout (3*Nc-1 parameters) is built from
    classes, n_ranges and the molar-fraction bounds.

    Returns
    -------
    chain : ndarray, shape (t_conv, n_params, N_chains)
        Physical parameter chains (molar fractions, nC, eta) in the column
        order of the layout; n_params = 3*Nc-1 for the standard layout.
    posterior_pdf : ndarray, shape (t_conv, N_chains)
        Log-posterior values.
    chain_reshaped : ndarray, shape (N_post * N_chains, n_params)
        Post-burnin samples from all chains concatenated.
    posterior_pdf_reshaped : ndarray, shape (N_post * N_chains,)
        Log-posterior for reshaped chain.
    AR : ndarray, shape (t_conv, N_chains)
        Acceptance rates (%).
    R_hat : ndarray, shape (n_checks, n_params)
        R-hat statistics.
    t_convergence : int
        Iteration at which convergence was declared.
    """
    rng = np.random.RandomState(seed)  # MT19937 — matches MATLAB rng(seed,'twister')
    if layout is None:
        layout = SurrogateLayout(classes, n_ranges, lower_bound_x, upper_bound_x)
    n_params = layout.n_params
    n_units = len(layout.units)

    # Start a persistent worker pool.  Passing *classes* to the initializer
    # pre-loads the species database in each worker so it is not re-pickled
    # on every pool.map() call.
    init_pool(classes=classes)

    fid = open(log_file, "w") if log_file is not None else None
    _hdr = (
        "    ____                  _____ ___    ______\n"
        "   / __ )____ ___  _____ / ___//   |  / ____/\n"
        "  / __  / __ `/ / / / _ \\__ \\/ /| | / /_\n"
        " / /_/ / /_/ / /_/ /  __/__/ / ___ |/ __/\n"
        "/_____/\\__,_/\\__, /\\___/____/_/  |_/_/\n"
        "            /____/\n\n"
        "BayeSAF: Emulation and Design of Sustainable Alternative Fuels\n"
        "via Bayesian Inference and Descriptors-Based Machine Learning\n\n"
        "Contributors / Copyright Notice\n"
        "© 2026 Jacopo Liberatori — jacopo.liberatori@centralesupelec.fr\n"
        "Postdoctoral Researcher @ Laboratoire EM2C, CentraleSupélec (CNRS)\n\n"
        "© 2026 Davide Cavalieri — davide.cavalieri@uniroma1.it\n"
        "Postdoctoral Researcher @ Sapienza University of Rome,\n"
        "Department of Mechanical and Aerospace Engineering (DIMA)\n\n"
        "© 2026 Matteo Blandino, Ph.D.\n\n"
        "Reference:\n"
        "J. Liberatori, D. Cavalieri, M. Blandino, M. Valorani, and P.P. Ciottoli.\n"
        "BayeSAF: Emulation and Design of Sustainable Alternative Fuels via Bayesian\n"
        "Inference and Descriptors-Based Machine Learning. Fuel 419, 138835 (2026).\n"
        "Available at: https://doi.org/10.1016/j.fuel.2026.138835.\n\n"
        "-----------------------------------------------------------------------\n\n"
        f"+++++ Differential evolution Markov chain (DE-MC) algorithm +++++\n"
        f"Maximum number of iterations: {max_iterations}\n"
        f"Number of chains: {N_chains}\n"
    )
    if fid is not None:
        fid.write(_hdr)

    # ── Preallocate ─────────────────────────────────────────────────────────
    # phys_chain stores physical values (molar fractions, nC, eta_B_star_norm)
    # with the column layout of `layout`; x_archive stores normalised states
    # for the snooker jump.
    phys_chain = np.full((max_iterations, n_params, N_chains), np.nan)
    p_x = np.full((max_iterations, N_chains), np.nan)
    accept_arr = np.full((max_iterations, N_chains), np.nan)
    AR = np.full((max_iterations, N_chains), np.nan)
    R_hat = np.full((max_iterations // 2, n_params), 1e18)
    x_archive = np.zeros((max_iterations, n_params, N_chains))

    # ── Initialise population ───────────────────────────────────────────────
    if x_init is not None:
        # External initial states (benchmark mode) — bypass random sampling.
        X, phys, p_X = layout.decode_init(x_init, posterior_fn)
    else:
        X, phys, p_X = layout.init_population(N_chains, posterior_fn, rng)
    # p_X is a 1-D array (N_chains,)

    t = 0  # 0-based index
    phys_chain[t] = phys.T
    p_x[t] = p_X
    accept_arr[t] = 1.0
    AR[t] = 100.0

    # R-matrix: each chain's "other chains" indices
    R_mat = np.array([np.delete(np.arange(N_chains), i) for i in range(N_chains)])

    # Temperature ladder
    if PT_switch == "True":
        if T_ladder == "linear":
            beta_arr = np.linspace(1.0, beta_min, N_chains)
        else:
            beta_arr = beta_min ** (np.arange(N_chains) / (N_chains - 1))
            beta_arr = np.sort(beta_arr)[::-1]
    else:
        beta_arr = np.ones(N_chains)
        swap_freq = max_iterations + 1  # effectively disable swaps

    convergence = False
    t_convergence = max_iterations
    counter_check = 0

    # ── R-hat skip mask ───────────────────────────────────────────────────────
    # Parameters whose value is structurally constant across all chains
    # (e.g. eta for n-paraffins, which have exactly one isomer per nC) will
    # always have W=0, producing NaN R-hat.  We identify them once here so we
    # can:
    #   • store NaN in R_hat (shows as a gap in the diagnostic plot)
    #   • treat them as trivially converged in the convergence test
    #
    # A parameter is constant if all species in that family share a single
    # distinct eta_B_star_norm value (regardless of nC).
    _rhat_skip = layout.rhat_skip()

    if fid is not None and _rhat_skip.any():
        skipped = [families[k] for unit in layout.units for k in np.where(unit.constant_eta())[0]]
        fid.write(
            f"R-hat skipped (constant eta) for families: {skipped}\n"
        )

    # ── Main loop ────────────────────────────────────────────────────────────
    for _iter in range(1, max_iterations):
        t = _iter

        if t > t_burnin:
            beta_arr = np.ones(N_chains)
            swap_freq_eff = max_iterations + 1
        else:
            swap_freq_eff = swap_freq

        lambda_vec = rng.uniform(-0.1, 0.1, N_chains)
        draw = np.argsort(rng.random_sample((N_chains - 1, N_chains)), axis=0)

        Xp = X.copy()
        phys_new = phys.copy()

        J = np.zeros(N_chains)
        log_pi1_x_vec = np.zeros(N_chains)
        log_pi1_y_vec = np.zeros(N_chains)

        # ── Per-chain proposal ─────────────────────────────────────────────
        for i in range(N_chains):
            r_gibbs = rng.random_sample()

            Xp_temp = Xp[i].copy()
            phys_temp = phys[i].copy()

            # Choose block
            if r_gibbs > p_gibbs:
                block = layout.blocks[rng.randint(0, len(layout.blocks))]
            else:
                block = layout.fractions_block(0 if n_units == 1 else rng.randint(0, n_units))
            unit = layout.units[block.unit]
            Nc = unit.Nc

            # Scaling factors and jump partner indices
            # MATLAB: D = randsample(1:n_pairs, 1, 'true');
            D = int(rng.randint(1, n_pairs + 1))

            # MATLAB:
            # a = R(i,draw(1:D,i));
            # b = R(i,draw(D+1:2*D,i));
            a_idx = R_mat[i][draw[:D, i]]
            b_idx = R_mat[i][draw[D:2 * D, i]]

            gamma_x = scaling_factor_x * 2.38 / np.sqrt(2 * D * (Nc - 1))
            gamma_nc = scaling_factor_nc * 2.38 / np.sqrt(2 * D * Nc)
            gamma_eta = scaling_factor_eta * 2.38 / np.sqrt(2 * D * Nc)

            g_x = gamma_x if rng.random_sample() < 0.9 else 1.0

            # MATLAB treats gamma_nc elementwise if it is a vector
            if np.isscalar(gamma_nc):
                g_nc = gamma_nc
                if rng.random_sample() >= 0.9:
                    g_nc = 1.0
            else:
                g_nc = np.asarray(gamma_nc, dtype=float).copy()
                mask = rng.random_sample(g_nc.shape) < 0.9
                g_nc[~mask] = 1.0

            g_eta = gamma_eta if rng.random_sample() < 0.9 else 1.0

            # ── Snooker-jump pre-computation ───────────────────────────────
            # Active only for the fractions block, after burn-in and after a
            # warm-up of 10*n_params steps (mirrors the MATLAB condition).
            r_move = rng.random_sample()
            use_snooker = (
                p_snooker > 0.0
                and r_move < p_snooker
                and t > 10 * n_params
                and t > t_burnin
                and block.kind == "fractions"
            )
            xR_snooker = None
            z_snooker = None
            z_proj = None
            g_snooker = 1.0

            if use_snooker:
                # Filter archive: discard all-zero time slices
                arc_norm = np.abs(x_archive).sum(axis=(1, 2))
                x_arc_filt = x_archive[arc_norm > 0]   # (T_filt, n_params, N_chains)
                T_filt = x_arc_filt.shape[0]

                if T_filt >= 3 and N_chains >= 4:
                    # Three distinct chains different from i
                    other_chains = np.delete(np.arange(N_chains), i)
                    r1, r2, r3 = other_chains[rng.permutation(len(other_chains))[:3]]

                    # Three independent archive rows
                    idx_samp = rng.randint(0, T_filt, 3)
                    xR_snooker = x_arc_filt[idx_samp[0], :, r1]
                    z_snooker  = xR_snooker - X[i]
                    v          = (x_arc_filt[idx_samp[1], :, r2]
                                  - x_arc_filt[idx_samp[2], :, r3])

                    alpha  = np.dot(z_snooker, v) / (np.dot(z_snooker, z_snooker) + 1e-8)
                    z_proj = alpha * z_snooker
                    g_snooker = rng.uniform(1.2, 2.2)
                else:
                    use_snooker = False   # not enough archive data yet

            cols = block.cols
            if block.kind == "fractions":
                if use_snooker:
                    # Snooker proposal (molar-fraction subspace only)
                    Xp_temp[cols] = (
                        X[i, cols]
                        + g_snooker * z_proj[cols]
                        + noise_x * _mt_randn_vec(rng, len(cols))
                    )
                    # Jacobian correction (ter Braak & Vrugt 2008, eq. 6)
                    d    = len(cols)
                    zp   = xR_snooker[cols] - Xp_temp[cols]
                    epsJ = 1e-12
                    J[i] = (d - 1) * (
                        np.log(np.linalg.norm(z_snooker[cols]) + epsJ)
                        - np.log(np.linalg.norm(zp) + epsJ)
                    )
                else:
                    # Standard DE move
                    for j in cols:
                        de_diff = np.sum(X[a_idx, j] - X[b_idx, j], axis=0)
                        Xp_temp[j] = (
                            X[i, j]
                            + (1 - lambda_vec[i]) * g_x * de_diff
                            + noise_x * _mt_randn(rng)
                        )

            elif block.kind == "nC":
                for k, j in enumerate(cols):
                    de_diff = np.sum(X[a_idx, j] - X[b_idx, j], axis=0)
                    if np.isscalar(g_nc):
                        jump_scale = g_nc
                    else:
                        jump_scale = g_nc[k]
                    Xp_temp[j] = (
                        X[i, j]
                        + (1 - lambda_vec[i]) * jump_scale * de_diff
                        + noise_nc * _mt_randn(rng)
                    )

            elif block.kind == "eta":
                for j in cols:
                    de_diff = np.sum(X[a_idx, j] - X[b_idx, j], axis=0)
                    Xp_temp[j] = (
                        X[i, j]
                        + (1 - lambda_vec[i]) * g_eta * de_diff
                        + noise_eta * _mt_randn(rng)
                        )

            # ── Boundary handling ──────────────────────────────────────────
            unit.fold_fractions(Xp_temp, phys_temp)

            # ── Gibbs-like update for discrete parameters ──────────────────
            last_comp_sum_ok = unit.fractions_sum_ok(phys_temp) and r_gibbs <= p_gibbs

            if last_comp_sum_ok:
                # Gibbs step: candidate construction
                n_comp_update = min(Nc, Nc if t > t_burnin else 1)
                max_comb = 10 * n_comp_update + (N_chains - 1) if t > t_burnin else np.inf
                components = rng.permutation(Nc)[:n_comp_update]

                mol_all, nC_all, eta_all = _build_gibbs_candidates(
                    i=i,
                    t=t,
                    t_burnin=t_burnin,
                    Nc=Nc,
                    N_chains=N_chains,
                    components=np.asarray(components, dtype=int),
                    classes=unit.classes,
                    n_ranges_list=unit.n_ranges,
                    mol_temp=phys_temp[unit.frac_cols],
                    nC_temp=phys_temp[unit.nc_cols],
                    eta_temp=phys_temp[unit.eta_cols],
                    nC_vals=phys[:, unit.nc_cols],
                    eta_vals=phys[:, unit.eta_cols],
                    max_comb=max_comb,
                    rng=rng,
                )

                # Embed the unit-level candidates into full state rows
                phys_cand = np.tile(phys_temp, (mol_all.shape[0], 1))
                phys_cand[:, unit.frac_cols] = mol_all
                phys_cand[:, unit.nc_cols] = nC_all
                phys_cand[:, unit.eta_cols] = eta_all

                log_p_all = posterior_cheap_fn(*layout.posterior_args(phys_cand))
                max_lp = np.max(log_p_all)
                if not np.isfinite(max_lp):
                    # All candidates have -inf posterior (e.g. all flash points ≤ 0 °F);
                    # fall back to uniform weights so the chain can still move.
                    w = np.ones(len(log_p_all)) / len(log_p_all)
                else:
                    w = np.exp(log_p_all - max_lp)
                    w /= np.sum(w)

                idx_sel = int(rng.choice(len(w), p=w))
                log_pi1_y_vec[i] = log_p_all[idx_sel]

                for cc in range(n_comp_update):
                    j_comp = int(components[cc])
                    nC_sel = int(nC_all[idx_sel, j_comp])
                    eta_sel = eta_all[idx_sel, j_comp]
                    phys_temp[unit.nc_cols[j_comp]] = nC_sel
                    phys_temp[unit.eta_cols[j_comp]] = eta_sel
                    # Update normalized X for nC and eta
                    unit.encode_discrete(Xp_temp, j_comp, nC_sel, eta_sel)

                log_pi1_x_vec[i] = posterior_cheap_fn(
                    *layout.posterior_args(phys[i:i+1])
                )[0]

            elif block.kind in ("nC", "eta"):
                # Continuous-to-discrete mapping with triangle-fold
                unit.decode_discrete(Xp_temp, phys_temp)

            Xp[i] = Xp_temp
            phys_new[i] = phys_temp

        # ── Delayed-acceptance ─────────────────────────────────────────────
        log_alpha1 = beta_arr * (log_pi1_y_vec - log_pi1_x_vec) + J
        u1 = np.log(rng.random_sample(N_chains))
        pass1 = u1 < log_alpha1

        accept_arr[t] = accept_arr[t - 1].copy()

        if pass1.any():
            idx_pass = np.where(pass1)[0]

            p_y_full = posterior_fn(*layout.posterior_args(phys_new[idx_pass]))

            delta_full = p_y_full - p_X[idx_pass]
            delta_cheap = log_pi1_y_vec[idx_pass] - log_pi1_x_vec[idx_pass]
            log_alpha2 = beta_arr[idx_pass] * (delta_full - delta_cheap)

            u2 = np.log(rng.random_sample(len(idx_pass)))
            accept2 = u2 < log_alpha2

            if accept2.any():
                idx_acc = idx_pass[accept2]
                X[idx_acc] = Xp[idx_acc]
                phys[idx_acc] = phys_new[idx_acc]
                p_X[idx_acc] = p_y_full[accept2]
                accept_arr[t, idx_acc] = accept_arr[t - 1, idx_acc] + 1

        AR[t] = 100.0 * accept_arr[t] / t

        # ── Parallel tempering swaps ───────────────────────────────────────
        if swap_freq_eff <= max_iterations and (t % swap_freq_eff) == 0:
            swap_idx_int = t // swap_freq_eff
            pair_start = 0 if (swap_idx_int % 2 == 0) else 1
            for k in range(pair_start, N_chains - 1, 2):
                Delta = (beta_arr[k] - beta_arr[k + 1]) * (p_X[k + 1] - p_X[k])
                if Delta >= 0 or np.log(rng.random_sample()) < Delta:
                    X[[k, k + 1]] = X[[k + 1, k]]
                    p_X[[k, k + 1]] = p_X[[k + 1, k]]
                    # Keep physical representations consistent with the swapped
                    # normalised state so that chain_reshaped (which is built from
                    # phys_chain) always matches p_x.
                    phys[[k, k + 1]] = phys[[k + 1, k]]

        # ── Outlier detection during burn-in ──────────────────────────────
        if t in (int(0.25 * t_burnin), int(0.5 * t_burnin),
                 int(0.75 * t_burnin), t_burnin):
            outliers = _detect_outliers(p_x[:t], t, outlier_method)
            best = int(np.argmax(p_X))
            for j_out in outliers:
                X[j_out] = X[best].copy()
                p_X[j_out] = p_X[best]
                phys[j_out] = phys[best].copy()
            if outliers and fid is not None:
                fid.write(f"Chains {outliers} detected as outliers at t={t}. "
                          f"Reset to chain {best}.\n")

        # ── Store current state ────────────────────────────────────────────
        phys_chain[t] = phys.T
        p_x[t] = p_X

        if t % 10 == 0:
            x_archive[t] = X.T

        # ── R-hat convergence check ────────────────────────────────────────
        if t > t_check and t % 2 == 0:
            # Use physical parameter chains for R-hat
            rhat = _r_hat(phys_chain, t_burnin, t)

            # Parameters with structurally constant eta (e.g. n-paraffins):
            # force NaN in storage (gap in plot) and 1.0 in the convergence
            # test (trivially converged — no mixing needed for a fixed param).
            rhat[_rhat_skip] = np.nan
            R_hat[counter_check] = rhat

            # For convergence test: NaN/inf from other discrete params → large
            # sentinel; constant-eta params (already NaN) → 1.0.
            rhat_test = np.where(_rhat_skip, 1.0, rhat)
            rhat_test = np.where(np.isfinite(rhat_test), rhat_test, 1e18)
            if np.all(rhat_test <= R_hat_threshold):
                t_convergence = t
                convergence = True
                if fid is not None:
                    fid.write(
                        f"R-hat below {R_hat_threshold} for all parameters: "
                        f"convergence after {t_convergence} iterations.\n"
                    )
                break

            counter_check += 1

        if t % 100 == 0:
            print(f"\r  DE-MC iteration {t}/{max_iterations}  AR={AR[t].mean():.1f}%", end="", flush=True)

    print()  # newline after progress

    if fid is not None:
        if not convergence:
            fid.write(f"Convergence not reached. Using all {max_iterations} samples.\n")
        fid.close()

    shutdown_pool()

    # ── Trim arrays to t_convergence ─────────────────────────────────────
    tc = t_convergence
    x_chain = phys_chain[:tc]

    p_x = p_x[:tc]
    AR = AR[:tc]
    R_hat = R_hat[:counter_check + 1]   # +1: include the row that triggered convergence

    # ── Post-burnin samples ───────────────────────────────────────────────
    N_post = tc - t_burnin
    chain_reshaped = np.vstack(
        [x_chain[t_burnin:tc, :, ch] for ch in range(N_chains)]
    )
    posterior_reshaped = np.concatenate(
        [p_x[t_burnin:tc, ch] for ch in range(N_chains)]
    )

    return (
        x_chain,
        p_x,
        chain_reshaped,
        posterior_reshaped,
        AR,
        R_hat,
        t_convergence,
    )
