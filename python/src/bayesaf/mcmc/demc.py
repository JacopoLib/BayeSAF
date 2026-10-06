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
acceptance, snooker jumps, blocked updates, an exact multiple-try Gibbs move
on carbon atoms and isomers with a uniform prior over isomers, and an adaptive
burn-in with outlier-chain resets. Convergence is monitored via the classic
Gelman-Rubin R-hat statistic and the effective sample size.

References:
ter Braak, C.J.F. (2006). A Markov chain Monte Carlo version of the genetic
algorithm differential evolution: easy Bayesian computing for real parameter
spaces. Stat. Comput. 16, 239-249.

Vrugt, J.A. et al. (2008, 2009, 2016). DREAM toolbox.

Liu, J.S., Liang, F., Wong, W.H. (2000). The multiple-try method and local
optimization in Metropolis sampling. J. Am. Stat. Assoc. 95, 121-134.

Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., Bürkner, P.-C. (2021).
Rank-normalization, folding, and localization: an improved R-hat for
assessing convergence of MCMC. Bayesian Anal. 16, 667-718.

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
11) t_burnin            : number of iterations of the parallel-tempering phase
                          of the burn-in
12) scaling_factor_X    : jump-rate scaling for molar fractions
13) scaling_factor_nc   : jump-rate scaling for numbers of carbon atoms
14) scaling_factor_eta  : jump-rate scaling for topochemical atom indices
15) noise_X             : noise parameter for molar-fraction proposals
16) noise_nc            : noise parameter for carbon-atom proposals
17) noise_eta           : noise parameter for topochemical-index proposals
18) N_chains            : number of chains
19) beta_min            : minimum inverse temperature during the burn-in
20) T_ladder            : temperature ladder form ('linear' or 'geometric')
21) swap_freq           : parallel tempering swap frequency during the burn-in
22) n_pairs             : maximum number of chain pairs used to propose the new sample
23) p_gibbs             : probability of a Gibbs move
24) p_snooker           : probability of a snooker jump
25) gibbs_tries         : number of tries of the multiple-try Gibbs move
                          (None: enumeration of all candidates)
26) n_cold              : number of chains at beta = 1 after the burn-in
27) beta_min_sampling   : minimum inverse temperature after the burn-in
28) swap_every_sampling : parallel tempering swap frequency after the burn-in
29) outlier_every       : iterations between two outlier-chain checks
30) burnin_cap          : maximum burn-in length, as a fraction of maxIterations
31) rhat_stop           : R-hat convergence threshold
32) rhat_checks         : number of consecutive R-hat checks below rhat_stop
33) rhat_spacing        : iterations between two R-hat checks
34) rhat_min_samples    : minimum number of post-burn-in iterations
35) ess_min             : minimum effective sample size of the molar fractions

Outputs:
1) x            : (t_convergence x 3*Nc-1 x N_chains) array of posterior
                  samples (molar fractions, nC, eta_B_star_norm)
2) p_x          : (t_convergence x N_chains) array of log-posterior values
3) x_nb         : post-burnin samples, shape ((t_convergence-t_burnin)*N_chains
                  x 3*Nc-1)
4) p_x_nb       : post-burnin log-posterior values
5) AR            : (t_convergence x N_chains) acceptance rate array
6) R_hat         : classic R-hat statistic array, shape
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

def _r_hat_classic(chains: np.ndarray, t_burnin: int, t: int) -> np.ndarray:
    """Classic Gelman-Rubin R-hat (plain chain means) over the second half of [t_burnin, t]."""
    idx_start = t_burnin + (t - t_burnin) // 2
    x = chains[idx_start: t + 1]              # (n, n_params, N_chains)
    n = x.shape[0]
    if n <= 2:
        return np.full(x.shape[1], np.nan)
    W = x.var(axis=0, ddof=1).mean(axis=1)
    B = n * x.mean(axis=0).var(axis=1, ddof=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.sqrt(((n - 1) / n * W + B / n) / W)


def _ess_bulk(draws: np.ndarray) -> float:
    """Bulk effective sample size of draws (n_iterations, n_chains) (Vehtari et al., 2021)."""
    from scipy.stats import norm, rankdata
    d = np.asarray(draws, dtype=float)
    if np.allclose(d, d.flat[0]):
        return float("inf")
    d = norm.ppf((rankdata(d, method="average").reshape(d.shape) - 0.375) / (d.size + 0.25))
    n = d.shape[0] // 2
    d = np.hstack([d[:n], d[n:2 * n]])
    m = d.shape[1]
    x = d - d.mean(0)
    f = np.fft.rfft(x, n=2 * n, axis=0)
    acov = np.fft.irfft(f * np.conj(f), axis=0)[:n] / n          # (n, m)
    W = (acov[0] * n / (n - 1)).mean(); B = n * d.mean(0).var(ddof=1)
    var_plus = (n - 1) / n * W + B / n
    rho = 1.0 - (W - acov.mean(1)) / var_plus
    rho[0] = 1.0
    P = rho[:-1:2] + rho[1::2]
    k = int(np.argmax(P < 0)) if np.any(P < 0) else len(P)
    P = np.minimum.accumulate(P[:k]) if k > 0 else np.array([1.0])
    tau = -1.0 + 2.0 * P.sum()
    return float(m * n / max(tau, 1.0 / np.log10(m * n)))


# ---------------------------------------------------------------------------
# Prior over isomers
# ---------------------------------------------------------------------------

def _uniform_isomer_prior(layout, posterior_fn, posterior_cheap_fn):
    """Posteriors with a uniform prior over the isomers of each component given
    its number of carbon atoms, in place of the decoder-cell measure."""
    unit = layout.units[0]
    table = []                                   # per component: {(nC, eta): log correction}
    for k in range(unit.Nc):
        corr = {}
        for nC in unit.n_ranges[k]:
            cells = {e: m for e, (m, _) in unit.eta_cells(k, int(nC)).items() if m > 0}
            for e, m in cells.items():
                corr[(float(nC), e)] = -np.log(len(cells)) - np.log(m)
        table.append(corr)

    def correction(nC, eta):
        out = np.zeros(nC.shape[0])
        for r in range(nC.shape[0]):
            for k in range(unit.Nc):
                out[r] += table[k].get((float(nC[r, k]), float(eta[r, k])), -np.inf)
        return out

    def wrap(f):
        def g(X, nC, eta):
            return f(X, nC, eta) + correction(nC, eta)
        return g
    return wrap(posterior_fn), wrap(posterior_cheap_fn)


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
    scaling_factor_x: float,
    scaling_factor_nc: float | np.ndarray,
    scaling_factor_eta: float,
    noise_x: float,
    noise_nc: float,
    noise_eta: float,
    N_chains: int,
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
    info: dict | None = None,
    gibbs_tries: int | None = None,
    n_cold: int | None = None,
    beta_min_sampling: float = 0.05,
    swap_every_sampling: int = 1,
    outlier_every: int | None = None,
    burnin_cap: float = 0.75,
    rhat_stop: float = 1.1,
    rhat_checks: int = 3,
    rhat_spacing: int = 500,
    rhat_min_samples: int = 2000,
    ess_min: float = 400,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    """
    Run the DE-MC algorithm with delayed acceptance and parallel tempering.

    Parameters mirror the MATLAB function signature, plus the optional
    arguments below.

    layout : ParameterLayout describing the state vector. If None, the
        standard SurrogateLayout (3*Nc-1 parameters) is built from classes,
        n_ranges and the molar-fraction bounds.
    info : dict receiving t_burnin, n_resets, swap_rates, betas and
        ESS_min_at_stop.
    gibbs_tries : number of uniform tries of the multiple-try Gibbs move on
        the (nC, eta) of one component; None enumerates all the candidates.
        The molar fractions first take a symmetric DE jump.
    n_cold, beta_min_sampling, swap_every_sampling : parallel tempering after
        the burn-in. The first n_cold chains (default N_chains // 2) stay at
        beta = 1 and are the only ones returned.
    outlier_every, burnin_cap : burn-in. Parallel tempering up to t_burnin;
        then, every outlier_every iterations, the chains whose mean
        log-posterior lies more than 5*sqrt(n_params/2) below the best chain
        are reset to good chains drawn at random. The burn-in ends at the
        first check without resets after 2*t_burnin, or at
        burnin_cap*max_iterations.
    rhat_stop, rhat_checks, rhat_spacing, rhat_min_samples, ess_min : the run
        stops once the classic Gelman-Rubin R-hat is below rhat_stop at
        rhat_checks checks rhat_spacing iterations apart, after at least
        rhat_min_samples post-burn-in iterations, and the bulk effective
        sample size (Vehtari et al., 2021) of every molar fraction is at
        least ess_min.

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
        Classic R-hat statistics, every 2 iterations after the burn-in.
    t_convergence : int
        Iteration at which convergence was declared.
    """
    rng = np.random.RandomState(seed)  # MT19937 — matches MATLAB rng(seed,'twister')
    if layout is None:
        layout = SurrogateLayout(classes, n_ranges, lower_bound_x, upper_bound_x)
    posterior_fn, posterior_cheap_fn = _uniform_isomer_prior(layout, posterior_fn, posterior_cheap_fn)
    n_params = layout.n_params

    # Start a persistent worker pool.  Passing *classes* to the initializer
    # pre-loads the species database in each worker so it is not re-pickled
    # on every pool.map() call.
    init_pool(classes=classes, batch_size=N_chains, nc=len(classes))

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
    rc_hist: dict = {}            # t -> max classic R-hat
    _ess_cols = sorted(int(c) for u in layout.units for c in u.frac_cols)
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

    # Parallel tempering after the burn-in: cold group + hot ladder
    n_cold = n_cold or N_chains // 2
    n_hot = N_chains - n_cold
    if n_cold < 4 or n_hot < 4:
        raise ValueError("at least 4 cold and 4 hot chains are needed")
    beta_prod = np.concatenate([np.ones(n_cold),
                                beta_min_sampling ** (np.arange(1, n_hot + 1) / n_hot)])
    cold_idx = np.arange(n_cold)
    group_of = np.where(np.arange(N_chains) < n_cold, 0, 1)
    groups = [cold_idx, np.arange(n_cold, N_chains)]
    swap_try = np.zeros(n_hot)
    swap_acc = np.zeros(n_hot)

    # Temperature ladder during the burn-in
    if T_ladder == "linear":
        beta_arr = np.linspace(1.0, beta_min, N_chains)
    else:
        beta_arr = beta_min ** (np.arange(N_chains) / (N_chains - 1))
        beta_arr = np.sort(beta_arr)[::-1]

    convergence = False
    t_convergence = max_iterations
    counter_check = 0

    # Effective burn-in
    t_burn = t_burnin
    burnin_done = False
    n_resets = 0
    every = outlier_every or max(50, t_burnin // 4)
    margin_best = 5.0 * np.sqrt(n_params / 2.0)

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

        in_prod = burnin_done
        if in_prod:
            beta_arr = beta_prod
            swap_freq_eff = max_iterations + 1
        elif t > t_burnin:
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
        stage1_override = np.full(N_chains, np.nan)   # Gibbs move: log alpha1

        # ── Per-chain proposal ─────────────────────────────────────────────
        for i in range(N_chains):
            r_gibbs = rng.random_sample()

            Xp_temp = Xp[i].copy()
            phys_temp = phys[i].copy()

            # ── Gibbs move on one component's (nC, eta) ───────────────────
            if r_gibbs <= p_gibbs:
                unit = layout.units[0]
                k = rng.randint(0, unit.Nc)
                cn, ce = unit.nc_cols[k], unit.eta_cols[k]
                cx = posterior_cheap_fn(*layout.posterior_args(phys[i:i + 1]))[0]
                moved = False
                if len(unit.frac_cols) > 0:
                    # Symmetric DE jump of the unit's molar fractions, x -> x'
                    if in_prod:
                        partners = groups[group_of[i]][groups[group_of[i]] != i]
                        pick = partners[rng.permutation(len(partners))[:2]]
                        a1, b1 = pick[:1], pick[1:2]
                    else:
                        a1, b1 = R_mat[i][draw[:1, i]], R_mat[i][draw[1:2, i]]
                    g1 = scaling_factor_x * 2.38 / np.sqrt(2 * len(unit.frac_cols))
                    g1 = g1 if rng.random_sample() < 0.9 else 1.0
                    for j in unit.frac_cols:
                        Xp_temp[j] = (X[i, j] + (1 - lambda_vec[i]) * g1 * np.sum(X[a1, j] - X[b1, j])
                                      + noise_x * _mt_randn(rng))
                    unit.fold_fractions(Xp_temp, phys_temp)
                    moved = True
                    if not unit.fractions_sum_ok(phys_temp):
                        stage1_override[i] = -np.inf
                        Xp[i], phys_new[i] = Xp_temp, phys_temp
                        continue
                if gibbs_tries is None:
                    cand = unit.candidates(k)
                    rows = np.tile(phys_temp, (len(cand), 1))
                    rows[:, cn], rows[:, ce] = cand[:, 0], cand[:, 1]
                    cu = posterior_cheap_fn(*layout.posterior_args(rows))
                    lcell = layout.log_cell(rows)
                    lw = cu + lcell
                    if moved:
                        rows_r = np.tile(phys[i], (len(cand), 1))
                        rows_r[:, cn], rows_r[:, ce] = cand[:, 0], cand[:, 1]
                        lw_r = posterior_cheap_fn(*layout.posterior_args(rows_r)) + lcell
                    finite = np.isfinite(lw)
                    if finite.any():
                        m_f = lw[finite].max()
                        w = np.where(finite, np.exp(lw - m_f), 0.0)
                        sel = int(rng.choice(len(w), p=w / w.sum()))
                        phys_temp[cn], phys_temp[ce] = cand[sel]
                        Xp_temp[cn], Xp_temp[ce] = unit.sample_u(k, cand[sel, 0], cand[sel, 1], rng)
                        log_pi1_y_vec[i], log_pi1_x_vec[i] = cu[sel], cx
                        a1_log = (beta_arr[i] - 1.0) * (cu[sel] - cx)
                        if moved:
                            fr = np.isfinite(lw_r)
                            log_z_f = m_f + np.log(w.sum())
                            log_z_r = (lw_r[fr].max() + np.log(np.sum(np.exp(lw_r[fr] - lw_r[fr].max())))
                                       if fr.any() else np.inf)
                            a1_log += log_z_f - log_z_r
                        stage1_override[i] = a1_log
                    else:
                        stage1_override[i] = -np.inf
                else:
                    kt = int(gibbs_tries)
                    def _tries(us, base):
                        rows = np.tile(base, (len(us), 1))
                        for r, (a, b) in enumerate(us):
                            rows[r, cn], rows[r, ce] = unit.decode_component(k, a, b)
                        return rows, posterior_cheap_fn(*layout.posterior_args(rows))
                    u_c = rng.random_sample((kt, 2))
                    rows_c, cu_c = _tries(u_c, phys_temp)   # tries at x'
                    lw_c = beta_arr[i] * cu_c
                    if np.isfinite(lw_c).any():
                        m = lw_c[np.isfinite(lw_c)].max()
                        w = np.where(np.isfinite(lw_c), np.exp(lw_c - m), 0.0)
                        sel = int(rng.choice(kt, p=w / w.sum()))
                        u_r = np.vstack([rng.random_sample((kt - 1, 2)), [[X[i, cn], X[i, ce]]]])
                        _, cu_r = _tries(u_r[:-1], phys[i])   # references at x
                        lw_r = beta_arr[i] * np.append(cu_r, cx)
                        m_r = lw_r[np.isfinite(lw_r)].max()
                        log_sum_c = m + np.log(w.sum())
                        log_sum_r = m_r + np.log(np.sum(np.exp(lw_r[np.isfinite(lw_r)] - m_r)))
                        phys_temp[cn], phys_temp[ce] = rows_c[sel, cn], rows_c[sel, ce]
                        Xp_temp[cn], Xp_temp[ce] = u_c[sel]
                        log_pi1_y_vec[i], log_pi1_x_vec[i] = cu_c[sel], cx
                        stage1_override[i] = log_sum_c - log_sum_r
                    else:
                        stage1_override[i] = -np.inf
                Xp[i] = Xp_temp
                phys_new[i] = phys_temp
                continue

            # Choose block
            block = layout.blocks[rng.randint(0, len(layout.blocks))]
            unit = layout.units[block.unit]
            Nc = unit.Nc

            # Scaling factors and jump partner indices
            # MATLAB: D = randsample(1:n_pairs, 1, 'true');
            D = int(rng.randint(1, n_pairs + 1))

            # MATLAB:
            # a = R(i,draw(1:D,i));
            # b = R(i,draw(D+1:2*D,i));
            if in_prod:
                partners = groups[group_of[i]][groups[group_of[i]] != i]
                pick = partners[rng.permutation(len(partners))[:2 * D]]
                a_idx, b_idx = pick[:D], pick[D:2 * D]
            else:
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
                    if in_prod:
                        other_chains = groups[group_of[i]][groups[group_of[i]] != i]
                    else:
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
            if block.kind in ("nC", "eta"):
                # Continuous-to-discrete mapping with triangle-fold
                unit.decode_discrete(Xp_temp, phys_temp)

            Xp[i] = Xp_temp
            phys_new[i] = phys_temp

        # ── Delayed-acceptance ─────────────────────────────────────────────
        log_alpha1 = beta_arr * (log_pi1_y_vec - log_pi1_x_vec) + J
        exact_rows = ~np.isnan(stage1_override)
        log_alpha1[exact_rows] = stage1_override[exact_rows]
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

        # ── Parallel tempering swaps after the burn-in ─────────────────────
        if in_prod and (t % swap_every_sampling) == 0:
            start = (t // swap_every_sampling) % 2
            for lev in range(start, n_hot, 2):
                i_s = int(rng.choice(cold_idx)) if lev == 0 else n_cold + lev - 1
                j_s = n_cold + lev
                Delta = (beta_arr[i_s] - beta_arr[j_s]) * (p_X[j_s] - p_X[i_s])
                swap_try[lev] += 1
                if Delta >= 0 or np.log(rng.random_sample()) < Delta:
                    X[[i_s, j_s]] = X[[j_s, i_s]]
                    p_X[[i_s, j_s]] = p_X[[j_s, i_s]]
                    phys[[i_s, j_s]] = phys[[j_s, i_s]]
                    swap_acc[lev] += 1

        # ── Outlier detection during burn-in ──────────────────────────────
        if not burnin_done and t > t_burnin and (t - t_burnin) % every == 0:
            omega = p_x[t - every:t].mean(axis=0)
            thr = omega.max() - margin_best
            outliers = [j for j in range(N_chains) if omega[j] < thr]
            good = np.array([j for j in range(N_chains) if omega[j] >= thr])
            if outliers:
                # Reset to good chains drawn at random
                donors = good[rng.randint(0, len(good), len(outliers))]
                for j_out, j_don in zip(outliers, donors):
                    X[j_out] = X[j_don].copy()
                    p_X[j_out] = p_X[j_don]
                    phys[j_out] = phys[j_don].copy()
                n_resets += len(outliers)
                if fid is not None:
                    fid.write(f"Chains {outliers} below the others (mean log-post < {thr:.1f}) "
                              f"at t={t}. Reset to random good chains {donors.tolist()} "
                              f"(good: {good.tolist()}).\n")
            elif t >= 2 * t_burnin:
                burnin_done = True
                t_burn = t
                if fid is not None:
                    fid.write(f"No outlier chains at t={t}: burn-in ends here.\n")
            elif fid is not None:
                fid.write(f"No outlier chains at t={t}; burn-in continues to at least {2 * t_burnin}.\n")
            if not burnin_done and t >= burnin_cap * max_iterations:
                burnin_done = True
                t_burn = t
                if fid is not None:
                    fid.write(f"WARNING: burn-in cap reached at t={t} with outlier chains still present.\n")

        # ── Store current state ────────────────────────────────────────────
        phys_chain[t] = phys.T
        p_x[t] = p_X

        if t % 10 == 0:
            x_archive[t] = X.T

        # ── R-hat convergence check ────────────────────────────────────────
        # Parameters with structurally constant eta (e.g. n-paraffins): NaN in
        # storage (gap in plot), trivially converged in the test.
        if burnin_done and t > t_burn and (t - t_burn) % 2 == 0:
            rc = _r_hat_classic(phys_chain[:, :, :n_cold], t_burn, t)
            rc[_rhat_skip] = np.nan
            R_hat[counter_check] = rc
            counter_check += 1

        if t % 100 == 0:
            diag = ""
            if burnin_done and t - t_burn >= 20:
                rc = _r_hat_classic(phys_chain[:, :, :n_cold], t_burn, t)
                rc[_rhat_skip] = np.nan
                rc_t = np.where(np.isnan(rc), 1.0, rc)
                rc_hist[t] = np.where(np.isfinite(rc_t), rc_t, 1e18).max()
                diag = f"  R-hat max={np.nanmax(rc):.3f}"
                ts = [t - j * rhat_spacing for j in range(rhat_checks)]
                stop = (t - t_burn >= rhat_min_samples
                        and all(tt in rc_hist and rc_hist[tt] <= rhat_stop for tt in ts))
                if stop:
                    win = phys_chain[t_burn:t + 1, :, :n_cold]
                    ess_now = min(_ess_bulk(win[:, c, :]) for c in _ess_cols)
                    diag += f"  ESS min={ess_now:.0f}"
                    stop = ess_now >= ess_min
                    if info is not None:
                        info["ESS_min_at_stop"] = ess_now
                if stop:
                    t_convergence = t
                    convergence = True
                    if fid is not None:
                        fid.write(f"R-hat <= {rhat_stop} at t = {ts[::-1]} and effective sample size "
                                  f">= {ess_min}: convergence after {t_convergence} iterations.\n")
                    print(f"\r  DE-MC iteration {t}/{max_iterations}  AR={AR[t].mean():.1f}%{diag}"
                          f"  -> converged", flush=True)
                    break
            print(f"\r  DE-MC iteration {t}/{max_iterations}  AR={AR[t].mean():.1f}%{diag}", end="", flush=True)

    print()  # newline after progress

    if fid is not None:
        if not convergence:
            fid.write(f"Convergence not reached. Using all {max_iterations} samples.\n")
        fid.close()

    shutdown_pool()

    # ── Trim arrays to t_convergence ─────────────────────────────────────
    tc = t_convergence
    keep = np.arange(n_cold)
    x_chain = phys_chain[:tc][:, :, keep]

    p_x = p_x[:tc][:, keep]
    AR = AR[:tc][:, keep]
    R_hat = R_hat[:counter_check]

    # ── Post-burnin samples ───────────────────────────────────────────────
    if not burnin_done:
        t_burn = tc - 1
    chain_reshaped = np.vstack(
        [x_chain[t_burn:tc, :, ch] for ch in range(len(keep))]
    )
    posterior_reshaped = np.concatenate(
        [p_x[t_burn:tc, ch] for ch in range(len(keep))]
    )
    if info is not None:
        info["t_burnin"] = t_burn
        info["n_resets"] = n_resets
        info["swap_rates"] = swap_acc / np.maximum(swap_try, 1)
        info["betas"] = beta_prod

    return (
        x_chain,
        p_x,
        chain_reshaped,
        posterior_reshaped,
        AR,
        R_hat,
        t_convergence,
    )
