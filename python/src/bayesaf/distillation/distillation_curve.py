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
The distillation_curve module computes the atmospheric (or at a specified
pressure) distillation curve of the surrogate mixture. At each evaporation
step (0 % to 100 % recovered volume) the bubble temperature is found by
solving the Rachford-Rice equation (Brent method, equivalent to MATLAB fzero)
for a vanishingly small vapour fraction beta. Species data are cached once
before the stepping loop to avoid repeated look-ups (optimisation present in
DistillationCurve_fzero.m). Batch variants operate over multiple samples in
parallel using shared memory and multiprocessing.

Inputs:
1) X0         : (1 x Nc-1) array containing the initial liquid-phase molar
                fractions of the first Nc-1 surrogate components
2) classes    : (1 x Nc) list of lists of Species — candidate species per
                surrogate component
3) index_n_eta: (1 x Nc) array with the species index for each component
4) pressure   : pressure at which the distillation curve is computed  [Pa]

Outputs:
1) volume_fraction  : volume fraction evaporated  [%], shape (101,)
2) T_distillation   : boiling temperature at each step  [K], shape (101,)
------------------------------------------------------------------------
"""

from __future__ import annotations

import math
import multiprocessing as mp
from multiprocessing.shared_memory import SharedMemory
from typing import Sequence

import numpy as np
from scipy.optimize import brentq

from bayesaf.thermo_transport.hydrocarbons import Species
from bayesaf.thermo_transport.properties import liquid_property
from .flash import flash_residual

# ---------------------------------------------------------------------------
# Optional Numba JIT backend
# ---------------------------------------------------------------------------
# When Numba is available, brentq and the Rachford-Rice residual are fused
# into a single compiled function (_find_T_nb).  This eliminates the thousands of
# Python-level callbacks that scipy.brentq would otherwise make per curve,
# yielding a further 3-5× speedup on top of the pre-extracted coefficient
# optimisation.  Falls back silently to pure-numpy if Numba is absent.

try:
    import numba as _nb
    _NUMBA_AVAILABLE = True
except ImportError:
    _NUMBA_AVAILABLE = False
    _nb = None


if _NUMBA_AVAILABLE:

    @_nb.njit(cache=True)
    def _rr_residual_nb(T, beta, pAp, pBp, pCp, pDp, pEp, X_full, pressure):
        """Rachford-Rice residual — scalar loop, no numpy overhead."""
        lt = math.log10(T)
        s = 0.0
        _eps = 1e-12
        for i in range(len(X_full)):
            psat_i = (10.0 ** (pAp[i] + pBp[i] / T + pCp[i] * lt
                               + pDp[i] * T + pEp[i] * T * T) * 133.322387415)
            K_i = psat_i / pressure
            den = 1.0 + beta * (K_i - 1.0)
            if den >= 0.0 and den < _eps:
                den = _eps
            elif den < 0.0 and den > -_eps:
                den = -_eps
            s += X_full[i] * (K_i - 1.0) / den
        return s

    @_nb.njit(cache=True)
    def _find_T_nb(T_guess, beta, pAp, pBp, pCp, pDp, pEp, X, pressure):
        """
        Full root finder (bracket widening + Brent's) compiled to native code.

        Parameters
        ----------
        X : 1-D array, shape (Nc-1,)  — liquid molar fractions (last omitted)
        """
        # Build X_full = [X | 1 - sum(X)]
        Nc = len(X) + 1
        X_full = np.empty(Nc)
        X_sum = 0.0
        for i in range(Nc - 1):
            X_full[i] = X[i]
            X_sum += X[i]
        X_full[Nc - 1] = 1.0 - X_sum

        # ── Bracket widening (mirrors _find_root in Python path) ─────────────
        T_lo = max(200.0, T_guess - 100.0)
        T_hi = min(1200.0, T_guess + 100.0)
        fa = _rr_residual_nb(T_lo, beta, pAp, pBp, pCp, pDp, pEp, X_full, pressure)
        fb = _rr_residual_nb(T_hi, beta, pAp, pBp, pCp, pDp, pEp, X_full, pressure)
        if fa * fb > 0.0:
            T_lo = max(150.0, T_lo - 100.0)
            T_hi = min(1500.0, T_hi + 100.0)
            fa = _rr_residual_nb(T_lo, beta, pAp, pBp, pCp, pDp, pEp, X_full, pressure)
            fb = _rr_residual_nb(T_hi, beta, pAp, pBp, pCp, pDp, pEp, X_full, pressure)
            if fa * fb > 0.0:
                T_lo = max(150.0, T_lo - 100.0)
                T_hi = min(1500.0, T_hi + 100.0)
                fa = _rr_residual_nb(T_lo, beta, pAp, pBp, pCp, pDp, pEp, X_full, pressure)
                fb = _rr_residual_nb(T_hi, beta, pAp, pBp, pCp, pDp, pEp, X_full, pressure)
                if fa * fb > 0.0:
                    return T_guess          # fallback: no bracket found

        # ── Brent's method (xtol = 0.01 K) ──────────────────────────────────
        xtol = 0.01

        # Ensure b = better endpoint (smaller |f|)
        if abs(fa) < abs(fb):
            T_lo, T_hi = T_hi, T_lo
            fa, fb = fb, fa

        T_c, fc = T_lo, fa
        mflag = True
        T_d = T_lo      # previous T_c (only used when mflag is False)

        for _ in range(50):
            if abs(T_hi - T_lo) < xtol or fb == 0.0:
                break

            if fa != fc and fb != fc:
                # Inverse quadratic interpolation
                T_s = (T_lo * fb * fc / ((fa - fb) * (fa - fc))
                       + T_hi * fa * fc / ((fb - fa) * (fb - fc))
                       + T_c  * fa * fb / ((fc - fa) * (fc - fb)))
            else:
                # Secant method
                T_s = T_hi - fb * (T_hi - T_lo) / (fb - fa)

            # Conditions to fall back to bisection
            mid3 = (3.0 * T_lo + T_hi) / 4.0
            lo3, hi3 = (mid3, T_hi) if mid3 < T_hi else (T_hi, mid3)
            c1 = T_s < lo3 or T_s > hi3
            c2 = mflag       and abs(T_s - T_hi) >= 0.5 * abs(T_hi - T_c)
            c3 = (not mflag) and abs(T_s - T_hi) >= 0.5 * abs(T_c - T_d)
            c4 = mflag       and abs(T_hi - T_c) < xtol
            c5 = (not mflag) and abs(T_c - T_d) < xtol

            if c1 or c2 or c3 or c4 or c5:
                T_s = (T_lo + T_hi) * 0.5
                mflag = True
            else:
                mflag = False

            fs = _rr_residual_nb(T_s, beta, pAp, pBp, pCp, pDp, pEp, X_full, pressure)
            T_d = T_c
            T_c, fc = T_hi, fb

            if fa * fs < 0.0:
                T_hi, fb = T_s, fs
            else:
                T_lo, fa = T_s, fs

            if abs(fa) < abs(fb):
                T_lo, T_hi = T_hi, T_lo
                fa, fb = fb, fa

        return T_hi


# ---------------------------------------------------------------------------
# Persistent worker pool
# ---------------------------------------------------------------------------
# One module-level pool that lives for the duration of the program.
# Call init_pool() once before the MCMC loop and shutdown_pool() afterwards.

_pool: mp.pool.Pool | None = None
_pool_n_workers: int | None = None
_pool_has_classes: bool = False

# Module-level storage for classes shared with workers via pool initializer.
_POOL_CLASSES: list | None = None

# ---------------------------------------------------------------------------
# Shared-memory state
# ---------------------------------------------------------------------------
# When init_pool() is called with batch_size and nc, four SharedMemory blocks
# are allocated: input (X_all, idx_all) and output (vol_frac, T_dist).
# Per-iteration I/O then reduces to a fast memcpy on the main side and a
# direct array read/write on the worker side — no pickling of array data.
#
# Workers receive only (start, end, pressure) per chunk instead of the full
# numpy slices, cutting per-task pickle payload and
# reducing the number of IPC round-trips from N_samples to n_workers.

# Main-process handles
_shm_X:   SharedMemory | None = None   # X_all         float64 (N, Nc-1)
_shm_idx: SharedMemory | None = None   # idx_all        int64   (N, Nc)
_shm_vf:  SharedMemory | None = None   # vol_frac_all   float64 (N, 101)
_shm_T:   SharedMemory | None = None   # T_dist_all     float64 (N, 101)
_shm_X_arr:   np.ndarray | None = None
_shm_idx_arr: np.ndarray | None = None
_shm_vf_arr:  np.ndarray | None = None
_shm_T_arr:   np.ndarray | None = None
_shm_batch_cap: int | None = None      # allocated row capacity

# Worker-side handles (set via pool initializer)
_WORKER_SHM_X:   SharedMemory | None = None
_WORKER_SHM_IDX: SharedMemory | None = None
_WORKER_SHM_VF:  SharedMemory | None = None
_WORKER_SHM_T:   SharedMemory | None = None
_WORKER_X_ARR:   np.ndarray | None = None
_WORKER_IDX_ARR: np.ndarray | None = None
_WORKER_VF_ARR:  np.ndarray | None = None
_WORKER_T_ARR:   np.ndarray | None = None


def _pool_initializer(classes: list) -> None:
    """Called once per worker process to pre-load the species database."""
    global _POOL_CLASSES
    _POOL_CLASSES = classes


def _pool_initializer_shm(
    classes: list,
    shm_X_name: str,   X_shape: tuple,
    shm_idx_name: str, idx_shape: tuple,
    shm_vf_name: str,  vf_shape: tuple,
    shm_T_name: str,   T_shape: tuple,
) -> None:
    """Pool initializer: pre-loads classes AND attaches to shared memory."""
    global _POOL_CLASSES
    global _WORKER_SHM_X, _WORKER_SHM_IDX, _WORKER_SHM_VF, _WORKER_SHM_T
    global _WORKER_X_ARR, _WORKER_IDX_ARR, _WORKER_VF_ARR, _WORKER_T_ARR

    _POOL_CLASSES = classes
    _WORKER_SHM_X   = SharedMemory(name=shm_X_name)
    _WORKER_SHM_IDX = SharedMemory(name=shm_idx_name)
    _WORKER_SHM_VF  = SharedMemory(name=shm_vf_name)
    _WORKER_SHM_T   = SharedMemory(name=shm_T_name)
    _WORKER_X_ARR   = np.ndarray(X_shape,   dtype=np.float64, buffer=_WORKER_SHM_X.buf)
    _WORKER_IDX_ARR = np.ndarray(idx_shape, dtype=np.int64,   buffer=_WORKER_SHM_IDX.buf)
    _WORKER_VF_ARR  = np.ndarray(vf_shape,  dtype=np.float64, buffer=_WORKER_SHM_VF.buf)
    _WORKER_T_ARR   = np.ndarray(T_shape,   dtype=np.float64, buffer=_WORKER_SHM_T.buf)


def _teardown_shared_memory() -> None:
    """Release shared memory blocks (main process)."""
    global _shm_X, _shm_idx, _shm_vf, _shm_T
    global _shm_X_arr, _shm_idx_arr, _shm_vf_arr, _shm_T_arr, _shm_batch_cap
    for shm in (_shm_X, _shm_idx, _shm_vf, _shm_T):
        if shm is not None:
            try:
                shm.close()
                shm.unlink()
            except Exception:
                pass
    _shm_X = _shm_idx = _shm_vf = _shm_T = None
    _shm_X_arr = _shm_idx_arr = _shm_vf_arr = _shm_T_arr = None
    _shm_batch_cap = None


def init_pool(
    n_workers: int | None = None,
    classes: list | None = None,
    batch_size: int | None = None,
    nc: int | None = None,
) -> mp.pool.Pool:
    """Start (or return) the persistent worker pool.

    Parameters
    ----------
    n_workers : int, optional
        Number of worker processes.  Defaults to ``os.cpu_count()``.
    classes : list of list[Species], optional
        Species databases.  Pre-loaded into each worker to avoid per-task
        pickling of the large species database.
    batch_size : int, optional
        Number of MCMC samples per batch (= N_chains).  When provided
        together with *nc*, shared memory blocks are allocated for zero-copy
        I/O between the main process and worker processes.
    nc : int, optional
        Number of surrogate components.  Required alongside *batch_size*.
    """
    global _pool, _pool_n_workers, _pool_has_classes
    global _shm_X, _shm_idx, _shm_vf, _shm_T
    global _shm_X_arr, _shm_idx_arr, _shm_vf_arr, _shm_T_arr, _shm_batch_cap

    n = n_workers or mp.cpu_count()
    use_shm = (batch_size is not None and nc is not None and classes is not None)
    rebuild = _pool is None or _pool_n_workers != n

    if rebuild:
        if _pool is not None:
            _pool.terminate()
            _teardown_shared_memory()

        if use_shm:
            N, Nc = batch_size, nc
            X_shape   = (N, Nc - 1)
            idx_shape = (N, Nc)
            vf_shape  = (N, 101)
            T_shape   = (N, 101)
            _shm_X   = SharedMemory(create=True, size=int(np.prod(X_shape))   * 8)
            _shm_idx = SharedMemory(create=True, size=int(np.prod(idx_shape)) * 8)
            _shm_vf  = SharedMemory(create=True, size=int(np.prod(vf_shape))  * 8)
            _shm_T   = SharedMemory(create=True, size=int(np.prod(T_shape))   * 8)
            _shm_X_arr   = np.ndarray(X_shape,   dtype=np.float64, buffer=_shm_X.buf)
            _shm_idx_arr = np.ndarray(idx_shape, dtype=np.int64,   buffer=_shm_idx.buf)
            _shm_vf_arr  = np.ndarray(vf_shape,  dtype=np.float64, buffer=_shm_vf.buf)
            _shm_T_arr   = np.ndarray(T_shape,   dtype=np.float64, buffer=_shm_T.buf)
            _shm_batch_cap = N
            _pool = mp.Pool(
                processes=n,
                initializer=_pool_initializer_shm,
                initargs=(
                    classes,
                    _shm_X.name,   X_shape,
                    _shm_idx.name, idx_shape,
                    _shm_vf.name,  vf_shape,
                    _shm_T.name,   T_shape,
                ),
            )
            _pool_has_classes = True
        elif classes is not None:
            _pool = mp.Pool(
                processes=n,
                initializer=_pool_initializer,
                initargs=(classes,),
            )
            _pool_has_classes = True
        else:
            _pool = mp.Pool(processes=n)
            _pool_has_classes = False
        _pool_n_workers = n
    return _pool


def shutdown_pool() -> None:
    """Terminate the persistent worker pool and release shared memory."""
    global _pool, _pool_has_classes, _pool_n_workers
    if _pool is not None:
        _pool.terminate()
        _pool = None
        _pool_n_workers = None
        _pool_has_classes = False
    _teardown_shared_memory()


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _find_root(
    T_guess: float,
    beta: float,
    classes: list[list[Species]],
    index_n_eta: np.ndarray,
    X: np.ndarray,
    pressure: float,
) -> float:
    """Bracket-adaptive brentq wrapper (mirrors find_root_fzero in MATLAB)."""
    T_low = max(200.0, T_guess - 100.0)
    T_high = min(1200.0, T_guess + 100.0)
    F = lambda T: flash_residual(T, beta, classes, index_n_eta, X, pressure)

    for _ in range(3):
        f_low, f_high = F(T_low), F(T_high)
        if np.sign(f_low) != np.sign(f_high):
            break
        T_low = max(150.0, T_low - 100.0)
        T_high = min(1500.0, T_high + 100.0)
    else:
        return T_guess  # fallback

    try:
        return brentq(F, T_low, T_high, xtol=1e-6)
    except Exception:
        return T_guess


def _eval_rho(
    T: float,
    classes: list[list[Species]],
    index_n_eta: np.ndarray,
    pressure: float,
) -> np.ndarray:
    Nc = len(classes)
    out = np.empty(Nc, dtype=float)
    for j in range(Nc):
        out[j] = liquid_property("rho", T, classes[j][index_n_eta[j]], pressure)
    return out


def _eval_psat(
    T: float,
    classes: list[list[Species]],
    index_n_eta: np.ndarray,
    pressure: float,
) -> np.ndarray:
    Nc = len(classes)
    out = np.empty(Nc, dtype=float)
    for j in range(Nc):
        out[j] = liquid_property("vaporPressure", T, classes[j][index_n_eta[j]], pressure)
    return out


# ---------------------------------------------------------------------------
# Single distillation curve
# ---------------------------------------------------------------------------

def distillation_curve(
    X0: np.ndarray,
    classes: list[list[Species]],
    index_n_eta: Sequence[int] | np.ndarray,
    pressure: float = 101325.0,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute the distillation curve for a single surrogate mixture.

    Parameters
    ----------
    X0 : ndarray, shape (Nc-1,)
        Molar fractions of the first Nc-1 components.
    classes : list of list[Species]
        Species lists for each component.
    index_n_eta : array-like of int, shape (Nc,)
        Species index per component.
    pressure : float
        Pressure [Pa].

    Returns
    -------
    volume_fraction : ndarray, shape (101,)
        Volume fractions evaporated (0–100 in steps of 1).
    T_distillation : ndarray, shape (101,)
        Bubble/dew temperatures at each volume fraction [K].
    """
    Nc = len(classes)
    index_n_eta = np.asarray(index_n_eta, dtype=int)
    X0 = np.asarray(X0, dtype=float)

    n_steps = 101

    T_dist = np.zeros(n_steps)
    vol_frac = np.zeros(n_steps)

    # ── Pre-extract all species data ONCE (avoids repeated attribute lookups
    #    inside the brentq inner loop, which calls flash_residual ~30×/step) ──
    W        = np.array([classes[j][index_n_eta[j]].mol_weight for j in range(Nc)], dtype=float)
    psat_c   = np.array([classes[j][index_n_eta[j]].coeff_psat for j in range(Nc)], dtype=float)  # (Nc, 5)
    rho_c    = np.array([classes[j][index_n_eta[j]].coeff_rho  for j in range(Nc)], dtype=float)  # (Nc, 3)
    Tc_arr   = np.array([classes[j][index_n_eta[j]].Tc          for j in range(Nc)], dtype=float)

    # Unpack coefficient columns for fully vectorised operations
    pAp, pBp, pCp, pDp, pEp = psat_c[:, 0], psat_c[:, 1], psat_c[:, 2], psat_c[:, 3], psat_c[:, 4]
    rA, rB, rC = rho_c[:, 0], rho_c[:, 1], rho_c[:, 2]

    # Reusable full-composition buffer (avoids per-call allocation inside brentq)
    _Xbuf = np.empty(Nc, dtype=float)

    def _psat_v(T: float) -> np.ndarray:
        """Vapour pressure [Pa] for all Nc species at temperature T (scalar)."""
        lt = math.log10(T)
        return 10.0 ** (pAp + pBp / T + pCp * lt + pDp * T + pEp * T * T) * 133.322387415

    def _rho_v(T: float) -> np.ndarray:
        """Liquid density [kg/m³] for all Nc species at temperature T."""
        tau = np.maximum(0.0, 1.0 - T / Tc_arr)
        return 1000.0 * rA * rB ** (-(tau) ** rC)

    def _flash_v(T: float, beta: float, X: np.ndarray) -> float:
        """Rachford-Rice residual (vectorised, no attribute lookups)."""
        _Xbuf[:Nc - 1] = X
        _Xbuf[Nc - 1] = 1.0 - X.sum()
        K = _psat_v(T) / pressure
        den = 1.0 + beta * (K - 1.0)
        eps = 1e-12
        den = np.where(np.abs(den) < eps, np.where(den >= 0, eps, -eps), den)
        return float(np.dot(_Xbuf, (K - 1.0) / den))

    def _root_v(T_guess: float, beta: float, X: np.ndarray) -> float:
        """brentq root finder using the fast flash residual."""
        T_lo = max(200.0, T_guess - 100.0)
        T_hi = min(1200.0, T_guess + 100.0)
        F = lambda T: _flash_v(T, beta, X)
        for _ in range(3):
            if math.copysign(1, F(T_lo)) != math.copysign(1, F(T_hi)):
                break
            T_lo = max(150.0, T_lo - 100.0)
            T_hi = min(1500.0, T_hi + 100.0)
        else:
            return T_guess
        try:
            # xtol=0.01 K: more than sufficient for distillation curves.
            return brentq(F, T_lo, T_hi, xtol=0.01)
        except Exception:
            return T_guess

    # ── Step 0: initial temperature (bubble point at beta → 0) ──────────────
    beta = 1e-6
    X = X0.copy()
    if _NUMBA_AVAILABLE:
        temp = _find_T_nb(300.0, beta, pAp, pBp, pCp, pDp, pEp, X, pressure)
    else:
        temp = _root_v(300.0, beta, X)
    T_dist[0] = temp

    rho0 = _rho_v(temp)
    vol_i = X0 * W[:Nc - 1] / rho0[:Nc - 1]
    V0 = vol_i.sum() + (1.0 - X0.sum()) * W[Nc - 1] / rho0[Nc - 1]

    # ── Steps 1–100 ──────────────────────────────────────────────────────────
    if _NUMBA_AVAILABLE:
        # Hot path: brentq + flash fused in native code — no Python callbacks.
        for ii in range(n_steps - 1):
            beta = 1.0 / (100.0 - ii)
            temp = _find_T_nb(temp, beta, pAp, pBp, pCp, pDp, pEp, X, pressure)
            T_dist[ii + 1] = temp

            rho_i = _rho_v(temp)
            K = _psat_v(temp) / pressure

            den = 1.0 + beta * (K[:Nc - 1] - 1.0)
            X = X / den
            Y = K[:Nc - 1] * X
            Y_last = 1.0 - Y.sum()

            Vtilde = (Y * W[:Nc - 1] / rho_i[:Nc - 1]).sum() + Y_last * W[Nc - 1] / rho_i[Nc - 1]
            vol_frac[ii + 1] = vol_frac[ii] + Vtilde / V0 / (n_steps - 1)
    else:
        for ii in range(n_steps - 1):
            beta = 1.0 / (100.0 - ii)
            temp = _root_v(temp, beta, X)
            T_dist[ii + 1] = temp

            rho_i = _rho_v(temp)
            K = _psat_v(temp) / pressure

            den = 1.0 + beta * (K[:Nc - 1] - 1.0)
            X = X / den                          # update liquid composition (Nc-1)
            Y = K[:Nc - 1] * X
            Y_last = 1.0 - Y.sum()

            Vtilde = (Y * W[:Nc - 1] / rho_i[:Nc - 1]).sum() + Y_last * W[Nc - 1] / rho_i[Nc - 1]
            # Each of the 100 flash steps evaporates exactly 1/100 of the initial
            # moles (β_k × remaining = 1/100 for all k).  Vtilde/V0 ≈ 1 per step
            # (molar volumes are similar for all species), so without the 1/100
            # normalisation the sum would reach ~100 before the final ×100 scaling,
            # giving a vf range of ~10 000 instead of 0–100.
            vol_frac[ii + 1] = vol_frac[ii] + Vtilde / V0 / (n_steps - 1)

    # Rescale to 0–100
    vol_frac = vol_frac * 100.0

    return vol_frac, T_dist


# ---------------------------------------------------------------------------
# Batch version (parallelised over samples)
# ---------------------------------------------------------------------------

def _worker(args):
    X_row, classes, idx_row, pressure = args
    return distillation_curve(X_row, classes, idx_row, pressure)


def _worker_shared(args):
    """Lightweight worker: uses classes pre-loaded via pool initializer."""
    X_row, idx_row, pressure = args
    return distillation_curve(X_row, _POOL_CLASSES, idx_row, pressure)


def _worker_chunk_shm(args):
    """Process samples [start:end] using shared memory; return None (results in shm)."""
    start, end, pressure = args
    for s in range(start, end):
        X_row   = _WORKER_X_ARR[s].copy()
        idx_row = _WORKER_IDX_ARR[s].copy()
        vf, T   = distillation_curve(X_row, _POOL_CLASSES, idx_row, pressure)
        _WORKER_VF_ARR[s] = vf
        _WORKER_T_ARR[s]  = T


def distillation_curve_batch(
    X_all: np.ndarray,
    classes: list[list[Species]],
    index_n_eta_all: np.ndarray,
    pressure: float = 101325.0,
    n_workers: int | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute distillation curves for multiple MCMC samples in parallel.

    Parameters
    ----------
    X_all : ndarray, shape (N_samples, Nc-1)
    classes : list of list[Species]
    index_n_eta_all : ndarray of int, shape (N_samples, Nc)
    pressure : float
    n_workers : int, optional
        Number of worker processes.  Defaults to ``os.cpu_count()``.

    Returns
    -------
    volume_fraction_all : ndarray, shape (N_samples, 101)
    T_distillation_all  : ndarray, shape (N_samples, 101)
    """
    N_samples = X_all.shape[0]

    # ── Shared memory fast path ───────────────────────────────────────────────
    if (_pool is not None and n_workers is None
            and _shm_X is not None and N_samples <= _shm_batch_cap):
        # Zero-copy input: write proposals directly into shared memory.
        np.copyto(_shm_X_arr[:N_samples], X_all)
        _shm_idx_arr[:N_samples] = index_n_eta_all   # handles int32→int64 cast
        # Coarse tasks: one chunk per worker → n_workers round-trips instead of N_samples.
        n_w = _pool_n_workers
        chunk = max(1, (N_samples + n_w - 1) // n_w)
        tasks = [
            (i * chunk, min((i + 1) * chunk, N_samples), pressure)
            for i in range(n_w) if i * chunk < N_samples
        ]
        _pool.map(_worker_chunk_shm, tasks)
        # Zero-copy output: read results directly from shared memory.
        return _shm_vf_arr[:N_samples].copy(), _shm_T_arr[:N_samples].copy()

    # ── Fallback paths ────────────────────────────────────────────────────────
    args = [(X_all[s], classes, index_n_eta_all[s], pressure) for s in range(N_samples)]
    if _pool is not None and n_workers is None:
        if _pool_has_classes:
            slim_args = [(X_all[s], index_n_eta_all[s], pressure) for s in range(N_samples)]
            results = _pool.map(_worker_shared, slim_args)
        else:
            results = _pool.map(_worker, args)
    elif N_samples >= (n_workers or mp.cpu_count()):
        with mp.Pool(processes=n_workers) as pool:
            results = pool.map(_worker, args)
    else:
        results = [_worker(a) for a in args]

    vol_all = np.array([r[0] for r in results])
    T_all   = np.array([r[1] for r in results])
    return vol_all, T_all


# ---------------------------------------------------------------------------
# Cheap single-sample variant (early termination at vol1)
# ---------------------------------------------------------------------------

def distillation_curve_cheap(
    X0: np.ndarray,
    classes: list[list[Species]],
    index_n_eta: Sequence[int] | np.ndarray,
    pressure: float = 101325.0,
    vol1: float = 10.0,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Like :func:`distillation_curve` but stops once ``volume_fraction >= vol1``.

    Remaining entries are NaN (for compatibility with the full-curve API).

    Parameters
    ----------
    vol1 : float
        Volume-fraction threshold for early termination (0–100 scale).
    """
    Nc = len(classes)
    index_n_eta = np.asarray(index_n_eta, dtype=int)
    X0 = np.asarray(X0, dtype=float)

    n_steps = 101
    T_dist = np.full(n_steps, np.nan)
    vol_frac = np.arange(n_steps, dtype=float)  # 0…100

    # ── Pre-extract species data once (same optimisation as distillation_curve)
    W        = np.array([classes[j][index_n_eta[j]].mol_weight for j in range(Nc)], dtype=float)
    psat_c   = np.array([classes[j][index_n_eta[j]].coeff_psat for j in range(Nc)], dtype=float)
    rho_c    = np.array([classes[j][index_n_eta[j]].coeff_rho  for j in range(Nc)], dtype=float)
    Tc_arr   = np.array([classes[j][index_n_eta[j]].Tc          for j in range(Nc)], dtype=float)

    pAp, pBp, pCp, pDp, pEp = psat_c[:, 0], psat_c[:, 1], psat_c[:, 2], psat_c[:, 3], psat_c[:, 4]
    rA, rB, rC = rho_c[:, 0], rho_c[:, 1], rho_c[:, 2]
    _Xbuf = np.empty(Nc, dtype=float)

    def _psat_v(T: float) -> np.ndarray:
        lt = math.log10(T)
        return 10.0 ** (pAp + pBp / T + pCp * lt + pDp * T + pEp * T * T) * 133.322387415

    def _rho_v(T: float) -> np.ndarray:
        return 1000.0 * rA * rB ** (-(1.0 - T / Tc_arr) ** rC)

    def _flash_v(T: float, beta: float, X: np.ndarray) -> float:
        _Xbuf[:Nc - 1] = X
        _Xbuf[Nc - 1] = 1.0 - X.sum()
        K = _psat_v(T) / pressure
        den = 1.0 + beta * (K - 1.0)
        eps = 1e-12
        den = np.where(np.abs(den) < eps, np.where(den >= 0, eps, -eps), den)
        return float(np.dot(_Xbuf, (K - 1.0) / den))

    def _root_v(T_guess: float, beta: float, X: np.ndarray) -> float:
        T_lo = max(200.0, T_guess - 100.0)
        T_hi = min(1200.0, T_guess + 100.0)
        F = lambda T: _flash_v(T, beta, X)
        for _ in range(3):
            if math.copysign(1, F(T_lo)) != math.copysign(1, F(T_hi)):
                break
            T_lo = max(150.0, T_lo - 100.0)
            T_hi = min(1500.0, T_hi + 100.0)
        else:
            return T_guess
        try:
            return brentq(F, T_lo, T_hi, xtol=0.01)
        except Exception:
            return T_guess

    beta = 1e-6
    X = X0.copy()
    if _NUMBA_AVAILABLE:
        temp = _find_T_nb(300.0, beta, pAp, pBp, pCp, pDp, pEp, X, pressure)
    else:
        temp = _root_v(300.0, beta, X)
    T_dist[0] = temp

    rho0 = _rho_v(temp)
    V0 = (X0 * W[:Nc - 1] / rho0[:Nc - 1]).sum() + (1.0 - X0.sum()) * W[Nc - 1] / rho0[Nc - 1]

    for ii in range(n_steps - 1):
        beta = 1.0 / (100.0 - ii)
        if _NUMBA_AVAILABLE:
            temp = _find_T_nb(temp, beta, pAp, pBp, pCp, pDp, pEp, X, pressure)
        else:
            temp = _root_v(temp, beta, X)
        T_dist[ii + 1] = temp

        rho_i = _rho_v(temp)
        K = _psat_v(temp) / pressure

        den = 1.0 + beta * (K[:Nc - 1] - 1.0)
        X = X / den
        Y = K[:Nc - 1] * X
        Y_last = 1.0 - Y.sum()

        Vtilde = (Y * W[:Nc - 1] / rho_i[:Nc - 1]).sum() + Y_last * W[Nc - 1] / rho_i[Nc - 1]
        # Same 1/100 normalisation as distillation_curve; since this cheap
        # version has no final ×100 scaling, Vtilde/V0/(n_steps-1)×100 = Vtilde/V0.
        vol_frac[ii + 1] = vol_frac[ii] + Vtilde / V0

        if vol_frac[ii + 1] >= vol1:
            T_dist[ii + 2:] = np.nan
            break

    return vol_frac, T_dist


def _worker_cheap_shared(args):
    """Lightweight cheap worker: uses classes pre-loaded via pool initializer."""
    X_row, idx_row, pressure, vol1 = args
    return distillation_curve_cheap(X_row, _POOL_CLASSES, idx_row, pressure, vol1)


def _worker_cheap(args):
    X_row, classes, idx_row, pressure, vol1 = args
    return distillation_curve_cheap(X_row, classes, idx_row, pressure, vol1)


def _worker_chunk_shm_cheap(args):
    """Cheap variant of _worker_chunk_shm (early termination at vol1)."""
    start, end, pressure, vol1 = args
    for s in range(start, end):
        X_row   = _WORKER_X_ARR[s].copy()
        idx_row = _WORKER_IDX_ARR[s].copy()
        vf, T   = distillation_curve_cheap(X_row, _POOL_CLASSES, idx_row, pressure, vol1)
        _WORKER_VF_ARR[s] = vf
        _WORKER_T_ARR[s]  = T


def distillation_curve_batch_cheap(
    X_all: np.ndarray,
    classes: list[list[Species]],
    index_n_eta_all: np.ndarray,
    pressure: float = 101325.0,
    vol1: float = 10.0,
    n_workers: int | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Parallelised cheap distillation curves (early termination at *vol1*)."""
    N_samples = X_all.shape[0]

    # ── Shared memory fast path ───────────────────────────────────────────────
    if (_pool is not None and n_workers is None
            and _shm_X is not None and N_samples <= _shm_batch_cap):
        np.copyto(_shm_X_arr[:N_samples], X_all)
        _shm_idx_arr[:N_samples] = index_n_eta_all
        n_w = _pool_n_workers
        chunk = max(1, (N_samples + n_w - 1) // n_w)
        tasks = [
            (i * chunk, min((i + 1) * chunk, N_samples), pressure, vol1)
            for i in range(n_w) if i * chunk < N_samples
        ]
        _pool.map(_worker_chunk_shm_cheap, tasks)
        return _shm_vf_arr[:N_samples].copy(), _shm_T_arr[:N_samples].copy()

    # ── Fallback paths ────────────────────────────────────────────────────────
    args = [(X_all[s], classes, index_n_eta_all[s], pressure, vol1) for s in range(N_samples)]
    if _pool is not None and n_workers is None:
        if _pool_has_classes:
            slim_args = [(X_all[s], index_n_eta_all[s], pressure, vol1) for s in range(N_samples)]
            results = _pool.map(_worker_cheap_shared, slim_args)
        else:
            results = _pool.map(_worker_cheap, args)
    elif N_samples >= (n_workers or mp.cpu_count()):
        with mp.Pool(processes=n_workers) as pool:
            results = pool.map(_worker_cheap, args)
    else:
        results = [_worker_cheap(a) for a in args]

    vol_all = np.array([r[0] for r in results])
    T_all   = np.array([r[1] for r in results])
    return vol_all, T_all
