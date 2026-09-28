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
The layout module describes how the DE-MC state vector is organised, so that
the sampler never hard-codes column indices.

A state vector is made of one or more *composition units*. A unit is one
surrogate mixture of Nc components, described by (Nc-1) molar fractions, Nc
numbers of carbon atoms, and Nc normalized topochemical atom indices. Each
unit exposes three proposal *blocks* (fractions, nC, eta). The sampler keeps
two arrays with the same column layout:
- X    : normalised coordinates in [0,1], where the DE jumps are made;
- phys : physical values (molar fractions, nC, eta_B_star_norm), which are
         what the posterior is evaluated on and what is stored in the chain.

SurrogateLayout is the standard BayeSAF case: a single unit whose columns are
[mol(Nc-1) | nC(Nc) | eta(Nc)], i.e. 3*Nc-1 parameters. Other models (e.g. a
hierarchical one with a unit per fuel) provide their own layout and
posterior-argument convention.
------------------------------------------------------------------------
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Sequence

import numpy as np
from scipy.spatial.distance import cdist

from bayesaf.thermo_transport.hydrocarbons import Species


# ---------------------------------------------------------------------------
# Parameter encoding/decoding
# ---------------------------------------------------------------------------

def _decode_nc(u: float, n_range: list[int]) -> int:
    """Map normalised u ∈ [0,1] to a discrete carbon count."""
    k = len(n_range)
    idx = min(int(np.floor(u * k)), k - 1)
    return n_range[idx]


def _decode_eta(
    u: float,
    sp_list: list[Species],
    nC_val: int,
) -> float:
    """Map normalised u ∈ [0,1] to the closest eta_B_star_norm for nC_val.

    Replicates the MATLAB formula exactly:
        idx_eta = round(1 + u*(N-1))   [1-based, round-half-away-from-zero]
    converted to 0-based Python indexing:
        idx_eta = floor(1 + u*(N-1) + 0.5) - 1
    """
    eta_all = np.array([sp.eta_B_star for sp in sp_list])
    eta_sorted = np.sort(eta_all)
    N = len(eta_sorted)
    # MATLAB-compatible rounding: floor(x + 0.5) = round-half-away-from-zero
    idx_eta = int(np.floor(1.0 + u * (N - 1) + 0.5)) - 1
    idx_eta = max(0, min(idx_eta, N - 1))
    eta_target = eta_sorted[idx_eta]

    # Restrict to species with matching nC
    nC_arr = np.array([sp.nC for sp in sp_list], dtype=int)
    eta_norm_arr = np.array([sp.eta_B_star_norm for sp in sp_list])
    match = np.where(nC_arr == nC_val)[0]
    eta_match = eta_all[match]
    eta_norm_match = eta_norm_arr[match]
    return float(eta_norm_match[np.argmin(np.abs(eta_match - eta_target))])


def _triangle_fold(x: float) -> float:
    """Reflect x into [0,1] with a triangle-wave mapping."""
    return 1.0 - abs(1.0 - (x % 2))


# ---------------------------------------------------------------------------
# Composition unit and proposal blocks
# ---------------------------------------------------------------------------

@dataclass
class CompositionUnit:
    """One surrogate mixture inside the state vector.

    frac_cols, nc_cols and eta_cols are the columns of X/phys holding the
    (Nc-1) molar fractions, Nc carbon numbers and Nc eta_B_star_norm values.
    """

    frac_cols: np.ndarray
    nc_cols: np.ndarray
    eta_cols: np.ndarray
    classes: list[list[Species]]
    n_ranges: list[list[int]]
    lower_bound_x: np.ndarray
    upper_bound_x: np.ndarray

    @property
    def Nc(self) -> int:
        return len(self.nc_cols)

    def fold_fractions(self, Xp_row: np.ndarray, phys_row: np.ndarray) -> None:
        """Boundary handling for the molar fractions; updates both rows in place."""
        for j, c in enumerate(self.frac_cols):
            if Xp_row[c] < 0:
                Xp_row[c] = max(1 - abs(Xp_row[c]), 0)
            elif Xp_row[c] > 1:
                Xp_row[c] = min(abs(Xp_row[c] - 1), 1)
            phys_row[c] = (Xp_row[c] * (self.upper_bound_x[j] - self.lower_bound_x[j])
                           + self.lower_bound_x[j])

    def fractions_sum_ok(self, phys_row: np.ndarray) -> bool:
        """True if the implied last molar fraction lies within its bounds."""
        s = phys_row[self.frac_cols].sum()
        return s <= 1.0 - self.lower_bound_x[-1] and s >= 1.0 - self.upper_bound_x[-1]

    def decode_discrete(self, Xp_row: np.ndarray, phys_row: np.ndarray) -> None:
        """Triangle-fold the nC/eta coordinates and decode them in place."""
        for k, c in enumerate(self.nc_cols):
            xj = _triangle_fold(Xp_row[c])
            Xp_row[c] = xj
            phys_row[c] = _decode_nc(xj, self.n_ranges[k])
        for k, c in enumerate(self.eta_cols):
            xj = _triangle_fold(Xp_row[c])
            Xp_row[c] = xj
            phys_row[c] = _decode_eta(xj, self.classes[k], int(phys_row[self.nc_cols[k]]))

    def encode_discrete(self, Xp_row: np.ndarray, k: int, nC: float, eta: float) -> None:
        """Write normalised coordinates for component k after a Gibbs move."""
        n_rng = self.n_ranges[k]
        Xp_row[self.nc_cols[k]] = (nC - min(n_rng)) / (max(n_rng) - min(n_rng) + 1e-12)

        sp_list = self.classes[k]
        eta_sorted = np.sort(np.array([sp.eta_B_star for sp in sp_list]))
        eta_list_nc = np.array([sp.eta_B_star for sp in sp_list if sp.nC == nC])
        if len(eta_list_nc) > 0:
            eta_unnorm = min(eta_list_nc) + eta * (max(eta_list_nc) - min(eta_list_nc))
            idx_neighb = int(np.argmin(np.abs(eta_sorted - eta_unnorm)))
        else:
            idx_neighb = 0
        Xp_row[self.eta_cols[k]] = idx_neighb / max(1, len(eta_sorted) - 1)

    def decode(self, X_rows: np.ndarray, phys_rows: np.ndarray) -> None:
        """Decode normalised rows into physical rows (no folding)."""
        for r in range(X_rows.shape[0]):
            for j, c in enumerate(self.frac_cols):
                phys_rows[r, c] = (X_rows[r, c] * (self.upper_bound_x[j] - self.lower_bound_x[j])
                                   + self.lower_bound_x[j])
            for k, c in enumerate(self.nc_cols):
                phys_rows[r, c] = _decode_nc(X_rows[r, c], self.n_ranges[k])
            for k, c in enumerate(self.eta_cols):
                phys_rows[r, c] = _decode_eta(X_rows[r, c], self.classes[k],
                                              int(phys_rows[r, self.nc_cols[k]]))

    def constant_eta(self) -> np.ndarray:
        """Per component: True if all species share one eta_B_star_norm value."""
        return np.array([len({sp.eta_B_star_norm for sp in cls}) == 1 for cls in self.classes])


@dataclass(frozen=True)
class Block:
    """A group of columns updated together by one DE proposal."""

    kind: str          # "fractions", "nC" or "eta"
    unit: int          # index into layout.units
    cols: np.ndarray


# ---------------------------------------------------------------------------
# Layouts
# ---------------------------------------------------------------------------

class ParameterLayout:
    """Base class: units, blocks and the posterior-argument convention.

    Subclasses must set ``units``, ``blocks`` and ``n_params`` and implement
    ``posterior_args`` and ``init_population``.
    """

    units: list[CompositionUnit]
    blocks: list[Block]
    n_params: int

    def posterior_args(self, phys_rows: np.ndarray) -> tuple:
        """Arguments passed to the posterior callables for a batch of rows."""
        raise NotImplementedError

    def init_population(
        self,
        N_chains: int,
        posterior_fn: Callable,
        rng: np.random.RandomState,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return (X, phys, log-posterior) for the initial population."""
        raise NotImplementedError

    def decode_init(
        self,
        x_init: np.ndarray,
        posterior_fn: Callable,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Decode an external normalised initial population (benchmark mode)."""
        X = x_init.copy()
        phys = np.zeros_like(X, dtype=float)
        for unit in self.units:
            unit.decode(X, phys)
        return X, phys, posterior_fn(*self.posterior_args(phys))

    def fractions_block(self, unit: int) -> Block:
        return next(b for b in self.blocks if b.kind == "fractions" and b.unit == unit)

    def rhat_skip(self) -> np.ndarray:
        """Columns that are structurally constant (skipped by the R-hat test)."""
        skip = np.zeros(self.n_params, dtype=bool)
        for unit in self.units:
            skip[unit.eta_cols[unit.constant_eta()]] = True
        return skip


class SurrogateLayout(ParameterLayout):
    """Standard BayeSAF layout: one unit, columns [mol(Nc-1) | nC(Nc) | eta(Nc)]."""

    def __init__(
        self,
        classes: list[list[Species]],
        n_ranges: list[Sequence[int]],
        lower_bound_x: np.ndarray,
        upper_bound_x: np.ndarray,
    ) -> None:
        Nc = len(classes)
        self.n_params = 3 * Nc - 1
        unit = CompositionUnit(
            frac_cols=np.arange(0, Nc - 1),
            nc_cols=np.arange(Nc - 1, 2 * Nc - 1),
            eta_cols=np.arange(2 * Nc - 1, 3 * Nc - 1),
            classes=classes,
            n_ranges=[sorted({sp.nC for sp in cls}) for cls in classes],
            lower_bound_x=np.asarray(lower_bound_x, dtype=float),
            upper_bound_x=np.asarray(upper_bound_x, dtype=float),
        )
        self.units = [unit]
        self.blocks = [
            Block("fractions", 0, unit.frac_cols),
            Block("nC", 0, unit.nc_cols),
            Block("eta", 0, unit.eta_cols),
        ]

    def posterior_args(self, phys_rows: np.ndarray) -> tuple:
        u = self.units[0]
        return phys_rows[:, u.frac_cols], phys_rows[:, u.nc_cols], phys_rows[:, u.eta_cols]

    def init_population(
        self,
        N_chains: int,
        posterior_fn: Callable,
        rng: np.random.RandomState,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Sample an initial population on the parameter space.

        Dirichlet(1) molar fractions within bounds and uniform nC/eta
        coordinates; N_chains diverse points are then picked by greedy
        max-min distance in physical space.
        """
        unit = self.units[0]
        Nc = unit.Nc
        n_params = self.n_params
        N_needed = 10 * N_chains
        batch = max(1000, 5 * N_needed)
        K_rest = 2 * Nc   # nC (Nc) + eta (Nc)
        lb, ub = unit.lower_bound_x, unit.upper_bound_x

        valid_samples = np.empty((0, n_params), dtype=float)

        while valid_samples.shape[0] < N_needed:
            # Dirichlet(1) samples
            Y = -np.log(rng.random_sample((batch, Nc)))
            dir_samples = Y / Y.sum(axis=1, keepdims=True)

            ok = np.all((dir_samples >= lb) & (dir_samples <= ub), axis=1)
            val_dir = dir_samples[ok]
            if val_dir.shape[0] == 0:
                continue

            # Normalise molar fractions to [0,1]
            val_norm = (val_dir[:, :Nc - 1] - lb[:Nc - 1]) / (ub[:Nc - 1] - lb[:Nc - 1])
            n_valid = val_norm.shape[0]
            rest = rng.random_sample((n_valid, K_rest))
            batch_samples = np.hstack([val_norm, rest])
            valid_samples = np.vstack([valid_samples, batch_samples])

        # Decode physical parameters for all valid samples
        X_lhs = valid_samples.copy()
        phys_lhs = np.zeros_like(X_lhs)
        unit.decode(valid_samples, phys_lhs)

        # Compute posterior for all samples
        p_lhs = posterior_fn(*self.posterior_args(phys_lhs))

        # Select N_chains diverse initial points (max-min-distance greedy)
        N = phys_lhs.shape[0]
        selected = np.zeros(N_chains, dtype=int)
        selected[0] = rng.randint(N)
        for k in range(1, N_chains):
            remaining = np.setdiff1d(np.arange(N), selected[:k])
            dists = cdist(phys_lhs[remaining], phys_lhs[selected[:k]]).min(axis=1)
            selected[k] = remaining[dists.argmax()]

        perm = rng.permutation(N_chains)
        selected = selected[perm]

        return X_lhs[selected], phys_lhs[selected], p_lhs[selected]
