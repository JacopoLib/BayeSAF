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
The map_finder module returns the molar fractions, numbers of carbon atoms,
and topochemical atom indices characterising the maximum a posteriori (MAP)
surrogate. It identifies the post-burnin sample with the highest log-
posterior value and extracts the corresponding surrogate parameters.

Auxiliary parameters:
Nc           : number of surrogate mixture components  [-]
numIterations: number of iterations for each chain in the DE-MC algorithm
t_convergence: number of iterations at which convergence is declared
t_burnin     : number of burn-in iterations discarded
N_chains     : number of chains in the DE-MC algorithm

Inputs:
1) posterior_pdf  : (N_post,) array of log-posterior values for each
                    post-burnin sample
2) chain_reshaped : (N_post x 3*Nc-1) array of post-burnin samples
                    (molar fractions, nC, eta_B_star_norm)
3) num_components : number of surrogate mixture components  [-]

Outputs:
1) x_MAP          : molar fractions of the MAP surrogate components
2) nc_MAP         : numbers of carbon atoms of the MAP surrogate components
3) eta_B_star_MAP : topochemical atom indices of the MAP surrogate components
------------------------------------------------------------------------
"""

from __future__ import annotations

import numpy as np


def find_map(
    posterior_pdf: np.ndarray,
    chain_reshaped: np.ndarray,
    num_components: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Return the MAP surrogate parameters from the post-burnin chain.

    Parameters
    ----------
    posterior_pdf : ndarray, shape (N_post,)
        Log-posterior values for each post-burnin sample.
    chain_reshaped : ndarray, shape (N_post, 3*Nc-1)
        Post-burnin samples (molar fractions, nC, eta_B_star_norm).
    num_components : int
        Number of surrogate mixture components.

    Returns
    -------
    x_MAP : ndarray, shape (Nc,)
        Molar fractions of the MAP surrogate.
    nc_MAP : ndarray, shape (Nc,)
        Carbon-atom counts of the MAP surrogate.
    eta_B_star_MAP : ndarray, shape (Nc,)
        Normalised topochemical atom indices of the MAP surrogate.
    """
    Nc = num_components
    idx = int(np.argmax(posterior_pdf))

    x_partial = chain_reshaped[idx, :Nc - 1]
    x_MAP = np.append(x_partial, 1.0 - x_partial.sum())
    nc_MAP = chain_reshaped[idx, Nc - 1 : 2 * Nc - 1]
    eta_B_star_MAP = chain_reshaped[idx, 2 * Nc - 1 : 3 * Nc - 1]

    return x_MAP, nc_MAP, eta_B_star_MAP
