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
The flash module returns the left-hand side of the Rachford-Rice flash
distillation equation for the given surrogate mixture at a fixed temperature,
recovered mole fraction, liquid-phase molar composition, and pressure.
Vapour pressures of all components are evaluated via Yaws' polynomial
correlations and used to compute the equilibrium ratios K_i = p_sat,i / p.

Auxiliary parameters:
Nc: number of surrogate mixture components  [-]

Inputs:
1) T          : temperature  [K]
2) beta       : recovered mole fraction (vapour fraction)  [-]
3) classes    : (1 x Nc) list of lists of Species — candidate species per
                surrogate component
4) index_n_eta: (1 x Nc) array, with the i-th element representing the index
                of the species in the i-th list of classes displaying the
                closest topochemical atom index to the index currently
                investigated for the i-th surrogate mixture component
5) X          : (1 x Nc-1) array containing the molar fractions of the
                surrogate components in the liquid phase
6) pressure   : pressure the distillation curve is computed at  [Pa]

Outputs:
1) F: left-hand side of the Rachford-Rice flash distillation equation
------------------------------------------------------------------------
"""

from __future__ import annotations

import numpy as np

from bayesaf.thermo_transport.hydrocarbons import Species
from bayesaf.thermo_transport.properties import liquid_property


def flash_residual(
    T: float,
    beta: float,
    classes: list[list[Species]],
    index_n_eta: list[int] | np.ndarray,
    X: np.ndarray,
    pressure: float,
) -> float:
    """
    Rachford-Rice flash equation residual at fixed *T* and vapour fraction *beta*.

    Parameters
    ----------
    T : float
        Temperature [K].
    beta : float
        Vapour fraction (mol/mol).
    classes : list of list[Species]
        Candidate species for each surrogate component.
    index_n_eta : array-like of int, shape (Nc,)
        Species index per component.
    X : ndarray, shape (Nc-1,)
        Molar fractions of the first Nc-1 liquid components.
    pressure : float
        Pressure [Pa].

    Returns
    -------
    float
        Residual of the Rachford-Rice equation (should be zero at equilibrium).
    """
    Nc = len(classes)
    psat = np.empty(Nc, dtype=float)
    for i in range(Nc):
        sp = classes[i][index_n_eta[i]]
        psat[i] = liquid_property("vaporPressure", T, sp, pressure)

    K = psat / pressure

    X_full = np.empty(Nc, dtype=float)
    X_full[:Nc - 1] = X
    X_full[Nc - 1] = 1.0 - X.sum()

    return float(np.sum(X_full * (K - 1.0) / (1.0 + beta * (K - 1.0))))
