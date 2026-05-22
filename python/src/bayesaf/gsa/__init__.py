"""
    ____                  _____ ___    ______
   / __ )____ ___  _____ / ___//   |  / ____/
  / __  / __ `/ / / / _ \\__ \\/ /| | / /_
 / /_/ / /_/ / /_/ /  __/__/ / ___ |/ __/
/_____/\\__,_\\__, /\\___/____/_/  |_/_/
            /____/

BayeSAF: Emulation and Design of Sustainable Alternative Fuels
via Bayesian Inference and Descriptors-Based Machine Learning

© 2026 Jacopo Liberatori, Davide Cavalieri, Matteo Blandino
Reference: https://doi.org/10.1016/j.fuel.2026.138835

Global sensitivity analysis via grouped first-order Sobol' indices.
"""
from .sobol import sobol_indices, sobol_distillation
from .sobol_lumped import sobol_lumped
__all__ = ["sobol_indices", "sobol_distillation", "sobol_lumped"]
