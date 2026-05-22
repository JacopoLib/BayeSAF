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

Differential Evolution Markov Chain (DE-MC) sampler with parallel
tempering and delayed acceptance.
"""
from .demc import run_demc
__all__ = ["run_demc"]
