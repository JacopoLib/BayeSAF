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

Rachford-Rice flash and distillation curve computation.
"""
from .flash import flash_residual
from .distillation_curve import (
    distillation_curve,
    distillation_curve_cheap,
    distillation_curve_batch,
    distillation_curve_batch_cheap,
    init_pool,
    shutdown_pool,
)
__all__ = [
    "flash_residual",
    "distillation_curve",
    "distillation_curve_cheap",
    "distillation_curve_batch",
    "distillation_curve_batch_cheap",
    "init_pool",
    "shutdown_pool",
]
