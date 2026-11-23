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

Post-processing: MAP estimation, posterior pushforward,
and visualisation of Bayesian inference results.
"""
from .map_finder import find_map
from .map_writer import write_map
from .props_pushforward import props_pushforward
from .visualize import visualize_results
__all__ = [
    "find_map",
    "write_map",
    "props_pushforward",
    "visualize_results",
]
