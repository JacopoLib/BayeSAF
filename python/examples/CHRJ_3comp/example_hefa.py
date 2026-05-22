r"""
    ____                  _____ ___    ______
   / __ )____ ___  _____ / ___//   |  / ____/
  / __  / __ `/ / / / _ \\__ \/ /| | / /_
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
The present example illustrates how to develop a three-component surrogate
mixture emulating the HEFA-SPK POSF-6152 (CHRJ) fuel.

Usage
-----
    cd python
    pip install -e .
    cd examples/CHRJ_3comp
    python example_hefa.py
"""

from __future__ import annotations

import numpy as np

from bayesaf.thermo_transport.hydrocarbons import load_hydrocarbons
from bayesaf.utilities.misc import p_distillation
from bayesaf.utilities.data_import import data_import
from bayesaf.pdfs.pdf_build import build_pdfs, build_pdfs_cheap
from bayesaf.mcmc.demc import run_demc
from bayesaf.postprocessing.map_finder import find_map
from bayesaf.postprocessing.map_writer import write_map
from bayesaf.postprocessing.props_pushforward import props_pushforward
from bayesaf.postprocessing.visualize import visualize_results

# ---------------------------------------------------------------------------
# Real fuel metadata
# ---------------------------------------------------------------------------
DATASET = ["rho", "nu", "distillation", "deltaT_dist",
          "molWeight", "HC", "DCN", "LHV", "flash", "freezing"]
DATA_DIR = "../../../exp/CHRJ_POSF6152"   # directory containing <name>.csv files
FUEL_NAME = "CHRJ POSF-6152"
T_BUBBLE = 445.00       # bubble temperature of the real fuel [K] (if known, otherwise set it to 1e+18)

# ---------------------------------------------------------------------------
# Surrogate configuration
# ---------------------------------------------------------------------------

NUM_COMPONENTS = 3

FAMILIES = ["nparaffins", "isoparaffins", "isoparaffins"]

N_RANGES = [
    list(range(8, 16)),    # nparaffins  C8–C15
    list(range(9, 13)),    # isoparaffins C9–C12
    list(range(13, 17)),    # isoparaffins C13–C16
]

LOWER_ETA = np.zeros(NUM_COMPONENTS)
UPPER_ETA = np.ones(NUM_COMPONENTS)

LOWER_MOL_FRAC = np.array([0.1, 0.1, 0.1])
UPPER_MOL_FRAC = np.array([0.75, 0.75, 0.75])

ALPHA_DIRICHLET = np.ones(NUM_COMPONENTS)   # symmetric uniform Dirichlet

# ---------------------------------------------------------------------------
# DE-MC parameters
# ---------------------------------------------------------------------------

MAX_ITERATIONS   = 5000
T_BURNIN         = 500
T_CHECK          = 2 * T_BURNIN
SCALING_X        = 0.1
SCALING_NC       = 1.0
SCALING_ETA      = 1.0
NOISE_X          = 1e-6
NOISE_NC         = 1e-6
NOISE_ETA        = 1e-6
N_CHAINS         = 2 * (3 * NUM_COMPONENTS - 1)   # 28 chains
OUTLIER_METHOD   = "mad"
R_HAT_THRESHOLD  = 1.2

# Parallel tempering
PT_SWITCH  = "True"
BETA_MIN   = 1e-4
T_LADDER   = "geometric"
SWAP_FREQ  = 20

# Maximum number of chain pairs for proposal
N_PAIRS  = 1

# Proposal probabilities
P_GIBBS    = 0.1
P_SNOOKER  = 0.1

# Post-processing
BAND_PERCENTILES  = True
CONFIDENCE_WIDTH  = 0.95


def main() -> None:
    # ------------------------------------------------------------------
    # 1. Build species databases (classes)
    # ------------------------------------------------------------------
    print("Loading hydrocarbon databases...")
    classes = [load_hydrocarbons(fam, n_ranges) for fam, n_ranges in zip(FAMILIES, N_RANGES)]

    # ------------------------------------------------------------------
    # 2. Import experimental data
    # ------------------------------------------------------------------
    print("Importing experimental data...")
    pressure_distillation = p_distillation(DATASET, data_dir=DATA_DIR)

    full_data, min_T_array, max_T_array, variable_names = data_import(
        DATASET, FAMILIES, N_RANGES, T_BUBBLE, NUM_COMPONENTS, pressure_distillation,
        data_dir=DATA_DIR,
    )

    # ------------------------------------------------------------------
    # 3. Build posterior PDF handles
    # ------------------------------------------------------------------
    print("Building PDF handles...")
    prior_fn, likelihood_fn, posterior_fn = build_pdfs(
        full_data, FAMILIES, classes, variable_names, pressure_distillation,
        LOWER_MOL_FRAC, UPPER_MOL_FRAC, ALPHA_DIRICHLET, N_RANGES,
        LOWER_ETA, UPPER_ETA,
    )

    # Cheap posterior (excludes deltaT_dist)
    idx_dist = [i for i, v in enumerate(variable_names) if v == "deltaT_dist"]
    full_data_nd    = [d for i, d in enumerate(full_data)      if i not in idx_dist]
    variable_names_nd = [v for i, v in enumerate(variable_names) if i not in idx_dist]
    _, _, posterior_cheap_fn = build_pdfs_cheap(
        full_data_nd, FAMILIES, classes, variable_names_nd, pressure_distillation,
        LOWER_MOL_FRAC, UPPER_MOL_FRAC, ALPHA_DIRICHLET, N_RANGES,
        LOWER_ETA, UPPER_ETA,
    )
    
    # ------------------------------------------------------------------
    # 4. Run DE-MC
    # ------------------------------------------------------------------
    print("Running DE-MC...")
    (chain, posterior_pdf, chain_reshaped, posterior_reshaped,
     AR, R_hat, t_convergence) = run_demc(
        families=FAMILIES,
        posterior_fn=posterior_fn,
        posterior_cheap_fn=posterior_cheap_fn,
        classes=classes,
        lower_bound_x=LOWER_MOL_FRAC,
        upper_bound_x=UPPER_MOL_FRAC,
        n_ranges=N_RANGES,
        lower_bound_eta=LOWER_ETA,
        upper_bound_eta=UPPER_ETA,
        max_iterations=MAX_ITERATIONS,
        t_burnin=T_BURNIN,
        t_check=T_CHECK,
        scaling_factor_x=SCALING_X,
        scaling_factor_nc=SCALING_NC,
        scaling_factor_eta=SCALING_ETA,
        noise_x=NOISE_X,
        noise_nc=NOISE_NC,
        noise_eta=NOISE_ETA,
        N_chains=N_CHAINS,
        outlier_method=OUTLIER_METHOD,
        R_hat_threshold=R_HAT_THRESHOLD,
        PT_switch=PT_SWITCH,
        beta_min=BETA_MIN,
        T_ladder=T_LADDER,
        swap_freq=SWAP_FREQ,
        n_pairs=N_PAIRS,
        p_gibbs=P_GIBBS,
        p_snooker=P_SNOOKER,
        log_file="demc_log.txt",
        seed=42,
    )

    converged = t_convergence < MAX_ITERATIONS
    print(f"\nDE-MC {'converged' if converged else 'did NOT converge'} at t = {t_convergence}.")
    
    if R_hat.size > 0:

        rhat_last = R_hat[-1]
    
        # remove NaN entries
        rhat_valid = rhat_last[np.isfinite(rhat_last)]
    
        if rhat_valid.size == 0:
            print("    → R-hat not available (all parameters constant).")
    
        elif np.all(rhat_valid <= R_HAT_THRESHOLD):
            print(f"    → R-hat ≤ {R_HAT_THRESHOLD} for all non-constant parameters: convergence criterion satisfied.")
    
        else:
            print(f"    → R-hat > {R_HAT_THRESHOLD}: chains not yet converged. Increase MAX_ITERATIONS.")
    
    else:
    
        print("  R-hat = not computed (chain shorter than t_check — increase MAX_ITERATIONS)")
        print(f"  mean AR   = {AR[np.isfinite(AR)].mean():.3f} %")

    # ------------------------------------------------------------------
    # 5. MAP estimate
    # ------------------------------------------------------------------
    print("\nFinding MAP surrogate...")
    x_MAP, nc_MAP, eta_B_star_MAP = find_map(posterior_reshaped, chain_reshaped, NUM_COMPONENTS)

    write_map(
        FAMILIES, NUM_COMPONENTS,
        x_MAP, nc_MAP, eta_B_star_MAP,
        FUEL_NAME, output_file="MAP.txt",
    )

    # ------------------------------------------------------------------
    # 6. Pushforward & visualisation
    # ------------------------------------------------------------------
    print("\nComputing property pushforward distributions...")
    pf = props_pushforward(
        full_data, FAMILIES, classes, chain_reshaped,
        NUM_COMPONENTS, variable_names,
        x_MAP, nc_MAP, eta_B_star_MAP,
        min_T_array, max_T_array,
        pressure_distillation,
        confidence_width=CONFIDENCE_WIDTH,
    )

    print("\nGenerating figures...")
    visualize_results(
        full_data, FAMILIES, classes,
        chain, AR, R_hat, chain_reshaped,
        T_BURNIN, R_HAT_THRESHOLD, N_RANGES, NUM_COMPONENTS, variable_names,
        pf,
        FUEL_NAME,
        confidence_width=CONFIDENCE_WIDTH,
        band_percentiles=BAND_PERCENTILES,
        output_dir=".",
    )

    print("\nDone. Results written to MAP.txt and ./Figures/.")


if __name__ == "__main__":
    main()
