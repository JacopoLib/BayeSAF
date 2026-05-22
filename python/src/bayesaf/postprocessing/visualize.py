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
The visualize module provides tools to analyze and visualize results from
the Bayesian inference analysis. It generates MCMC diagnostic plots (chain
traces, acceptance rates, R-hat convergence) and posterior PDF plots
(marginal distributions of molar fractions, numbers of carbon atoms, and
topochemical atom indices), as well as property pushforward plots showing
MAP predictions and posterior confidence bands versus experimental data.
All figures are saved as 600-dpi PNGs via matplotlib.

Auxiliary parameters:
Nd           : number of experimental measurements per property  [-]
Np           : number of thermophysical properties targeted  [-]
N_chains     : number of DE-MC chains
maxIterations: maximum number of iterations for each chain in DE-MC

Inputs:
1) fullData              : (1 x Np) list of (Nd x 3) or (Nd x 4) arrays of
                           experimental measurements
2) families              : (1 x numComponents) list of hydrocarbon family
                           strings
3) classes               : (1 x numComponents) list of lists of Species
4) chain                 : (t_convergence x 3*numComponents-1 x N_chains)
                           full posterior sample array
5) AR                    : (t_convergence x N_chains) acceptance rate array
6) R_hat                 : R-hat statistic array,
                           shape (floor((t_convergence-t_burnin)/2) x
                           3*numComponents-1)
7) chain_reshaped        : ((t_convergence - t_burnin)*N_chains x
                           3*numComponents-1) post-burnin sample array
8) t_convergence         : iteration at which convergence was declared
9) t_burnin              : number of burn-in iterations discarded
10) n_ranges             : (1 x numComponents) carbon-atom ranges per
                           component
11) numComponents        : number of surrogate mixture components  [-]
12) variable_names       : (1 x Np) list of property name strings
13) x_MAP                : (1 x numComponents) MAP molar fractions
14) nc_MAP               : (1 x numComponents) MAP numbers of carbon atoms
15) eta_B_star_MAP       : (1 x numComponents) MAP topochemical indices
16) minT_array           : (1 x Np) minimum temperatures per property
17) maxT_array           : (1 x Np) maximum temperatures per property
18) pressure_distillation: pressure for distillation curve evaluation  [Pa]
19) fuel_name            : legend string for real fuel experimental data
20) band_percentiles     : whether to colour the confidence band with
                           percentiles ('True' or 'False')
21) confidence_width     : width of the confidence interval
------------------------------------------------------------------------
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")   # non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.stats import gaussian_kde

from bayesaf.thermo_transport.hydrocarbons import Species

# ---------------------------------------------------------------------------
# Global style — match MATLAB (Times, fontsize 20, box on, white bg, 600 dpi)
# ---------------------------------------------------------------------------

matplotlib.rcParams.update({
    "font.family":       "serif",
    "font.serif":        ["Times New Roman", "Times", "DejaVu Serif"],
    "mathtext.fontset":  "stix",        # Times-like math (≈ MATLAB LaTeX)
    "figure.facecolor":  "white",
    "axes.facecolor":    "white",
    "savefig.facecolor": "white",
})

_FS = 20          # font size  (MATLAB: fontsize = 20)
_FN = "serif"     # font name  (MATLAB: fontname = 'Times')
_LW = 1.3         # axis line width
_PIX = (8, 5.8)   # figure size in inches (≈ MATLAB pix × pix/1.38)
_BLUE   = (0.0039, 0.451,  0.741)
_GREEN  = (0.45,   0.80,   0.69)
_ORANGE = (0.8500, 0.3250, 0.0980)
_GRAY   = (0.7,    0.7,    0.7)
_RED    = (0.647,  0.090,  0.184)


def _style(ax: plt.Axes) -> None:
    """Apply MATLAB-matching axes style: box on, white bg, Times font."""
    ax.tick_params(labelsize=_FS)
    for spine in ax.spines.values():
        spine.set_linewidth(_LW)
        spine.set_visible(True)     # all four sides  (MATLAB: box on)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontfamily(_FN)
    ax.set_facecolor("white")


def _save(folder: str | Path, name: str, dpi: int = 600) -> None:
    p = Path(folder)
    p.mkdir(parents=True, exist_ok=True)
    plt.savefig(p / f"{name}.png", dpi=dpi, bbox_inches="tight",
                facecolor="white")
    plt.close("all")


# ---------------------------------------------------------------------------
# MCMC diagnostics
# ---------------------------------------------------------------------------

def plot_traces(
    chain: np.ndarray,      # (t_conv, 3*Nc-1, N_chains)
    num_components: int,
    output_dir: str = ".",
) -> None:
    """Trace plots for all chain parameters."""
    Nc = num_components
    N_chains = chain.shape[2]
    t_conv = chain.shape[0]
    colors = plt.cm.viridis(np.linspace(0, 1, N_chains))
    folder = Path(output_dir) / "Figures" / "MCMC_check"

    # Molar fractions
    for m in range(Nc - 1):
        fig, ax = plt.subplots(figsize=_PIX)
        for i in range(N_chains):
            ax.plot(chain[:, m, i], lw=0.75, color=colors[i])
        ax.set_xlabel("MCMC Iterations", fontsize=_FS, fontfamily=_FN)
        ax.set_ylabel(f"Trace of $X_{m+1}$", fontsize=_FS, fontfamily=_FN)
        _style(ax)
        _save(folder, f"Trace_X{m+1}")

    # nC
    for m in range(Nc):
        fig, ax = plt.subplots(figsize=_PIX)
        for i in range(N_chains):
            ax.plot(chain[:, Nc - 1 + m, i], lw=0.75, color=colors[i])
        ax.set_xlabel("MCMC Iterations", fontsize=_FS, fontfamily=_FN)
        ax.set_ylabel(f"$n_{{C,{m+1}}}$", fontsize=_FS, fontfamily=_FN)
        _style(ax)
        _save(folder, f"Trace_nc{m+1}")

    # eta
    for m in range(Nc):
        fig, ax = plt.subplots(figsize=_PIX)
        for i in range(N_chains):
            ax.plot(chain[:, 2 * Nc - 1 + m, i], lw=0.75, color=colors[i])
        ax.set_xlabel("MCMC Iterations", fontsize=_FS, fontfamily=_FN)
        ax.set_ylabel(f"Trace of $\\eta_{{B,{m+1}}}^{{*}}$", fontsize=_FS, fontfamily=_FN)
        _style(ax)
        _save(folder, f"Trace_eta{m+1}")


# ---------------------------------------------------------------------------
# Posterior PDFs
# ---------------------------------------------------------------------------

def plot_posterior_pdfs(
    chain_reshaped: np.ndarray,   # (N_post, 3*Nc-1)
    num_components: int,
    n_ranges: list[list[int]],
    output_dir: str = ".",
) -> None:
    """Marginal and joint PDFs for all parameter groups."""
    Nc = num_components
    folder_marg = Path(output_dir) / "Figures" / "PosteriorPDF" / "MarginalProbability"
    folder_joint = Path(output_dir) / "Figures" / "PosteriorPDF" / "Joint2D"
    colors = plt.cm.tab20(np.linspace(0, 1, max(Nc + 1, 2)))

    # --- Molar fractions: marginal KDE ---
    for m in range(Nc - 1):
        data_m = chain_reshaped[:, m]
        kde = gaussian_kde(data_m)
        xv = np.linspace(data_m.min(), data_m.max(), 200)
        fig, ax = plt.subplots(figsize=_PIX)
        ax.plot(xv, kde(xv), lw=2, color=_BLUE)
        ax.set_xlabel(f"$X_{m+1}$", fontsize=_FS, fontfamily=_FN)
        ax.set_ylabel("Marginal Probability", fontsize=_FS, fontfamily=_FN)
        xlo, xhi = ax.get_xlim()
        if xlo < 0:
            ax.set_xlim(0, xhi)
        _style(ax)
        _save(folder_marg, f"Marginal_X{m+1}")

    # --- Molar fractions: joint 2D KDE ---
    for m in range(Nc - 1):
        for p in range(m + 1, Nc - 1):
            dm = chain_reshaped[:, m]
            dp = chain_reshaped[:, p]
            xg = np.linspace(dm.min(), dm.max(), 100)
            yg = np.linspace(dp.min(), dp.max(), 100)
            XX, YY = np.meshgrid(xg, yg)
            pts = np.vstack([XX.ravel(), YY.ravel()])
            kde2 = gaussian_kde(np.vstack([dm, dp]))
            ZZ = kde2(pts).reshape(XX.shape)
            fig, ax = plt.subplots(figsize=_PIX)
            cf = ax.contourf(XX, YY, ZZ / ZZ.max(), levels=256, cmap="viridis_r")
            cb = plt.colorbar(cf, ax=ax)
            cb.set_label("Normalized Joint Probability", fontsize=_FS, fontfamily=_FN)
            cb.ax.tick_params(labelsize=_FS)
            ax.set_xlabel(f"$X_{m+1}$", fontsize=_FS, fontfamily=_FN)
            ax.set_ylabel(f"$X_{p+1}$", fontsize=_FS, fontfamily=_FN)
            _style(ax)
            _save(folder_joint, f"Joint2D_X{m+1}_X{p+1}")

    # --- nC: marginal histograms ---
    for m in range(Nc):
        data_m = chain_reshaped[:, Nc - 1 + m]
        nr = n_ranges[m]
        edges = np.arange(min(nr) - 0.5, max(nr) + 1.5)
        counts, _ = np.histogram(data_m, bins=edges, density=False)
        probs = counts / counts.sum()
        fig, ax = plt.subplots(figsize=_PIX)
        ax.bar(np.arange(1, len(nr) + 1), probs, color=colors[m % len(colors)])
        ax.set_xticks(np.arange(1, len(nr) + 1))
        ax.set_xticklabels([str(v) for v in nr], fontsize=_FS - 2)
        ax.set_xlabel(f"$n_{{C,{m+1}}}$", fontsize=_FS, fontfamily=_FN)
        ax.set_ylabel("Marginal Probability", fontsize=_FS, fontfamily=_FN)
        _style(ax)
        _save(folder_marg, f"Marginal_nc{m+1}")

    # --- nC: joint 2D histograms ---
    for m in range(Nc):
        for p in range(m + 1, Nc):
            dm = chain_reshaped[:, Nc - 1 + m]
            dp = chain_reshaped[:, Nc - 1 + p]
            nm, np_ = n_ranges[m], n_ranges[p]
            em = np.arange(min(nm) - 0.5, max(nm) + 1.5)
            ep = np.arange(min(np_) - 0.5, max(np_) + 1.5)
            H, _, _ = np.histogram2d(dm, dp, bins=[em, ep], density=False)
            H = H / H.sum()
            fig, ax = plt.subplots(figsize=_PIX)
            im = ax.imshow(H.T, origin="lower", cmap="binary", aspect="auto")
            cb = plt.colorbar(im, ax=ax)
            cb.set_label("Normalized Joint Probability", fontsize=_FS, fontfamily=_FN)
            cb.ax.tick_params(labelsize=_FS)
            ax.set_xticks(np.arange(len(nm)))
            ax.set_xticklabels([str(v) for v in nm], fontsize=_FS - 3)
            ax.set_yticks(np.arange(len(np_)))
            ax.set_yticklabels([str(v) for v in np_], fontsize=_FS - 3)
            ax.set_xlabel(f"$n_{{C,{m+1}}}$", fontsize=_FS, fontfamily=_FN)
            ax.set_ylabel(f"$n_{{C,{p+1}}}$", fontsize=_FS, fontfamily=_FN)
            _style(ax)
            _save(folder_joint, f"Joint2D_nc{m+1}_nc{p+1}")

    # --- eta: marginal histograms ---
    top_eta = {}
    for m in range(Nc):
        data_m = chain_reshaped[:, 2 * Nc - 1 + m]
        unique_vals = np.unique(data_m)
        # Build edges as midpoints between unique values
        if len(unique_vals) > 1:
            mids = (unique_vals[:-1] + unique_vals[1:]) / 2.0
            edges = np.concatenate([[unique_vals[0] - 1e-10], mids, [unique_vals[-1] + 1e-10]])
        else:
            edges = np.array([unique_vals[0] - 0.5, unique_vals[0] + 0.5])
        counts, _ = np.histogram(data_m, bins=edges, density=False)
        probs = counts / counts.sum()
        top_k = min(10, len(unique_vals))
        top_idx = np.argsort(probs)[-top_k:]
        top_eta[m] = unique_vals[top_idx]

        fig, ax = plt.subplots(figsize=_PIX)
        ax.bar(np.arange(len(unique_vals)), probs, color=colors[(2 * Nc + m) % len(colors)])
        step = max(1, round(len(unique_vals) / 5))
        tick_idx = np.arange(0, len(unique_vals), step)
        ax.set_xticks(tick_idx)
        ax.set_xticklabels([f"{unique_vals[i]:.2f}" for i in tick_idx], fontsize=_FS - 3)
        ax.set_xlabel(f"$\\eta_{{B,{m+1}}}^{{*}}$", fontsize=_FS, fontfamily=_FN)
        ax.set_ylabel("Marginal Probability", fontsize=_FS, fontfamily=_FN)
        _style(ax)
        _save(folder_marg, f"Marginal_eta{m+1}")

    print(f"Posterior PDF figures saved under {folder_marg.parent}.")


# ---------------------------------------------------------------------------
# Property pushforward plots
# ---------------------------------------------------------------------------

def plot_property_pushforward(
    full_data: list[np.ndarray],
    variable_names: list[str],
    temperature_ranges: list,
    property_MAP: list,
    volumeFraction_MAP,
    percentiles_list: np.ndarray,
    mean_model: list,
    percentiles_model: list,
    model_output_samples: list[np.ndarray],
    fuel_name: str,
    confidence_width: float = 0.9,
    band_percentiles: bool = False,
    output_dir: str = ".",
) -> None:
    """
    Plot property pushforward confidence bands and MAP predictions.

    Parameters
    ----------
    band_percentiles : bool
        If True colour each percentile slice individually (colorbar mode);
        if False use a single filled band.
    """
    folder = Path(output_dir) / "Figures" / "PropertyCheck_MAP"
    conf_pct = round(confidence_width * 100)
    pct_low  = 50 * (1 - confidence_width)
    pct_high = 100 - pct_low

    _LUMPED = {"molWeight", "HC", "DCN", "flash", "freezing", "LHV", "deltaT_dist"}
    _YLABELS = {
        "rho": "Density [kg/m$^3$]",
        "mu": "Dynamic Viscosity [$\\mu$Pa s]",
        "nu": "Kinematic Viscosity [mm$^2$/s]",
        "kappa": "Thermal Conductivity [W/m/K]",
        "specificHeat": "Specific Heat Capacity [kJ/kg/K]",
        "latentHeat": "Enthalpy of Vaporization [J/kg]",
        "vaporPressure": "Vapor Pressure [Pa]",
        "sigma": "Surface Tension [N/m]",
        "distillation": "Temperature [$^\\circ$C]",
        "molWeight": "MW [g/mol]",
        "HC": "H/C [-]",
        "DCN": "DCN [-]",
        "flash": "Flash Point [$^\\circ$C]",
        "freezing": "Freezing Point [$^\\circ$C]",
        "LHV": "Lower Heating Value [MJ/kg]",
    }

    n_pct = len(percentiles_list)
    if band_percentiles:
        cmap = cm.get_cmap("viridis", n_pct)

    for h, var in enumerate(variable_names):
        data = full_data[h]
        pct = percentiles_model[h]

        if var not in _LUMPED:
            t_range = temperature_ranges[h]
            prop_map = property_MAP[h]
            is_dist = (var == "distillation")

            fig, ax = plt.subplots(figsize=_PIX)
            xlabel = ("Recovered Volume Fraction [%]"
                      if is_dist else "Temperature [$^\\circ$C]")

            n_t = len(t_range)
            if band_percentiles:
                for i in range(n_t - 1):
                    for j in range(n_pct - 1):
                        ax.fill_between(
                            [t_range[i], t_range[i + 1]],
                            [pct[i, j], pct[i + 1, j]],
                            [pct[i, j + 1], pct[i + 1, j + 1]],
                            color=cmap(j), alpha=1.0,
                        )
                sm = plt.cm.ScalarMappable(cmap="viridis",
                                           norm=plt.Normalize(pct_low, pct_high))
                sm.set_array([])
                cb = plt.colorbar(sm, ax=ax)
                cb.set_label("Percentile [%]", fontsize=_FS, fontfamily=_FN)
                cb.ax.tick_params(labelsize=_FS)
            else:
                # Single band (outermost percentiles)
                lo_band = pct[:, 0]
                hi_band = pct[:, -1]
                ax.fill_between(t_range, lo_band, hi_band,
                                color=_BLUE, alpha=0.2,
                                label=f"{conf_pct}% Confidence Interval")

            ax.plot(t_range, prop_map, lw=2, color="black", label="MAP Surrogate")

            # Experimental data with diamond markers (MATLAB: 'd', MarkerSize 8)
            if is_dist:
                data_y = data[:, 1] - 273.15
                ax.plot(data[:, 0], data_y, "d",
                        ms=8, mfc=_GRAY, mec="k", lw=1, label=fuel_name)
            elif var in ("mu", "nu"):
                data_y = 1e6 * data[:, 1]
                ax.plot(data[:, 0] - 273.15, data_y, "d",
                        ms=8, mfc=_GRAY, mec="k", lw=1, label=fuel_name)
            elif var == "specificHeat":
                data_y = 1e-3 * data[:, 1]
                ax.plot(data[:, 0] - 273.15, data_y, "d",
                        ms=8, mfc=_GRAY, mec="k", lw=1, label=fuel_name)
            else:
                data_y = data[:, 1]
                ax.plot(data[:, 0] - 273.15, data_y, "d",
                        ms=8, mfc=_GRAY, mec="k", lw=1, label=fuel_name)

            # xlim / ylim matching MATLAB
            if is_dist:
                ax.set_xlim(t_range[0], t_range[-1])
                y_lo = 0.99 * min(pct[:, 0].min(), (data[:, 1] - 273.15).min())
                y_hi = 1.01 * max(pct[:, -1].max(), (data[:, 1] - 273.15).max())
                ax.set_ylim(y_lo, y_hi)
            else:
                # xlim: 5% padding on Kelvin scale, converted back to °C
                x_lo = 0.95 * (t_range[0]  + 273.15) - 273.15
                x_hi = 1.05 * (t_range[-1] + 273.15) - 273.15
                ax.set_xlim(x_lo, x_hi)
                y_lo = 0.98 * min(pct[:, 0].min(), data_y.min())
                y_hi = 1.02 * max(pct[:, -1].max(), data_y.max())
                ax.set_ylim(y_lo, y_hi)

            ax.set_xlabel(xlabel, fontsize=_FS, fontfamily=_FN)
            ax.set_ylabel(_YLABELS.get(var, var), fontsize=_FS, fontfamily=_FN)
            ax.legend(fontsize=_FS, loc="best")
            _style(ax)
            _save(folder, var)

        elif var == "deltaT_dist":
            # Two KDE plots: T50-T10 and T90-T10
            # MATLAB style: plain KDE line only (no fill, no MAP line, no exp data)
            # Temperature differences: K ≡ °C numerically — label as °C (MATLAB)
            samples_2d = model_output_samples[h]   # (N, 2)
            for col, name, xlabel in [
                (0, "T5010", "$T_{50-10}$ [$^\\circ$C]"),
                (1, "T9010", "$T_{90-10}$ [$^\\circ$C]"),
            ]:
                col_samples = samples_2d[:, col]
                finite = np.isfinite(col_samples)
                col_disp = col_samples[finite]

                fig, ax = plt.subplots(figsize=_PIX)
                if col_disp.size > 1 and np.std(col_disp) > 0:
                    kde = gaussian_kde(col_disp)
                    x_lo, x_hi = col_disp.min(), col_disp.max()
                    x_grid = np.linspace(x_lo, x_hi, 300)
                    ax.plot(x_grid, kde(x_grid), color=_BLUE, lw=2)
                    ax.set_xlim(x_lo, x_hi)
                ax.set_xlabel(xlabel, fontsize=_FS, fontfamily=_FN)
                ax.set_ylabel("Probability Density", fontsize=_FS, fontfamily=_FN)
                _style(ax)
                _save(folder, name)

        else:
            # Scalar lumped property: posterior KDE
            # MATLAB style: plain blue KDE line only — no fill, no MAP line, no exp
            _OFFSET = {"flash": -273.15, "freezing": -273.15}
            offset = _OFFSET.get(var, 0.0)

            samples = model_output_samples[h]           # (N,)
            finite = np.isfinite(samples)
            samples_disp = samples[finite] + offset

            fig, ax = plt.subplots(figsize=_PIX)
            if samples_disp.size > 1 and np.std(samples_disp) > 0:
                kde = gaussian_kde(samples_disp)
                x_lo, x_hi = samples_disp.min(), samples_disp.max()
                x_grid = np.linspace(x_lo, x_hi, 300)
                ax.plot(x_grid, kde(x_grid), color=_BLUE, lw=2)
                # xlim: 5% padding (MATLAB: xlim([0.95*min, 1.05*max]))
                ax.set_xlim(0.95 * x_lo, 1.05 * x_hi)
            ax.set_xlabel(_YLABELS.get(var, var), fontsize=_FS, fontfamily=_FN)
            ax.set_ylabel("Probability Density", fontsize=_FS, fontfamily=_FN)
            _style(ax)
            _save(folder, var)

    print(f"Property check figures saved under {folder}.")


# ---------------------------------------------------------------------------
# Sobol' index plots
# ---------------------------------------------------------------------------

def plot_sobol(
    variable_names: list[str],
    sobol_idx: list[np.ndarray],
    temperature_ranges: list,
    output_dir: str = ".",
) -> None:
    """Grouped Sobol' first-order index bar charts."""
    folder = Path(output_dir) / "Figures" / "GSA"
    _LUMPED = {"molWeight", "HC", "DCN", "flash", "freezing", "LHV", "deltaT_dist"}
    bar_colors = [_GREEN, _BLUE, _ORANGE]
    labels = ["$X_i$", "$n_{C,i}$", "$\\eta_{B,i}^{*}$"]

    lumped_vars: list[str] = []
    lumped_si: list[np.ndarray] = []

    for h, var in enumerate(variable_names):
        si = sobol_idx[h]  # (3, n_T) or (3,) or (3, 2)

        if var not in _LUMPED:
            t_range = temperature_ranges[h]
            is_dist = (var == "distillation")
            if is_dist:
                selected_idx = [1, 25, 50, 75, 99]
                x_vals = [1, 25, 50, 75, 99]
            else:
                selected_idx = [0, 6, 12, 18, 24]
                x_vals = [round(float(t_range[i]), 0) for i in selected_idx]
            sel_si = si[:, selected_idx]  # (3, 5)

            fig, ax = plt.subplots(figsize=_PIX)
            width = 0.25
            x = np.arange(5)
            for g in range(3):
                ax.bar(x + (g - 1) * width, sel_si[g], width,
                       color=bar_colors[g], edgecolor="black", lw=1,
                       label=labels[g])
            ax.set_xticks(x)
            ax.set_xticklabels([str(v) for v in x_vals], fontsize=_FS - 2)
            xlabel = ("Recovered Volume Fraction [%]"
                      if is_dist else "Temperature [$^\\circ$C]")
            ax.set_xlabel(xlabel, fontsize=_FS, fontfamily=_FN)
            ax.set_ylabel(f"$S_{{i,{var}}}$ [-]", fontsize=_FS, fontfamily=_FN)
            ax.set_ylim(0, 1.05)
            ax.legend(labels, fontsize=_FS, loc="upper center",
                      bbox_to_anchor=(0.5, 1.18), ncol=3)
            fig.subplots_adjust(top=0.82)
            _style(ax)
            _save(folder, f"Si_{var}")

        else:
            if var == "deltaT_dist":
                # si is (3, 2)
                lumped_vars += ["$T_{50-10}$", "$T_{90-10}$"]
                for col in range(2):
                    lumped_si.append(si[:, col])
            else:
                lumped_vars.append(var)
                lumped_si.append(si)  # (3,)

    # Lumped properties bar chart
    if lumped_si:
        si_mat = np.column_stack(lumped_si)  # (3, n_lumped)
        n_lump = si_mat.shape[1]
        fig, ax = plt.subplots(figsize=_PIX)
        width = 0.25
        x = np.arange(n_lump)
        for g in range(3):
            ax.bar(x + (g - 1) * width, si_mat[g], width,
                   color=bar_colors[g], edgecolor="black", lw=1,
                   label=labels[g])
        ax.set_xticks(x)
        ax.set_xticklabels(lumped_vars, fontsize=_FS - 2, rotation=45, ha="right")
        ax.set_ylabel("$S_{i,lumped}$ [-]", fontsize=_FS, fontfamily=_FN)
        ax.set_ylim(0, 1.05)
        ax.legend(labels, fontsize=_FS, loc="upper center",
                  bbox_to_anchor=(0.5, 1.18), ncol=3)
        fig.subplots_adjust(top=0.82)
        _style(ax)
        _save(folder, "Si_lumped")

    print(f"Sobol' index figures saved under {folder}.")


# ---------------------------------------------------------------------------
# Acceptance rate + R-hat
# ---------------------------------------------------------------------------

def plot_mcmc_diagnostics(
    AR: np.ndarray,        # (t_conv, N_chains)
    R_hat: np.ndarray,     # (n_iter, n_params)
    R_hat_threshold: float,
    t_burnin: int,
    num_components: int,
    output_dir: str = ".",
) -> None:
    """Plot acceptance rates and Gelman-Rubin R-hat statistic."""
    folder = Path(output_dir) / "Figures" / "MCMC_check"
    Nc = num_components
    N_chains = AR.shape[1]
    colors = plt.cm.viridis(np.linspace(0, 1, N_chains))

    # Acceptance rate
    fig, ax = plt.subplots(figsize=_PIX)
    for i in range(N_chains):
        ax.plot(AR[:, i], lw=0.75, color=colors[i])
    ax.axvline(t_burnin, color="red", ls="--", lw=1, label=f"Burn-in ({t_burnin})")
    ax.set_xlabel("MCMC Iterations", fontsize=_FS, fontfamily=_FN)
    ax.set_ylabel("Acceptance Rate [%]", fontsize=_FS, fontfamily=_FN)
    ax.legend(fontsize=_FS - 2)
    _style(ax)
    _save(folder, "AcceptanceRate")

    # R-hat (NaN values for discrete parameters plot as gaps — expected behaviour)
    if R_hat.size > 0:
        param_labels = (
            [f"$X_{m+1}$" for m in range(Nc - 1)]
            + [f"$n_{{C,{m+1}}}$" for m in range(Nc)]
            + [f"$\\eta_{{B,{m+1}}}^{{*}}$" for m in range(Nc)]
        )
        fig, ax = plt.subplots(figsize=_PIX)
        colors_p = plt.cm.tab10(np.linspace(0, 1, R_hat.shape[1]))
        # NaN → plot as gaps; any residual sentinels clamped to nan as well
        R_hat_plot = np.where(np.isfinite(R_hat), R_hat, np.nan)
        for p in range(R_hat_plot.shape[1]):
            ax.plot(R_hat_plot[:, p], lw=0.75, color=colors_p[p],
                    label=param_labels[p] if p < len(param_labels) else f"p{p}")
        ax.axhline(R_hat_threshold, color="red", ls="--", lw=1, label=f"R-hat = {R_hat_threshold}")
        finite_vals = R_hat_plot[np.isfinite(R_hat_plot)]
        if finite_vals.size > 0:
            ax.set_ylim(0, min(finite_vals.max() * 1.15, 5.0))
        ax.set_xlabel("MCMC Iteration (post-burnin)", fontsize=_FS, fontfamily=_FN)
        ax.set_ylabel("R-hat", fontsize=_FS, fontfamily=_FN)
        ax.legend(fontsize=_FS - 3, ncol=2, loc="upper right")
        _style(ax)
        _save(folder, "Rhat")
    else:
        print("  R-hat plot skipped (no convergence checks performed).")

    print(f"MCMC diagnostic figures saved under {folder}.")


# ---------------------------------------------------------------------------
# Convenience wrapper
# ---------------------------------------------------------------------------

def visualize_results(
    full_data: list[np.ndarray],
    families: list[str],
    classes: list[list[Species]],
    chain: np.ndarray,
    AR: np.ndarray,
    R_hat: np.ndarray,
    chain_reshaped: np.ndarray,
    t_burnin: int,
    R_hat_threshold: float,
    n_ranges: list[list[int]],
    num_components: int,
    variable_names: list[str],
    pushforward_results: dict,
    fuel_name: str,
    confidence_width: float = 0.9,
    band_percentiles: bool = False,
    output_dir: str = ".",
) -> None:
    """
    Generate all BayeSAF figures.

    Parameters
    ----------
    pushforward_results : dict
        Output of :func:`bayesaf.postprocessing.props_pushforward.props_pushforward`.
    """
    pr = pushforward_results

    print("Plotting MCMC diagnostics...")
    plot_mcmc_diagnostics(AR, R_hat, R_hat_threshold, t_burnin, num_components, output_dir)

    print("Plotting MCMC traces...")
    plot_traces(chain, num_components, output_dir)

    print("Plotting posterior PDFs...")
    plot_posterior_pdfs(chain_reshaped, num_components, n_ranges, output_dir)

    print("Plotting property pushforward...")
    plot_property_pushforward(
        full_data, variable_names,
        pr["temperature_ranges"], pr["property_MAP"], pr["volumeFraction_MAP"],
        pr["percentiles_list"], pr["mean_model"], pr["percentiles_model"],
        pr["model_output_samples"],
        fuel_name, confidence_width, band_percentiles, output_dir,
    )

    print("Plotting Sobol' indices...")
    plot_sobol(variable_names, pr["sobol_idx"], pr["temperature_ranges"], output_dir)

    print("All figures saved.")
