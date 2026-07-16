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
The map_writer module writes a text file containing information about the
maximum a posteriori (MAP) surrogate in the same format as the MATLAB
MAPwriter. For each surrogate component it identifies the chemical species
whose number of carbon atoms and normalised topochemical atom index most
closely match the MAP parameters, computes blend properties using the same
mixing rules as MATLAB, and writes the results to a summary file.

Inputs:
1) families      : (1 x numComponents) list of strings denoting the
                   hydrocarbon family of each surrogate component
2) numComponents : number of surrogate mixture components  [-]
3) x_MAP         : (1 x numComponents) array of molar fractions of the MAP
                   surrogate components
4) nc_MAP        : (1 x numComponents) array of numbers of carbon atoms of
                   the MAP surrogate components
5) eta_B_star_MAP: (1 x numComponents) array of normalised topochemical atom
                   indices of the MAP surrogate components
6) fuel_name     : string denoting the real fuel name

Outputs:
1) species_MAP: (1 x numComponents) list of species name strings for the
                MAP surrogate
------------------------------------------------------------------------
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd

from bayesaf.thermo_transport.hydrocarbons import _DB_ROOT, _FAMILY_CSV
from bayesaf.utilities.composition import mol_to_mass, mol_weight

_NAMES_CSV: dict[str, str] = {
    "nparaffins": "nparaffins/nparaffins_Names.csv",
    "isoparaffins": "isoparaffins/isoparaffins_Names.csv",
    "isoparaffins_mono_bis": "isoparaffins_mono_bis/isoparaffins_mono_bis_Names.csv",
    "cycloparaffins": "cycloparaffins/cycloparaffins_Names.csv",
    "dicycloparaffins": "dicycloparaffins/dicycloparaffins_Names.csv",
    "alkylbenzenes": "alkylbenzenes/alkylbenzenes_Names.csv",
    "alkylnaphthalenes": "alkylnaphthalenes/alkylnaphthalenes_Names.csv",
    "cycloaromatics": "cycloaromatics/cycloaromatics_Names.csv",
}

# H-atom count formula per family
_H_FORMULA: dict[str, str] = {
    "nparaffins":            "2n+2",
    "isoparaffins":          "2n+2",
    "isoparaffins_mono_bis": "2n+2",
    "cycloparaffins":        "2n",
    "dicycloparaffins":      "2n-2",
    "alkylbenzenes":         "2n-6",
    "alkylnaphthalenes":     "2n-12",
    "cycloaromatics":        "2n-8",
}

_R_UNIV = 8.314   # J/(mol·K)
_T_RHO  = 300.0   # K — reference temperature for density / volume fractions (matches MATLAB)


def _n_hydrogen(family: str, nC: int) -> int:
    formula = _H_FORMULA[family]
    if formula == "2n+2":
        return 2 * nC + 2
    elif formula == "2n":
        return 2 * nC
    elif formula == "2n-2":
        return 2 * nC - 2
    elif formula == "2n-6":
        return 2 * nC - 6
    elif formula == "2n-12":
        return 2 * nC - 12
    elif formula == "2n-8":
        return 2 * nC - 8
    return 0


def _lookup_species_full(family: str, nC: int, eta_B_star: float) -> dict:
    """
    Return a dict of all MAP-relevant properties for the best-matching species.

    The species is identified as the row in the database CSV whose nC matches
    and whose eta_B_star_norm is closest to *eta_B_star*.
    """
    df       = pd.read_csv(_DB_ROOT / _FAMILY_CSV[family], sep=";")
    nc_rows  = df[df["nC"] == nC]
    if nc_rows.empty:
        return {"name": "Unknown"}

    diff = np.abs(nc_rows["eta_B_star_norm"].values - eta_B_star)
    best     = nc_rows.iloc[np.argmin(diff)]
    row_idx  = best.name          # integer label in df

    names_df = pd.read_csv(_DB_ROOT / _NAMES_CSV[family], sep=None, engine="python")
    name     = str(names_df.iloc[row_idx]["Name"])

    return {
        "name":   name,
        "W":      float(best["W"]),      # g/mol
        "DCN":    float(best["DCN"]),
        "Tf":     float(best["Tf"]),     # flash point [K]
        "Tfz":    float(best["Tfz"]),    # freezing point [K]
        "Hc":     float(best["Hc"]),     # LHV [MJ/kg]
        "Tc":     float(best["Tc"]),     # critical temperature [K]
        "Pc":     float(best["Pc"]),     # critical pressure [bar]
        "Vc":     float(best["Vc"]),     # critical volume [m³/mol]
        "omega":  float(best["omega"]),  # acentric factor
        "Arho":   float(best["Arho"]),   # Yaws density coefficients
        "Brho":   float(best["Brho"]),
        "Crho":   float(best["Crho"]),
    }


def _yaws_rho(props: dict, T: float) -> float:
    """Yaws liquid density [kg/m³] at temperature T [K]."""
    return 1000.0 * props["Arho"] * props["Brho"] ** (-(1.0 - T / props["Tc"]) ** props["Crho"])


def _vec_fmt(values: np.ndarray, decimals: int = 5, sep: str = "     ") -> str:
    """
    Format a numeric array as space-separated values, mimicking MATLAB num2str.

    Integer-valued floats (e.g. 1.0) are printed without a decimal point;
    all others use *decimals* decimal places.
    """
    parts = []
    for v in values:
        if abs(v - round(v)) < 1e-9:
            parts.append(str(int(round(v))))
        else:
            parts.append(f"{v:.{decimals}f}")
    return sep.join(parts)


def write_map(
    families: list[str],
    num_components: int,
    x_MAP: np.ndarray,
    nc_MAP: np.ndarray,
    eta_B_star_MAP: np.ndarray,
    fuel_name: str,
    output_file: str = "MAP.txt",
) -> list[str]:
    """
    Write MAP surrogate information to *output_file* (MATLAB-compatible format)
    and return the list of species names.

    Parameters
    ----------
    families : list[str]
    num_components : int
    x_MAP : ndarray, shape (Nc,)
    nc_MAP : ndarray, shape (Nc,)
    eta_B_star_MAP : ndarray, shape (Nc,)
    fuel_name : str
    output_file : str

    Returns
    -------
    list[str]
        Species names of the MAP components.
    """
    Nc = num_components

    # ------------------------------------------------------------------
    # 1. Look up species and collect all properties
    # ------------------------------------------------------------------
    all_props: list[dict] = []
    for k in range(Nc):
        p = _lookup_species_full(families[k], int(nc_MAP[k]), float(eta_B_star_MAP[k]))
        all_props.append(p)

    species_names = [p["name"] for p in all_props]
    W_i    = np.array([p["W"]    for p in all_props])   # g/mol
    DCN_i  = np.array([p["DCN"]  for p in all_props])
    Tf_i   = np.array([p["Tf"]   for p in all_props])   # K
    Tfz_i  = np.array([p["Tfz"]  for p in all_props])   # K
    Hc_i   = np.array([p["Hc"]   for p in all_props])   # MJ/kg
    Tc_i   = np.array([p["Tc"]   for p in all_props])   # K
    Vc_i   = np.array([p["Vc"]   for p in all_props])   # m³/mol
    omega_i = np.array([p["omega"] for p in all_props])

    # ------------------------------------------------------------------
    # 2. Mixture molecular weight and composition
    # ------------------------------------------------------------------
    W_mean = mol_weight(x_MAP, W_i)           # g/mol
    y_MAP  = mol_to_mass(x_MAP, W_i)          # mass fractions

    # Match MATLAB: round(sum(x.*nc), 2) and round(nH, 2) then round(H/C, 3)
    nC_mix = round(float(np.dot(x_MAP, nc_MAP)), 2)
    nH_mix = round(
        sum(float(x_MAP[k]) * _n_hydrogen(families[k], int(nc_MAP[k]))
            for k in range(Nc)),
        2,
    )
    HC_ratio = round(nH_mix / nC_mix, 3) if nC_mix > 0 else 0.0

    # ------------------------------------------------------------------
    # 3. Density and volume fractions  (MATLAB: Vi_MAP = rho_mix*y_MAP./rho_i)
    # ------------------------------------------------------------------
    rho_i   = np.array([_yaws_rho(all_props[k], _T_RHO) for k in range(Nc)])  # kg/m³
    rho_mix = np.dot(x_MAP, W_i) / np.sum(x_MAP * W_i / rho_i)               # kg/m³
    V_MAP   = rho_mix * y_MAP / rho_i                                          # volume fractions

    # ------------------------------------------------------------------
    # 4. DCN  (volume-weighted)
    # ------------------------------------------------------------------
    DCN_blend = float(np.dot(V_MAP, DCN_i))

    # ------------------------------------------------------------------
    # 5. Flash point  (mass-weighted blending index, ASTM D-7215 approach)
    #    MATLAB: BI = 51708 * exp((log(Tf_F) - 2.6287)^2 / -0.91725)
    # ------------------------------------------------------------------
    Tf_F     = (Tf_i - 273.15) * 9.0 / 5.0 + 32.0          # K → °F
    BI_flash = 51708.0 * np.exp((np.log(Tf_F) - 2.6287)**2 / (-0.91725))
    BI_blend_flash   = float(np.dot(y_MAP, BI_flash))
    flash_point_F    = math.exp(math.sqrt(-0.91725 * math.log(BI_blend_flash / 51708.0)) + 2.6287)
    flash_point      = (flash_point_F - 32.0) * 5.0 / 9.0 + 273.15  # °F → K

    # ------------------------------------------------------------------
    # 6. Freezing point  (volume-weighted BI^(1/0.05))
    #    MATLAB: BI_i = Tfz_i^(1/0.05);  freezing = (sum(Vi * BI_i))^0.05
    # ------------------------------------------------------------------
    BI_freeze    = Tfz_i ** (1.0 / 0.05)
    freezing_point = float(np.dot(V_MAP, BI_freeze)) ** 0.05

    # ------------------------------------------------------------------
    # 7. LHV  (mass-weighted)
    # ------------------------------------------------------------------
    LHV = float(np.dot(y_MAP, Hc_i))

    # ------------------------------------------------------------------
    # 8. Lee-Kesler pseudocritical mixing rules
    #    Vc_ij = (1/8)*(Vc_i^(1/3) + Vc_j^(1/3))^3
    #    Tc_ij = sqrt(Tc_i * Tc_j)
    #    Vc_m  = sum_ij xi*xj*Vc_ij
    #    Tc_m  = sum_ij xi*xj*Vc_ij^(1/4)*Tc_ij  /  Vc_m^(1/4)
    #    omega_m = sum_i xi*omega_i
    #    Pc_m  = (0.2905 - 0.085*omega_m) * R * Tc_m / Vc_m / 1e5   [bar]
    # ------------------------------------------------------------------
    Vc_m     = 0.0
    Tc_m_num = 0.0
    for ii in range(Nc):
        for jj in range(Nc):
            Vc_ij     = (1.0 / 8.0) * (Vc_i[ii]**(1.0/3.0) + Vc_i[jj]**(1.0/3.0))**3
            Tc_ij     = math.sqrt(Tc_i[ii] * Tc_i[jj])
            Vc_m     += x_MAP[ii] * x_MAP[jj] * Vc_ij
            Tc_m_num += x_MAP[ii] * x_MAP[jj] * Vc_ij**(1.0/4.0) * Tc_ij
    Tc_m    = Tc_m_num / Vc_m**(1.0/4.0)
    omega_m = float(np.dot(x_MAP, omega_i))
    Pc_m    = (0.2905 - 0.085 * omega_m) * _R_UNIV * Tc_m / Vc_m / 1e5   # bar

    # ------------------------------------------------------------------
    # 9. Write output file  (MATLAB MAPwriter format)
    # ------------------------------------------------------------------
    classes_str = ", ".join(families)
    hydro_str   = ", ".join(species_names)
    nc_str      = "  ".join(str(int(nc)) for nc in nc_MAP)
    eta_str     = _vec_fmt(eta_B_star_MAP)
    x_str       = _vec_fmt(x_MAP)
    y_str       = _vec_fmt(y_MAP)
    V_str       = _vec_fmt(V_MAP)

    header = (
        "    ____                  _____ ___    ______                               \n"
        "   / __ )____ ___  _____ / ___//   |  / ____/                               \n"
        "  / __  / __ `/ / / / _ \\__ \\/ /| | / /_                                \n"
        " / /_/ / /_/ / /_/ /  __/__/ / ___ |/ __/                                   \n"
        "/_____/\\__,_/\\__, /\\___/____/_/  |_/_/                                   \n"
        "            /____/                                                          \n"
        "\n"
        "--------------------------------------------------------------------------  \n"
        "Contributors / Copyright Notice                                             \n"
        "© 2026 Jacopo Liberatori — jacopo.liberatori@centralesupelec.fr             \n"
        "Postdoctoral Researcher @ Laboratoire EM2C, CentraleSupélec (CNRS)          \n"
        "\n"
        "© 2026 Davide Cavalieri — davide.cavalieri@uniroma1.it                      \n"
        "Postdoctoral Researcher @ Sapienza University of Rome,                      \n"
        "Department of Mechanical and Aerospace Engineering (DIMA)                   \n"
        "\n"
        "© 2026 Matteo Blandino, Ph.D.                                               \n"
        "\n"
        "Reference:                                                                  \n"
        "J. Liberatori, D. Cavalieri, M. Blandino, M. Valorani, and P.P. Ciottoli.   \n"
        "BayeSAF: Emulation and Design of Sustainable Alternative Fuels via Bayesian \n"
        "Inference and Descriptors-Based Machine Learning. Fuel 419, 138835 (2026).  \n"
        "Available at: https://doi.org/10.1016/j.fuel.2026.138835.                    \n"
        "-----------------------------------------------------------------------     \n"
    )

    with open(output_file, "w") as fid:
        fid.write(header)
        fid.write(f"\n+++++ MAP surrogate for {fuel_name} +++++\n")
        fid.write(f"Number of MAP components: {Nc}\n")
        fid.write(f"MAP hydrocarbon classes: {classes_str}\n")
        fid.write(f"MAP hydrocarbons: {hydro_str}\n")
        fid.write(f"MAP molecular weight: {W_mean:.7g} g/mol\n")
        fid.write(f"MAP molecular formula: C{nC_mix:.7g} H{nH_mix:.7g}\n")
        fid.write(f"MAP hydrogen-to-carbon ratio: {HC_ratio:.7g}\n")
        fid.write(f"MAP derived cetane number: {DCN_blend:.7g}\n")
        fid.write(f"MAP flash point: {round(flash_point, 2):.7g} K\n")
        fid.write(f"MAP freezing point: {round(freezing_point, 2):.7g} K\n")
        fid.write(f"MAP lower heating value: {round(LHV, 2):.7g} MJ/kg\n")
        fid.write(f"MAP molar fractions: {x_str}\n")
        fid.write(f"MAP mass fractions: {y_str}\n")
        fid.write(f"MAP volume fractions: {V_str}\n")
        fid.write(f"MAP carbon atoms: {nc_str}\n")
        fid.write(f"MAP topochemical atom indices: {eta_str}\n")
        fid.write(
            f"\nMAP pseudocritical temperature according to Lee-Kesler mixing rule: "
            f"{Tc_m:.7g} K\n"
        )
        fid.write(
            f"MAP pseudocritical pressure according to Lee-Kesler mixing rule: "
            f"{Pc_m:.7g} bar\n"
        )

    print(f"MAP summary written to '{output_file}'.")
    return species_names
