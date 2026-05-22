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
The of_writer module writes the OpenFOAM liquidProperties class files for
the MAP surrogate components. Thermophysical properties of liquid species are
addressed through an in-house class exploiting Yaws' polynomials (YawsCp,
YawsD, YawsHvap, YawsKappa, YawsMu, YawsPsat, YawsRho, YawsSigma) that
replaces the built-in NSRDS class (except for liquid enthalpy h and second
virial coefficient B). Each surrogate component generates three files:
<species>.C (implementation), <species>.H (class declaration), and
<species>I.H (inline member functions).

Inputs:
1) families      : (1 x numComponents) list of strings denoting the
                   hydrocarbon family of each surrogate component
2) numComponents : number of surrogate mixture components  [-]
3) nc_MAP        : (1 x numComponents) array of numbers of carbon atoms of
                   the MAP surrogate components
4) eta_B_star_MAP: (1 x numComponents) array of normalised topochemical atom
                   indices of the MAP surrogate components
5) fuel_name     : string denoting the real fuel name
------------------------------------------------------------------------
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from bayesaf.thermo_transport.hydrocarbons import _DB_ROOT, _FAMILY_CSV
from bayesaf.thermo_transport.properties import liquid_property

_NAMES_CSV: dict[str, str] = {
    "nparaffins": "nparaffins/nparaffins_Names.csv",
    "isoparaffins": "isoparaffins/isoparaffins_Names.csv",
    "isoparaffins_mono_bis": "isoparaffins_mono_bis/isoparaffins_mono_bis_Names.csv",
    "cycloparaffins": "cycloparaffins/cycloparaffins_Names.csv",
    "dicycloparaffins": "dicycloparaffins/dicycloparaffins_Names.csv",
    "alkylbenzenes": "alkylbenzenes/alkylbenzenes_Names.csv",
    "alkylnaphtalenes": "alkylnaphtalenes/alkylnaphtalenes_Names.csv",
    "cycloaromatics": "cycloaromatics/cycloaromatics_Names.csv",
}

_DIGIT_WORDS = ["zero", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine"]

_OF_HEADER = """\
/*---------------------------------------------------------------------------\\
  =========                |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

"""


def _species_of_name(raw_name: str) -> str:
    """Convert a chemical name to a valid OpenFOAM/C++ identifier."""
    name = raw_name.replace("-", "").replace(",", "")
    for digit, word in enumerate(_DIGIT_WORDS):
        name = name.replace(str(digit), word)
    return name


def _lookup_row(family: str, nC: int, eta_B_star: float) -> tuple[int, pd.DataFrame]:
    """Return (absolute_row_index, full_dataframe) for the best-matching species."""
    df = pd.read_csv(_DB_ROOT / _FAMILY_CSV[family], sep=";")
    nc_rows = df[df["nC"] == nC]
    diff = np.abs(nc_rows["eta_B_star_norm"].values - eta_B_star)
    rel_idx = int(np.argmin(diff))
    abs_idx = nc_rows.index[rel_idx]
    return abs_idx, df


def _write_c_file(path: Path, sp: str, row: pd.Series, Tb: float) -> None:
    with open(path, "w") as f:
        f.write(_OF_HEADER)
        f.write(f'\\*---------------------------------------------------------------------------*/\n\n')
        f.write(f'#include "{sp}.H"\n')
        f.write('#include "addToRunTimeSelectionTable.H"\n\n')
        f.write('// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //\n\n')
        f.write('namespace Foam\n{\n')
        f.write(f'    defineTypeNameAndDebug({sp}, 0);\n')
        f.write(f'    addToRunTimeSelectionTable(liquidProperties, {sp},);\n')
        f.write(f'    addToRunTimeSelectionTable(liquidProperties, {sp}, dictionary);\n')
        f.write('}\n\n')
        f.write('// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //\n\n')
        f.write(f'Foam::{sp}::{sp}()\n:\n')
        f.write('    liquidProperties\n    (\n')
        f.write(f'        {row["W"]},\n')        # molecular weight [g/mol]
        f.write(f'        {row["Tc"]},\n')       # critical temperature [K]
        f.write(f'        {row["Pc"] * 1e5},\n') # critical pressure [Pa]
        f.write(f'        {row["Vc"] * 1e3},\n') # critical volume [m^3/kmol]
        f.write(f'        {row["Zc"]},\n')       # critical compressibility [-]
        f.write(f'        TBD,\n')               # triple point temperature
        f.write(f'        TBD,\n')               # triple point pressure
        f.write(f'        {Tb},\n')              # normal boiling temperature [K]
        f.write(f'        TBD,\n')               # dipole moment
        f.write(f'        {row["omega"]},\n')    # acentric factor
        f.write(f'        TBD\n')               # solubility parameter
        f.write('    ),\n')
        f.write(f'    rho_({row["Arho"]}, {row["Brho"]}, {row["Crho"]}, {row["Tc"]}),\n')
        f.write(f'    pv_({row["Asat"]}, {row["Bsat"]}, {row["Csat"]}, {row["Dsat"]}, {row["Esat"]}),\n')
        f.write(f'    hl_({row["Avap"]}, {row["Bvap"]}, {row["Tc"]}, {row["W"]}),\n')
        f.write(f'    Cp_({row["Ac"]}, {row["Bc"]}, {row["Cc"]}, {row["Dc"]}, {row["W"]}),\n')
        f.write('    h_(TBD1, TBD2, TBD3, TBD4, TBD5, TBD6),\n')
        f.write(f'    Cpg_(TBD1, TBD2, TBD3, TBD4, TBD5, {row["W"]}),\n')
        f.write('    B_(TBD1, TBD2, TBD3, TBD4, TBD5),\n')
        f.write(f'    mu_({row["Amu"]}, {row["Bmu"]}, {row["Cmu"]}, {row["Dmu"]}),\n')
        f.write('    mug_(TBD1, TBD2, TBD3),\n')
        f.write(f'    kappa_({row["Ak"]}, {row["Bk"]}, {row["Ck"]}),\n')
        f.write('    kappag_(TBD1, TBD2, TBD3),\n')
        f.write(f'    sigma_({row["Asigma"]}, {row["Bsigma"]}, {row["Tc"]}),\n')
        f.write('    D_(TBD1, TBD2, TBD3)\n')
        f.write('{}\n\n\n')
        f.write(f'Foam::{sp}::{sp}\n(\n')
        f.write('    const liquidProperties& l,\n')
        f.write('    const thermophysicalFunctions::YawsRho& density,\n')
        f.write('    const thermophysicalFunctions::YawsPsat& vapourPressure,\n')
        f.write('    const thermophysicalFunctions::YawsHvap& heatOfVapourisation,\n')
        f.write('    const thermophysicalFunctions::YawsCp& heatCapacity,\n')
        f.write('    const thermophysicalFunctions::NSRDS0& enthalpy,\n')
        f.write('    const thermophysicalFunctions::YawsCpg& idealGasHeatCapacity,\n')
        f.write('    const thermophysicalFunctions::NSRDS4& secondVirialCoeff,\n')
        f.write('    const thermophysicalFunctions::YawsMu& dynamicViscosity,\n')
        f.write('    const thermophysicalFunctions::YawsMug& vapourDynamicViscosity,\n')
        f.write('    const thermophysicalFunctions::YawsKappa& thermalConductivity,\n')
        f.write('    const thermophysicalFunctions::YawsKappa& vapourThermalConductivity,\n')
        f.write('    const thermophysicalFunctions::YawsSigma& surfaceTension,\n')
        f.write('    const thermophysicalFunctions::YawsD& vapourDiffusivity\n')
        f.write(')\n:\n')
        f.write('    liquidProperties(l),\n')
        f.write('    rho_(density),\n')
        f.write('    pv_(vapourPressure),\n')
        f.write('    hl_(heatOfVapourisation),\n')
        f.write('    Cp_(heatCapacity),\n')
        f.write('    h_(enthalpy),\n')
        f.write('    Cpg_(idealGasHeatCapacity),\n')
        f.write('    B_(secondVirialCoeff),\n')
        f.write('    mu_(dynamicViscosity),\n')
        f.write('    mug_(vapourDynamicViscosity),\n')
        f.write('    kappa_(thermalConductivity),\n')
        f.write('    kappag_(vapourThermalConductivity),\n')
        f.write('    sigma_(surfaceTension),\n')
        f.write('    D_(vapourDiffusivity)\n')
        f.write('{}\n\n\n')
        f.write(f'Foam::{sp}::{sp}(const dictionary& dict)\n:\n')
        f.write(f'    {sp}()\n{{\n')
        f.write('    readIfPresent(*this, dict);\n}\n\n\n')
        f.write('// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //\n\n')
        f.write(f'void Foam::{sp}::write(Ostream& os) const\n{{\n')
        f.write('    liquidProperties::write(*this, os);\n}\n\n\n')
        f.write('// ************************************************************************* //\n')


def _write_h_file(path: Path, sp: str) -> None:
    with open(path, "w") as f:
        f.write(_OF_HEADER)
        f.write(f'Class\n    Foam::{sp}\n\nDescription\n    {sp}\n')
        f.write(f'SourceFiles\n    {sp}.C\n\n\\*---------------------------------------------------------------------------*/\n\n')
        f.write(f'#ifndef {sp}_H\n#define {sp}_H\n\n')
        f.write('#include "liquidProperties.H"\n')
        f.write('#include "YawsRhoThermophysicalFunction.H"\n')
        f.write('#include "YawsCpThermophysicalFunction.H"\n')
        f.write('#include "YawsCpgThermophysicalFunction.H"\n')
        f.write('#include "YawsMuThermophysicalFunction.H"\n')
        f.write('#include "YawsMugThermophysicalFunction.H"\n')
        f.write('#include "YawsKappaThermophysicalFunction.H"\n')
        f.write('#include "YawsSigmaThermophysicalFunction.H"\n')
        f.write('#include "YawsDThermophysicalFunction.H"\n')
        f.write('#include "YawsPsatThermophysicalFunction.H"\n')
        f.write('#include "YawsHvapThermophysicalFunction.H"\n')
        f.write('#include "NSRDS0ThermophysicalFunction.H"\n')
        f.write('#include "NSRDS4ThermophysicalFunction.H"\n\n')
        f.write('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n')
        f.write('namespace Foam\n{\n\n')
        f.write(f'/*---------------------------------------------------------------------------\\\\\n')
        f.write(f'                           Class {sp} Declaration\n')
        f.write(f'\\\\---------------------------------------------------------------------------*/\n\n')
        f.write(f'class {sp}\n:\n    public liquidProperties\n{{\n')
        f.write('    // Private Data\n\n')
        f.write('        thermophysicalFunctions::YawsRho rho_;\n')
        f.write('        thermophysicalFunctions::YawsPsat pv_;\n')
        f.write('        thermophysicalFunctions::YawsHvap hl_;\n')
        f.write('        thermophysicalFunctions::YawsCp Cp_;\n')
        f.write('        thermophysicalFunctions::NSRDS0 h_;\n')
        f.write('        thermophysicalFunctions::YawsCpg Cpg_;\n')
        f.write('        thermophysicalFunctions::NSRDS4 B_;\n')
        f.write('        thermophysicalFunctions::YawsMu mu_;\n')
        f.write('        thermophysicalFunctions::YawsMug mug_;\n')
        f.write('        thermophysicalFunctions::YawsKappa kappa_;\n')
        f.write('        thermophysicalFunctions::YawsKappa kappag_;\n')
        f.write('        thermophysicalFunctions::YawsSigma sigma_;\n')
        f.write('        thermophysicalFunctions::YawsD D_;\n\n\n')
        f.write('public:\n\n')
        f.write('    friend class liquidProperties;\n\n')
        f.write('    //- Runtime type information\n')
        f.write(f'    TypeName("{sp}");\n\n\n')
        f.write('    // Constructors\n\n')
        f.write('        //- Construct null\n')
        f.write(f'        {sp}();\n\n')
        f.write('        //- Construct from components\n')
        f.write(f'        {sp}\n        (\n')
        f.write('            const liquidProperties& l,\n')
        f.write('            const thermophysicalFunctions::YawsRho& density,\n')
        f.write('            const thermophysicalFunctions::YawsPsat& vapourPressure,\n')
        f.write('            const thermophysicalFunctions::YawsHvap& heatOfVapourisation,\n')
        f.write('            const thermophysicalFunctions::YawsCp& heatCapacity,\n')
        f.write('            const thermophysicalFunctions::NSRDS0& enthalpy,\n')
        f.write('            const thermophysicalFunctions::YawsCpg& idealGasHeatCapacity,\n')
        f.write('            const thermophysicalFunctions::NSRDS4& secondVirialCoeff,\n')
        f.write('            const thermophysicalFunctions::YawsMu& dynamicViscosity,\n')
        f.write('            const thermophysicalFunctions::YawsMug& vapourDynamicViscosity,\n')
        f.write('            const thermophysicalFunctions::YawsKappa& thermalConductivity,\n')
        f.write('            const thermophysicalFunctions::YawsKappa& vapourThermalConductivity,\n')
        f.write('            const thermophysicalFunctions::YawsSigma& surfaceTension,\n')
        f.write('            const thermophysicalFunctions::YawsD& vapourDiffusivity\n')
        f.write('        );\n\n')
        f.write('        //- Construct from dictionary\n')
        f.write(f'        {sp}(const dictionary& dict);\n\n')
        f.write('        //- Construct and return clone\n')
        f.write('        virtual autoPtr<liquidProperties> clone() const\n        {\n')
        f.write(f'            return autoPtr<liquidProperties>(new {sp}(*this));\n        }}\n\n\n')
        f.write('    // Member Functions\n\n')
        for fn, comment in [
            ("rho",    "Liquid density [kg/m^3]"),
            ("pv",     "Vapour pressure [Pa]"),
            ("hl",     "Heat of vapourisation [J/kg]"),
            ("Cp",     "Liquid heat capacity [J/kg/K]"),
            ("h",      "Liquid enthalpy [J/kg]"),
            ("Cpg",    "Ideal gas heat capacity [J/kg/K]"),
            ("B",      "Second Virial Coefficient [m^3/kg]"),
            ("mu",     "Liquid viscosity [Pa s]"),
            ("mug",    "Vapour viscosity [Pa s]"),
            ("kappa",  "Liquid thermal conductivity [W/m/K]"),
            ("kappag", "Vapour thermal conductivity [W/m/K]"),
            ("sigma",  "Surface tension [N/m]"),
            ("D",      "Vapour diffusivity [m^2/s]"),
        ]:
            f.write(f'        //- {comment}\n')
            f.write(f'        inline scalar {fn}(scalar p, scalar T) const;\n\n')
        f.write('        //- Vapour diffusivity [m^2/s] with specified binary pair\n')
        f.write('        inline scalar D(scalar p, scalar T, scalar Wb) const;\n\n\n')
        f.write('    // I-O\n\n')
        f.write('        //- Write the function coefficients\n')
        f.write('        void write(Ostream& os) const;\n};\n\n\n')
        f.write('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n')
        f.write('} // End namespace Foam\n\n')
        f.write('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n')
        f.write(f'#include "{sp}I.H"\n\n')
        f.write('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n')
        f.write('#endif\n\n')
        f.write('// ************************************************************************* //\n')


def _write_ih_file(path: Path, sp: str) -> None:
    with open(path, "w") as f:
        f.write(_OF_HEADER)
        f.write('\\*---------------------------------------------------------------------------*/\n\n')
        for fn in ["rho", "pv", "hl", "Cp", "h", "Cpg", "B", "mu", "mug",
                   "kappa", "kappag", "sigma", "D"]:
            f.write(f'inline Foam::scalar Foam::{sp}::{fn}(scalar p, scalar T) const\n{{\n')
            f.write(f'    return {fn}_.f(p, T);\n}}\n\n\n')
        f.write(f'inline Foam::scalar Foam::{sp}::D(scalar p, scalar T, scalar Wb) const\n{{\n')
        f.write('    return D_.f(p, T, Wb);\n}\n\n\n')
        f.write('// ************************************************************************* //\n')


def write_openfoam(
    families: list[str],
    num_components: int,
    nc_MAP: np.ndarray,
    eta_B_star_MAP: np.ndarray,
    fuel_name: str,
    output_dir: str = ".",
) -> None:
    """
    Write OpenFOAM liquidProperties class files for each MAP surrogate component.

    Parameters
    ----------
    families : list[str]
    num_components : int
    nc_MAP : ndarray, shape (Nc,)
    eta_B_star_MAP : ndarray, shape (Nc,)
    fuel_name : str
    output_dir : str
        Root directory for output (default: current working directory).
    """
    fuel_id = fuel_name.replace(" ", "").replace("-", "")
    base = Path(output_dir) / "OpenFOAM" / fuel_id

    for l in range(num_components):
        abs_idx, df = _lookup_row(families[l], int(nc_MAP[l]), float(eta_B_star_MAP[l]))
        row = df.loc[abs_idx]

        # Look up species name
        names_df = pd.read_csv(_DB_ROOT / _NAMES_CSV[families[l]], sep=None, engine="python")
        raw_name = str(names_df.iloc[int(abs_idx)]["Name"])
        sp = _species_of_name(raw_name)

        # Compute normal boiling temperature
        from bayesaf.thermo_transport.hydrocarbons import Species
        sp_obj = Species(
            nC=int(row["nC"]),
            eta_B_star=float(row.get("eta_B_star", 0.0)),
            eta_B_star_norm=float(row["eta_B_star_norm"]),
            mol_weight=float(row["W"]),
            Tc=float(row["Tc"]),
            coeff_mu=np.array([row["Amu"], row["Bmu"], row["Cmu"], row["Dmu"]]),
            coeff_rho=np.array([row["Arho"], row["Brho"], row["Crho"]]),
            coeff_psat=np.array([row["Asat"], row["Bsat"], row["Csat"], row["Dsat"], row["Esat"]]),
            coeff_k=np.array([row["Ak"], row["Bk"], row["Ck"]]),
            coeff_cl=np.array([row["Ac"], row["Bc"], row["Cc"], row["Dc"]]),
            coeff_hv=np.array([row["Avap"], row["Bvap"]]),
            coeff_sigma=np.array([row["Asigma"], row["Bsigma"]]),
            DCN=float(row.get("DCN", 0.0)),
            Tf=float(row.get("Tf", 0.0)),
            Tfz=float(row.get("Tfz", 0.0)),
            Hc=float(row.get("Hc", 0.0)),
        )
        Tb = float(liquid_property("boilingTemperature", 0.0, sp_obj, 101325.0))

        folder = base / sp
        folder.mkdir(parents=True, exist_ok=True)

        _write_c_file(folder / f"{sp}.C", sp, row, Tb)
        _write_h_file(folder / f"{sp}.H", sp)
        _write_ih_file(folder / f"{sp}I.H", sp)

        print(f"  [{l+1}/{num_components}] OpenFOAM files written for {sp!r} → {folder}")

    print(f"OpenFOAM files written under '{base}'.")
