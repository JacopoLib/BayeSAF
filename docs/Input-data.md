# Input data

> **Names.** This page uses the Python names. In MATLAB: `DATASET` is `Dataset` (a cell array of names without extension), `DATA_DIR` does not exist (the files are found on the
> MATLAB path), `FAMILIES` is `families`, `N_RANGES` is `n_ranges`, `T_BUBBLE` is `Tbubble_realFuel`. The file formats are the same. The
> notes specific to MATLAB are at the end of the page.

## Species database (`database/`)

One directory per hydrocarbon family. The directory name and the family key used in the code are:

| Family key in `FAMILIES` | Directory | Notes |
|---|---|---|
| `nparaffins` | `nparaffins/` | |
| `isoparaffins` | `isoparaffins/` | |
| `cycloparaffins` | `cycloparaffins/` | |
| `dicycloparaffins` | `dicycloparaffins/` | |
| `alkylbenzenes` | `alkylbenzenes/` | |
| `alkylnaphthalenes` | `alkylnaphthalenes/` | |
| `cycloaromatics` | `cycloaromatics/` | |

Each directory contains three semicolon-separated CSV files:

| File | Content |
|---|---|
| `<family>_MatlabData.csv` | the data read by `load_hydrocarbons` |
| `<family>_Names.csv` | `Name; eta_B_star; eta_B_star_norm`, row-aligned with the data file |
| `<family>_SMILES.csv` | SMILES of each species |

### Columns of `<family>_MatlabData.csv`

| Column | Meaning |
|---|---|
| `nC` | number of carbon atoms |
| `eta_B_star` | topochemical atom index for branching |
| `eta_B_star_norm` | normalised index |
| `W` | molecular weight [g/mol] |
| `Tc`, `Pc`, `Vc`, `rhoc`, `Zc`, `omega` | critical temperature [K], pressure, molar volume, density, compressibility factor, acentric factor |
| `Tf` | flash point [K] |
| `Tfz` | freezing point [K] |
| `Hc` | lower heating value [MJ/kg] |
| `Amu, Bmu, Cmu, Dmu` | Yaws coefficients, dynamic viscosity |
| `Arho, Brho, Crho` | liquid density |
| `Asat, Bsat, Csat, Dsat, Esat` | vapour pressure |
| `Ac, Bc, Cc, Dc` | liquid specific heat |
| `Ak, Bk, Ck` | thermal conductivity |
| `Avap, Bvap` | latent heat of vaporization |
| `Asigma, Bsigma` | surface tension |
| `DCN` | derived cetane number |

`load_hydrocarbons` reads `nC`, `eta_B_star`, `eta_B_star_norm`, `W`, `Tc`, the Yaws columns, `DCN`,
`Tf`, `Tfz` and `Hc`. The columns `Pc`, `Vc`, `rhoc`, `Zc` and `omega` are used only by the MAP
writer for the pseudo-critical properties. The correspondence of the Yaws coefficients with the
formulas is given in [Thermophysical models](Thermophysical-models.md).

The row order of `_Names.csv` must match `_MatlabData.csv`, because the MAP writer
looks up the species name by row position.

## Experimental targets (`exp/<FUEL>/`)

One CSV file per property, with a header line and comma-separated numeric columns.
`data_import` reads `data_dir/<name>.csv` for each entry of `DATASET`.

### Mapping from file name to property

The property is detected from the **file name by substring**, in this order (first match wins):

`deltaT_dist`, `distillation`, `molWeight`, `specificHeat`, `latentHeat`, `vaporPressure`, `kappa`,
`sigma`, `rho`, `mu`, `nu`, `HC`, `DCN`, `flash`, `freezing`, `LHV`.

So `rho_expA` is read as `rho`. A name that does not contain any of these strings raises `ValueError`.
A name that contains two keys takes the first one in this list (for example `mu_nu` is `mu`). Keep the names exact.

### Column layout

The columns are positional, not read by name.

| Property | Columns |
|---|---|
| `molWeight`, `HC`, `DCN`, `LHV`, `flash`, `freezing` | `value, std` |
| `rho`, `mu`, `nu`, `kappa`, `specificHeat`, `latentHeat`, `vaporPressure`, `sigma` | `T [K], value, std` |
| `distillation` | `recovered volume [%], T [K], std [K], pressure [Pa]` |
| `deltaT_dist` | `range code, deltaT [K], std [K], pressure [Pa]` |

Units of the value column must match the model: `rho` in kg/m3, `mu` in Pa s, `nu` in m2/s,
`kappa` in W/(m K), `specificHeat` in J/(kg K), `latentHeat` in J/kg, `vaporPressure` in Pa,
`sigma` in N/m, `LHV` in MJ/kg, `molWeight` in g/mol.

The `std` of scalar properties (`molWeight`, `HC`, `DCN`, `flash`, `freezing`, `LHV`) is column 1; for the other
properties it is column 2 (`pdfs/pdf_build.py`, `_std_data`).

`deltaT_dist` range codes: `5010` is `T50 - T10`, `9010` is `T90 - T10`. Rows with other codes
produce a model value of zero.

`distillation` volumes should be in the range 0-100 and ordered by increasing volume. The cheap
likelihood uses the **first row** of the file, so put there the data point that should act as
the filter (in the CHRJ example it is the 5% point).

### Distillation pressure

`p_distillation` reads the pressure from column 4 of the first row of the `distillation` file (or, if
absent, of `deltaT_dist`). If neither is in `DATASET` it uses 101325 Pa.

### Temperature filtering

For the temperature-dependent properties (`rho`, `mu`, `nu`, `kappa`, `specificHeat`, `latentHeat`,
`vaporPressure`, `sigma`), `data_import` discards the rows with `T >= 0.999 * T_threshold`, where

```
T_threshold = min(T_BUBBLE, lowest normal boiling temperature among all candidate species)
```

computed at the distillation pressure. If you do not know the bubble temperature of the real fuel,
set `T_BUBBLE = 1e18`: the threshold is then the lowest boiling temperature among the candidates.

Scalar properties and the two distillation datasets are not filtered.

### Returned structures

`data_import(...)` returns four lists with one element per entry of `DATASET`:

| Return | Content |
|---|---|
| `full_data` | numeric array of the (filtered) file |
| `min_T_array`, `max_T_array` | minimum and maximum of column 0, or of column 1 (the temperature) for `distillation` and `deltaT_dist` |
| `variable_names` | detected property names |

## Example (`exp/CHRJ_POSF6152/`)

```
rho.csv            T [K], rho [kg/m^3], std [kg/m^3]
distillation.csv   beta [-], T [K], std [K], pressure [Pa]
deltaT_dist.csv    range, deltaT [K], std [K], pressure [Pa]
molWeight.csv      MW [g/mol], std [g/mol]
HC.csv             HC [-], std [-]
DCN.csv            DCN [-], std [-]
LHV.csv            LHV [MJ/kg], std [MJ/kg]
flash.csv          flash [K], std [K]
```

## Checklist for a new fuel

1. One CSV per property, in the units above, with a `std` for every row.
2. At least one `distillation` or `deltaT_dist` file if you need a pressure other than 101325 Pa.
3. `DATASET` lists only files that exist in `DATA_DIR`.
4. A sensible `T_BUBBLE`.
5. The family and carbon range of each component cover the chemistry of the fuel. A range too narrow cannot be recovered by the inference.
6. Informative bounds on the molar fractions, if the composition is roughly known (for example from GCxGC), speed up convergence.

## MATLAB implementation: notes

- **Where the files are found.** MATLAB looks for the CSV files on the path. The examples add the folders of the database
  and `exp/<FUEL>/` with `addpath` (relative paths, so they must be run from the folder of the example). `Dataset` contains
  the file names **without extension** (`'rho'`, `'distillation'`, ...); `readtable` finds `rho.csv`.
- **Detection of the property.** `data_import` tests the name of each entry with `contains`, in the same order as the Python implementation
  (the list above): the first match decides, and a name that matches none stops with an error.
- **Minimum and maximum.** As in Python: `minT_array` and `maxT_array` come from column 1 (temperature) of `distillation` and `deltaT_dist`
  data and from column 0 for the other properties.
- **Pressure.** `p_distillation(Dataset)` reads column 4 of the first row of the file whose name contains `distillation`, or else
  `deltaT_dist`, and uses 101325 Pa if there is none.
- **Temperature threshold.** Same rule as in Python (`threshold_temperature.m`); the boiling temperatures are computed with `fsolve`.
- **Species database.** `Hydrocarbons.m` reads `<family>_MatlabData.csv` with `readtable`, so the folder of each family must be on the path.
