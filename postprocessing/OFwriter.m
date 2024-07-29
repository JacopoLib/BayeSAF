function OFwriter(families, numComponents, nc_MAP, eta_B_star_MAP, fuel_name)

% DESCRIPTION:
% The OFwriter function writes the OpenFOAM liquidProperties class files for the MAP surrogate components.
% Note that thermophysical properties of the liquid species are addressed through an in-house class exploiting
% Yaws' polynomials (YawsCp, YawsD, YawsHvap, YawsKappa, YawsMu, YawsPsat, YawsRho, YawsSigma) that replaces
% the built-in NSRDS (National Standard Reference Data Series) class leveraging the NSRDS thermophysical functions
% (except for liquid enthalpy, h, and second virial coefficient, B).

% Inputs:
% 1) families      : (1 x numComponents) cell array, with the i-th cell containing a string
% that denotes the i-th hydrocarbon family under consideration
% 2) numComponents : number of surrogate mixture components  [-]
% 3) nc_MAP: (1 x numComponents) array containing the numbers of carbon atoms characterizing the MAP surrogate components
% 4) eta_B_star_MAP: (1 x numComponents) array containing the topochemical atom indices characterizing the MAP surrogate components
% 5) fuel_name     : string denoting the real fuel name

% ------------------------------------------------------------------------%
% Contributors/Copyright
% 2024 Jacopo Liberatori, jacopo.liberatori@uniroma1.it
% 2024 Davide Cavalieri,  davide.cavalieri@uniroma1.it
% Department of Mechanical and Aerospace Engineering (DIMA)
% Sapienza University of Rome
% ------------------------------------------------------------------------%

for l = 1:numComponents
    
    if strcmp(families(l), 'nparaffins')
        data = 'nparaffins/nparaffins_MatlabData.csv';
        table = readtable(data, 'VariableNamingRule', 'preserve');
        nc_idx = find(table.nC == nc_MAP(l));
        eta_B_star_diff = abs(table.eta_B_star_norm(nc_idx) - eta_B_star_MAP(l));
        min_eta_B_star_diff = min(abs(table.eta_B_star_norm(nc_idx) - eta_B_star_MAP(l)));
        row_MAP_idx = find(eta_B_star_diff == min_eta_B_star_diff, 1);

        names = 'nparaffins/nparaffins_Names.csv';
        table_names = readtable(names, 'VariableNamingRule', 'preserve');

        species_MAP = table_names.Name(nc_idx(1) + row_MAP_idx - 1);
    elseif strcmp(families(l), 'isoparaffins')
        data = 'isoparaffins/isoparaffins_MatlabData.csv';
        table = readtable(data, 'VariableNamingRule', 'preserve');
        nc_idx = find(table.nC == nc_MAP(l));
        eta_B_star_diff = abs(table.eta_B_star_norm(nc_idx) - eta_B_star_MAP(l));
        min_eta_B_star_diff = min(abs(table.eta_B_star_norm(nc_idx) - eta_B_star_MAP(l)));
        row_MAP_idx = find(eta_B_star_diff == min_eta_B_star_diff, 1);

        names = 'isoparaffins/isoparaffins_Names.csv';
        table_names = readtable(names, 'VariableNamingRule', 'preserve');

        species_MAP = table_names.Name(nc_idx(1) + row_MAP_idx - 1);
    elseif strcmp(families(l), 'cycloparaffins')
        data = 'cycloparaffins/cycloparaffins_MatlabData.csv';
        table = readtable(data, 'VariableNamingRule', 'preserve');
        nc_idx = find(table.nC == nc_MAP(l));
        eta_B_star_diff = abs(table.eta_B_star_norm(nc_idx) - eta_B_star_MAP(l));
        min_eta_B_star_diff = min(abs(table.eta_B_star_norm(nc_idx) - eta_B_star_MAP(l)));
        row_MAP_idx = find(eta_B_star_diff == min_eta_B_star_diff, 1);

        names = 'cycloparaffins/cycloparaffins_Names.csv';
        table_names = readtable(names, 'VariableNamingRule', 'preserve');

        species_MAP = table_names.Name(nc_idx(1) + row_MAP_idx - 1);
    elseif strcmp(families(l), 'dicycloparaffins')
        data = 'dicycloparaffins/dicycloparaffins_MatlabData.csv';
        table = readtable(data, 'VariableNamingRule', 'preserve');
        nc_idx = find(table.nC == nc_MAP(l));
        eta_B_star_diff = abs(table.eta_B_star_norm(nc_idx) - eta_B_star_MAP(l));
        min_eta_B_star_diff = min(abs(table.eta_B_star_norm(nc_idx) - eta_B_star_MAP(l)));
        row_MAP_idx = find(eta_B_star_diff == min_eta_B_star_diff, 1);

        names = 'dicycloparaffins/dicycloparaffins_Names.csv';
        table_names = readtable(names, 'VariableNamingRule', 'preserve');

        species_MAP = table_names.Name(nc_idx(1) + row_MAP_idx - 1);
    elseif strcmp(families(l), 'alkylbenzenes')
        data = 'alkylbenzenes/alkylbenzenes_MatlabData.csv';
        table = readtable(data, 'VariableNamingRule', 'preserve');
        nc_idx = find(table.nC == nc_MAP(l));
        eta_B_star_diff = abs(table.eta_B_star_norm(nc_idx) - eta_B_star_MAP(l));
        min_eta_B_star_diff = min(abs(table.eta_B_star_norm(nc_idx) - eta_B_star_MAP(l)));
        row_MAP_idx = find(eta_B_star_diff == min_eta_B_star_diff, 1);

        names = 'alkylbenzenes/alkylbenzenes_Names.csv';
        table_names = readtable(names, 'VariableNamingRule', 'preserve');

        species_MAP = table_names.Name(nc_idx(1) + row_MAP_idx - 1);
    elseif strcmp(families(l), 'alkylnaphtalenes')
        data = 'alkylnaphtalenes/alkylnaphtalenes_MatlabData.csv';
        table = readtable(data, 'VariableNamingRule', 'preserve');
        nc_idx = find(table.nC == nc_MAP(l));
        eta_B_star_diff = abs(table.eta_B_star_norm(nc_idx) - eta_B_star_MAP(l));
        min_eta_B_star_diff = min(abs(table.eta_B_star_norm(nc_idx) - eta_B_star_MAP(l)));
        row_MAP_idx = find(eta_B_star_diff == min_eta_B_star_diff, 1);

        names = 'alkylnaphtalenes/alkylnaphtalenes_Names.csv';
        table_names = readtable(names, 'VariableNamingRule', 'preserve');

        species_MAP = table.Name(nc_idx(1) + row_MAP_idx - 1);
    elseif strcmp(families(l), 'cycloaromatics')
        data = 'cycloaromatics/cycloaromatics_MatlabData.csv';
        table = readtable(data, 'VariableNamingRule', 'preserve');
        nc_idx = find(table.nC == nc_MAP(l));
        eta_B_star_diff = abs(table.eta_B_star_norm(nc_idx) - eta_B_star_MAP(l));
        min_eta_B_star_diff = min(abs(table.eta_B_star_norm(nc_idx) - eta_B_star_MAP(l)));
        row_MAP_idx = find(eta_B_star_diff == min_eta_B_star_diff, 1);

        names = 'cycloaromatics/cycloaromatics_Names.csv';
        table_names = readtable(names, 'VariableNamingRule', 'preserve');

        species_MAP = table_names.Name(nc_idx(1) + row_MAP_idx - 1);
    end

    species_MAP = strrep(species_MAP, '-', '');
    strWithoutCommas = strrep(species_MAP, ',', '');

    strWithWords = strWithoutCommas;
    digitWords = {'zero', 'one', 'two', 'three', 'four', 'five', 'six', 'seven', 'eight', 'nine'};

    for i = 1:10
        digit = num2str(i-1);
        word = digitWords{i};
        strWithWords = strrep(strWithWords, digit, word);
    end

    species_MAP = strWithWords;

    fuel_name = strrep(fuel_name, ' ', '');
    fuel_name = strrep(fuel_name, '-', '');

    % Specify the folder path and filename
    folderPathSpeciesl = strcat('./OpenFOAM/', fuel_name, '/', species_MAP{1}, '/');
    % Create the folder if it doesn't exist
    if ~isfolder(folderPathSpeciesl)
        mkdir(folderPathSpeciesl);
    end
    filenameSpeciesl_C = strcat(folderPathSpeciesl, species_MAP{1}, '.C');
    filenameSpeciesl_H = strcat(folderPathSpeciesl, species_MAP{1}, '.H');
    filenameSpeciesl_IH = strcat(folderPathSpeciesl, species_MAP{1}, 'I.H');

    fid = fopen(filenameSpeciesl_C, 'w'); % Open the .C file for writing

    % Write the file header
    fprintf(fid, '/*---------------------------------------------------------------------------*\\\n');
    fprintf(fid, '  =========                | \n');
    fprintf(fid, '  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox \n');
    fprintf(fid, '   \\    /   O peration     | Website:  https://openfoam.org \n');
    fprintf(fid, '    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation \n');
    fprintf(fid, '     \\/     M anipulation  | \n');
    fprintf(fid, '------------------------------------------------------------------------------- \n');
    fprintf(fid, 'License \n');
    fprintf(fid, '    This file is part of OpenFOAM. \n');
    fprintf(fid, ' \n');
    fprintf(fid, '    OpenFOAM is free software: you can redistribute it and/or modify it \n');
    fprintf(fid, '    under the terms of the GNU General Public License as published by \n');
    fprintf(fid, '    the Free Software Foundation, either version 3 of the License, or \n');
    fprintf(fid, '    (at your option) any later version. \n');
    fprintf(fid, ' \n');
    fprintf(fid, '    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT \n');
    fprintf(fid, '    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or \n');
    fprintf(fid, '    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License \n');
    fprintf(fid, '    for more details. \n');
    fprintf(fid, ' \n');
    fprintf(fid, '    You should have received a copy of the GNU General Public License \n');
    fprintf(fid, '    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>. \n');
    fprintf(fid, ' \n');
    fprintf(fid, '\\*---------------------------------------------------------------------------*/\n');
    fprintf(fid, ' \n');
    fprintf(fid, ['#include "', species_MAP{1}, '.H"\n']);
    fprintf(fid, '#include "addToRunTimeSelectionTable.H"\n');
    fprintf(fid, ' \n');
    fprintf(fid, '// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * // \n');
    fprintf(fid, ' \n');
    fprintf(fid, 'namespace Foam \n');
    fprintf(fid, '{ \n');
    fprintf(fid, ['    defineTypeNameAndDebug(', species_MAP{1}, ', 0); \n']);
    fprintf(fid, ['    addToRunTimeSelectionTable(liquidProperties, ' species_MAP{1}, ',); \n']);
    fprintf(fid, ['    addToRunTimeSelectionTable(liquidProperties, ', species_MAP{1}, ', dictionary); \n']);
    fprintf(fid, '} \n');
    fprintf(fid, ' \n');
    fprintf(fid, '// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * // \n');
    fprintf(fid, ' \n');
    fprintf(fid, ['Foam::', species_MAP{1}, '::', species_MAP{1}, '() \n']);
    fprintf(fid, ': \n');
    fprintf(fid, '    liquidProperties \n');
    fprintf(fid, '    ( \n');
    W = table.W(nc_idx(1) + row_MAP_idx - 1);
    fprintf(fid, ['        ', num2str(W), ',\n']); % molecular weight [g/mol]
    Tc = table.Tc(nc_idx(1) + row_MAP_idx - 1);
    fprintf(fid, ['        ', num2str(Tc), ',\n']); % critical temperature [K]
    Pc = table.Pc(nc_idx(1) + row_MAP_idx - 1)*1e+5;
    fprintf(fid, ['        ', num2str(Pc), ',\n']); % critical pressure [Pa]
    Vc = table.Vc(nc_idx(1) + row_MAP_idx - 1)*1e+3;
    fprintf(fid, ['        ', num2str(Vc), ',\n']); % critical volume [m^3/kmol]
    Zc = table.Zc(nc_idx(1) + row_MAP_idx - 1);
    fprintf(fid, ['        ', num2str(Zc), ',\n']); % critical compressibility factor [-]
    fprintf(fid, ['        ', 'TBD', ',\n']); % triple point temperature [K]
    fprintf(fid, ['        ', 'TBD', ',\n']); % triple point pressure [Pa]
    Arho = table.Arho(nc_idx(1) + row_MAP_idx - 1);
    Brho = table.Brho(nc_idx(1) + row_MAP_idx - 1);
    Crho = table.Crho(nc_idx(1) + row_MAP_idx - 1);
    Amu = table.Amu(nc_idx(1) + row_MAP_idx - 1);
    Bmu = table.Bmu(nc_idx(1) + row_MAP_idx - 1);
    Cmu = table.Cmu(nc_idx(1) + row_MAP_idx - 1);
    Dmu = table.Dmu(nc_idx(1) + row_MAP_idx - 1);
    Ak = table.Ak(nc_idx(1) + row_MAP_idx - 1);
    Bk = table.Bk(nc_idx(1) + row_MAP_idx - 1);
    Ck = table.Ck(nc_idx(1) + row_MAP_idx - 1);
    Ac = table.Ac(nc_idx(1) + row_MAP_idx - 1);
    Bc = table.Bc(nc_idx(1) + row_MAP_idx - 1);
    Cc = table.Cc(nc_idx(1) + row_MAP_idx - 1);
    Dc = table.Dc(nc_idx(1) + row_MAP_idx - 1);
    Avap = table.Avap(nc_idx(1) + row_MAP_idx - 1);
    Bvap = table.Bvap(nc_idx(1) + row_MAP_idx - 1);
    Asat = table.Asat(nc_idx(1) + row_MAP_idx - 1);
    Bsat = table.Bsat(nc_idx(1) + row_MAP_idx - 1);
    Csat = table.Csat(nc_idx(1) + row_MAP_idx - 1);
    Dsat = table.Dsat(nc_idx(1) + row_MAP_idx - 1);
    Esat = table.Esat(nc_idx(1) + row_MAP_idx - 1);
    Asigma = table.Asigma(nc_idx(1) + row_MAP_idx - 1);
    Bsigma = table.Bsigma(nc_idx(1) + row_MAP_idx - 1);
    Tb = ThermophysicalProperties_SingleLiquid('boilingTemperature', 0, W, Tc, Arho, Brho, Crho, Amu, Bmu, Cmu, Dmu, Ak, Bk, Ck, Ac, Bc, Cc, Dc, Avap, Bvap, Asat, Bsat, Csat, Dsat, Esat, Asigma, Bsigma, 101325);
    fprintf(fid, ['        ', num2str(Tb), ',\n']); % normal boiling temperature [K]
    fprintf(fid, ['        ', 'TBD', ',\n']); % dipole moment [-]
    omega = table.omega(nc_idx(1) + row_MAP_idx - 1);
    fprintf(fid, ['        ', num2str(omega), ',\n']); % Pitzer's acentric factor [-]
    fprintf(fid, ['        ', 'TBD \n']); % solubility parameter [(J/m^3)^(1/2)]
    fprintf(fid, '), \n');
    fprintf(fid, ['    rho_(', num2str(Arho), ', ', num2str(Brho), ', ', num2str(Crho), ', ', num2str(Tc), '), \n']); % liquid density [kg/m^3]
    fprintf(fid, ['    pv_(', num2str(Asat), ', ', num2str(Bsat), ', ', num2str(Csat), ', ', num2str(Dsat), ', ', num2str(Esat), '), \n']); % vapor pressure [Pa]
    fprintf(fid, ['    hl_(', num2str(Avap), ', ', num2str(Bvap), ', ', num2str(Tc),', ', num2str(W), '), \n']); % latent heat of vaporization [J/kg]
    fprintf(fid, ['    Cp_(', num2str(Ac), ', ', num2str(Bc), ', ', num2str(Cc), ', ', num2str(Dc), ', ', num2str(W), '), \n']); % liquid specific heat capacity [J/kg/K]
    fprintf(fid, '    h_(TBD1, TBD2, TBD3, TBD4, TBD5, TBD6), \n'); % liquid enthalpy - reference at 298.15 K [J/kg]
    fprintf(fid, ['    Cpg_(TBD1, TBD2, TBD3, TBD4, TBD5, ', num2str(W), '), \n']); % vapor specific heat capacity [J/kg/K]
    fprintf(fid, '    B_(TBD1, TBD2, TBD3, TBD4, TBD5), \n'); % second virial coefficient [m^3/kg]
    fprintf(fid, ['    mu_(', num2str(Amu), ', ', num2str(Bmu), ', ', num2str(Cmu), ', ', num2str(Dmu), '), \n']); % liquid dynamic viscosity [Pa]
    fprintf(fid, '    mug_(TBD1, TBD2, TBD3), \n'); % liquid enthalpy - reference at 298.15 K [J/kg]
    fprintf(fid, ['    kappa_(', num2str(Ak), ', ', num2str(Bk), ', ', num2str(Ck), '), \n']); % liquid thermal conductivity [W/m/K]
    fprintf(fid, '    kappag_(TBD1, TBD2, TBD3), \n'); % vapor thermal conductivity [W/m/K]
    fprintf(fid, ['    sigma_(', num2str(Asigma), ', ', num2str(Bsigma), ', ', num2str(Tc), '), \n']); % surface tension [N/m]
    fprintf(fid, '    D_(TBD1, TBD2, TBD3) \n'); % vapor diffusivity in air [m^2/s]
    fprintf(fid, '{} \n');
    fprintf(fid, ' \n');
    fprintf(fid, ' \n');
    fprintf(fid, ['Foam::', species_MAP{1}, '::', species_MAP{1}, '\n']);
    fprintf(fid, '( \n');
    fprintf(fid, '    const liquidProperties& l, \n');
    fprintf(fid, '    const thermophysicalFunctions::YawsRho& density, \n');
    fprintf(fid, '    const thermophysicalFunctions::YawsPsat& vapourPressure, \n');
    fprintf(fid, '    const thermophysicalFunctions::YawsHvap& heatOfVapourisation, \n');
    fprintf(fid, '    const thermophysicalFunctions::YawsCp& heatCapacity, \n');
    fprintf(fid, '    const thermophysicalFunctions::NSRDS0& enthalpy, \n');
    fprintf(fid, '    const thermophysicalFunctions::YawsCpg& idealGasHeatCapacity, \n');
    fprintf(fid, '    const thermophysicalFunctions::NSRDS4& secondVirialCoeff, \n');
    fprintf(fid, '    const thermophysicalFunctions::YawsMu& dynamicViscosity, \n');
    fprintf(fid, '    const thermophysicalFunctions::YawsMug& vapourDynamicViscosity, \n');
    fprintf(fid, '    const thermophysicalFunctions::YawsKappa& thermalConductivity, \n');
    fprintf(fid, '    const thermophysicalFunctions::YawsKappa& vapourThermalConductivity, \n');
    fprintf(fid, '    const thermophysicalFunctions::YawsSigma& surfaceTension, \n');
    fprintf(fid, '    const thermophysicalFunctions::YawsD& vapourDiffusivity \n');
    fprintf(fid, ') \n');
    fprintf(fid, ': \n');
    fprintf(fid, '    liquidProperties(l), \n');
    fprintf(fid, '    rho_(density), \n');
    fprintf(fid, '    pv_(vapourPressure), \n');
    fprintf(fid, '    hl_(heatOfVapourisation), \n');
    fprintf(fid, '    Cp_(heatCapacity), \n');
    fprintf(fid, '    h_(enthalpy), \n');
    fprintf(fid, '    Cpg_(idealGasHeatCapacity), \n');
    fprintf(fid, '    B_(secondVirialCoeff), \n');
    fprintf(fid, '    mu_(dynamicViscosity), \n');
    fprintf(fid, '    mug_(vapourDynamicViscosity), \n');
    fprintf(fid, '    kappa_(thermalConductivity), \n');
    fprintf(fid, '    kappag_(vapourThermalConductivity), \n');
    fprintf(fid, '    sigma_(surfaceTension), \n');
    fprintf(fid, '    D_(vapourDiffusivity) \n');
    fprintf(fid, '{} \n');
    fprintf(fid, ' \n');
    fprintf(fid, ' \n');
    fprintf(fid, ['Foam::', species_MAP{1}, '::', species_MAP{1}, '(const dictionary& dict) \n']);
    fprintf(fid, ': \n');
    fprintf(fid, ['    ', species_MAP{1}, '() \n']);
    fprintf(fid, '{ \n');
    fprintf(fid, '    readIfPresent(*this, dict); \n');
    fprintf(fid, '} \n');
    fprintf(fid, ' \n');
    fprintf(fid, ' \n');
    fprintf(fid, '// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * // \n');
    fprintf(fid, ' \n');
    fprintf(fid, ['void Foam::', species_MAP{1}, '::write(Ostream& os) const \n']);
    fprintf(fid, '{ \n');
    fprintf(fid, '    liquidProperties::write(*this, os); \n');
    fprintf(fid, '} \n');
    fprintf(fid, ' \n');
    fprintf(fid, ' \n');
    fprintf(fid, '// ************************************************************************* // \n');

    fclose(fid); % Close the .C file

    fid = fopen(filenameSpeciesl_H, 'w'); % Open the .H file for writing

    % Write the file header
    fprintf(fid, '/*---------------------------------------------------------------------------*\\\n');
    fprintf(fid, '  =========                | \n');
    fprintf(fid, '  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox \n');
    fprintf(fid, '   \\    /   O peration     | Website:  https://openfoam.org \n');
    fprintf(fid, '    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation \n');
    fprintf(fid, '     \\/     M anipulation  | \n');
    fprintf(fid, '------------------------------------------------------------------------------- \n');
    fprintf(fid, 'License \n');
    fprintf(fid, '    This file is part of OpenFOAM. \n');
    fprintf(fid, ' \n');
    fprintf(fid, '    OpenFOAM is free software: you can redistribute it and/or modify it \n');
    fprintf(fid, '    under the terms of the GNU General Public License as published by \n');
    fprintf(fid, '    the Free Software Foundation, either version 3 of the License, or \n');
    fprintf(fid, '    (at your option) any later version. \n');
    fprintf(fid, ' \n');
    fprintf(fid, '    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT \n');
    fprintf(fid, '    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or \n');
    fprintf(fid, '    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License \n');
    fprintf(fid, '    for more details. \n');
    fprintf(fid, ' \n');
    fprintf(fid, '    You should have received a copy of the GNU General Public License \n');
    fprintf(fid, '    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>. \n');
    fprintf(fid, ' \n');
    fprintf(fid, 'Class \n');
    fprintf(fid, ['    Foam::', species_MAP{1}, ' \n']);
    fprintf(fid, ' \n');
    fprintf(fid, 'Description \n');
    fprintf(fid, ['    ', species_MAP{1}, ' \n']);
    fprintf(fid, 'SourceFiles \n');
    fprintf(fid, ['    ', species_MAP{1}, '.C\n']);
    fprintf(fid, ' \n');
    fprintf(fid, '\\*---------------------------------------------------------------------------*/ \n');
    fprintf(fid, ' \n');
    fprintf(fid, ['#ifndef ', species_MAP{1}, '_H \n']);
    fprintf(fid, ['#define ', species_MAP{1}, '_H \n']);
    fprintf(fid, ' \n');
    fprintf(fid, '#include "liquidProperties.H" \n');
    fprintf(fid, '#include "YawsRhoThermophysicalFunction.H" \n');
    fprintf(fid, '#include "YawsCpThermophysicalFunction.H" \n');
    fprintf(fid, '#include "YawsCpgThermophysicalFunction.H" \n');
    fprintf(fid, '#include "YawsMuThermophysicalFunction.H" \n');
    fprintf(fid, '#include "YawsMugThermophysicalFunction.H" \n');
    fprintf(fid, '#include "YawsKappaThermophysicalFunction.H" \n');
    fprintf(fid, '#include "YawsSigmaThermophysicalFunction.H" \n');
    fprintf(fid, '#include "YawsDThermophysicalFunction.H" \n');
    fprintf(fid, '#include "YawsPsatThermophysicalFunction.H" \n');
    fprintf(fid, '#include "YawsHvapThermophysicalFunction.H" \n');
    fprintf(fid, '#include "NSRDS0ThermophysicalFunction.H" \n');
    fprintf(fid, '#include "NSRDS4ThermophysicalFunction.H" \n');
    fprintf(fid, ' \n');
    fprintf(fid, '// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // \n');
    fprintf(fid, ' \n');
    fprintf(fid, 'namespace Foam \n');
    fprintf(fid, '{ \n');
    fprintf(fid, ' \n');
    fprintf(fid, '/*---------------------------------------------------------------------------*\\ \n');
    fprintf(fid, ['                           Class ', species_MAP{1}, ' Declaration \n']);
    fprintf(fid, '\\*---------------------------------------------------------------------------*/ \n');
    fprintf(fid, ' \n');
    fprintf(fid, ['class ', species_MAP{1}, ' \n']);
    fprintf(fid, ': \n');
    fprintf(fid, '    public liquidProperties \n');
    fprintf(fid, '{ \n');
    fprintf(fid, '    // Private Data \n');
    fprintf(fid, ' \n');
    fprintf(fid, '        thermophysicalFunctions::YawsRho rho_; \n');
    fprintf(fid, '        thermophysicalFunctions::YawsPsat pv_; \n');
    fprintf(fid, '        thermophysicalFunctions::YawsHvap hl_; \n');
    fprintf(fid, '        thermophysicalFunctions::YawsCp Cp_; \n');
    fprintf(fid, '        thermophysicalFunctions::NSRDS0 h_; \n');
    fprintf(fid, '        thermophysicalFunctions::YawsCpg Cpg_; \n');
    fprintf(fid, '        thermophysicalFunctions::NSRDS4 B_; \n');
    fprintf(fid, '        thermophysicalFunctions::YawsMu mu_; \n');
    fprintf(fid, '        thermophysicalFunctions::YawsMug mug_; \n');
    fprintf(fid, '        thermophysicalFunctions::YawsKappa kappa_; \n');
    fprintf(fid, '        thermophysicalFunctions::YawsKappa kappag_; \n');
    fprintf(fid, '        thermophysicalFunctions::YawsSigma sigma_; \n');
    fprintf(fid, '        thermophysicalFunctions::YawsD D_; \n');
    fprintf(fid, ' \n');
    fprintf(fid, ' \n');
    fprintf(fid, 'public: \n');
    fprintf(fid, ' \n');
    fprintf(fid, '    friend class liquidProperties; \n');
    fprintf(fid, ' \n');
    fprintf(fid, '    //- Runtime type information \n');
    fprintf(fid, ['    TypeName("', species_MAP{1}, '"); \n']);
    fprintf(fid, ' \n');
    fprintf(fid, ' \n');
    fprintf(fid, '    // Constructors \n');
    fprintf(fid, ' \n');
    fprintf(fid, '        //- Construct null \n');
    fprintf(fid, ['        ', species_MAP{1}, '(); \n']);
    fprintf(fid, ' \n');
    fprintf(fid, '        //- Construct from components \n');
    fprintf(fid, ['        ', species_MAP{1}, ' \n']);
    fprintf(fid, '        ( \n');
    fprintf(fid, '            const liquidProperties& l, \n');
    fprintf(fid, '            const thermophysicalFunctions::YawsRho& density, \n');
    fprintf(fid, '            const thermophysicalFunctions::YawsPsat& vapourPressure, \n');
    fprintf(fid, '            const thermophysicalFunctions::YawsHvap& heatOfVapourisation, \n');
    fprintf(fid, '            const thermophysicalFunctions::YawsCp& heatCapacity, \n');
    fprintf(fid, '            const thermophysicalFunctions::NSRDS0& enthalpy, \n');
    fprintf(fid, '            const thermophysicalFunctions::YawsCpg& idealGasHeatCapacity, \n');
    fprintf(fid, '            const thermophysicalFunctions::NSRDS4& secondVirialCoeff, \n');
    fprintf(fid, '            const thermophysicalFunctions::YawsMu& dynamicViscosity, \n');
    fprintf(fid, '            const thermophysicalFunctions::YawsMug& vapourDynamicViscosity, \n');
    fprintf(fid, '            const thermophysicalFunctions::YawsKappa& thermalConductivity, \n');
    fprintf(fid, '            const thermophysicalFunctions::YawsKappa& vapourThermalConductivity, \n');
    fprintf(fid, '            const thermophysicalFunctions::YawsSigma& surfaceTension, \n');
    fprintf(fid, '            const thermophysicalFunctions::YawsD& vapourDiffusivity \n');
    fprintf(fid, '        ); \n');
    fprintf(fid, ' \n');
    fprintf(fid, '        //- Construct from dictionary \n');
    fprintf(fid, ['        ', species_MAP{1}, '(const dictionary& dict); \n']);
    fprintf(fid, ' \n');
    fprintf(fid, '        //- Construct and return clone \n');
    fprintf(fid, '        virtual autoPtr<liquidProperties> clone() const \n');
    fprintf(fid, '        { \n');
    fprintf(fid, ['            return autoPtr<liquidProperties>(new ', species_MAP{1}, '(*this)); \n']);
    fprintf(fid, '        } \n');
    fprintf(fid, ' \n');
    fprintf(fid, ' \n');
    fprintf(fid, '    // Member Functions \n');
    fprintf(fid, ' \n');
    fprintf(fid, '        //- Liquid density [kg/m^3] \n');
    fprintf(fid, '        inline scalar rho(scalar p, scalar T) const; \n');
    fprintf(fid, ' \n');
    fprintf(fid, '        //- Vapour pressure [Pa] \n');
    fprintf(fid, '        inline scalar pv(scalar p, scalar T) const; \n');
    fprintf(fid, ' \n');
    fprintf(fid, '        //- Heat of vapourisation [J/kg] \n');
    fprintf(fid, '        inline scalar hl(scalar p, scalar T) const; \n');
    fprintf(fid, ' \n');
    fprintf(fid, '        //- Liquid heat capacity [J/kg/K] \n');
    fprintf(fid, '        inline scalar Cp(scalar p, scalar T) const; \n');
    fprintf(fid, ' \n');
    fprintf(fid, '        //- Liquid enthalpy [J/kg] \n');
    fprintf(fid, '        inline scalar h(scalar p, scalar T) const; \n');
    fprintf(fid, ' \n');
    fprintf(fid, '        //- Ideal gas heat capacity [J/kg/K] \n');
    fprintf(fid, '        inline scalar Cpg(scalar p, scalar T) const; \n');
    fprintf(fid, ' \n');
    fprintf(fid, '        //- Second Virial Coefficient [m^3/kg] \n');
    fprintf(fid, '        inline scalar B(scalar p, scalar T) const; \n');
    fprintf(fid, ' \n');
    fprintf(fid, '        //- Liquid viscosity [Pa s] \n');
    fprintf(fid, '        inline scalar mu(scalar p, scalar T) const; \n');
    fprintf(fid, ' \n');
    fprintf(fid, '        //- Vapour viscosity [Pa s] \n');
    fprintf(fid, '        inline scalar mug(scalar p, scalar T) const; \n');
    fprintf(fid, ' \n');
    fprintf(fid, '        //- Liquid thermal conductivity  [W/m/K] \n');
    fprintf(fid, '        inline scalar kappa(scalar p, scalar T) const; \n');
    fprintf(fid, ' \n');
    fprintf(fid, '        //- Vapour thermal conductivity  [W/m/K] \n');
    fprintf(fid, '        inline scalar kappag(scalar p, scalar T) const; \n');
    fprintf(fid, ' \n');
    fprintf(fid, '        //- Surface tension [N/m] \n');
    fprintf(fid, '        inline scalar sigma(scalar p, scalar T) const; \n');
    fprintf(fid, ' \n');
    fprintf(fid, '        //- Vapour diffusivity [m^2/s] \n');
    fprintf(fid, '        inline scalar D(scalar p, scalar T) const; \n');
    fprintf(fid, ' \n');
    fprintf(fid, '        //- Vapour diffusivity [m^2/s] with specified binary pair \n');
    fprintf(fid, '        inline scalar D(scalar p, scalar T, scalar Wb) const; \n');
    fprintf(fid, ' \n');
    fprintf(fid, ' \n');
    fprintf(fid, '    // I-O \n');
    fprintf(fid, ' \n');
    fprintf(fid, '        //- Write the function coefficients \n');
    fprintf(fid, '        void write(Ostream& os) const; \n');
    fprintf(fid, '}; \n');
    fprintf(fid, ' \n');
    fprintf(fid, ' \n');
    fprintf(fid, '// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // \n');
    fprintf(fid, ' \n');
    fprintf(fid, '} // End namespace Foam \n');
    fprintf(fid, ' \n');
    fprintf(fid, '// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // \n');
    fprintf(fid, ' \n');
    fprintf(fid, ['#include "', species_MAP{1}, 'I.H" \n']);
    fprintf(fid, ' \n');
    fprintf(fid, '// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // \n');
    fprintf(fid, ' \n');
    fprintf(fid, '#endif \n');
    fprintf(fid, ' \n');
    fprintf(fid, '// ************************************************************************* // \n');
    fprintf(fid, ' \n');
    fprintf(fid, ' \n');
    fprintf(fid, ' \n');
    fprintf(fid, ' \n');
    fprintf(fid, ' \n');
    fprintf(fid, ' \n');
    fprintf(fid, ' \n');
    fprintf(fid, ' \n');

    fclose(fid); % Close the .H file

    fid = fopen(filenameSpeciesl_IH, 'w'); % Open the I.H file for writing

    % Write the file header
    fprintf(fid, '/*---------------------------------------------------------------------------*\\\n');
    fprintf(fid, '  =========                | \n');
    fprintf(fid, '  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox \n');
    fprintf(fid, '   \\    /   O peration     | Website:  https://openfoam.org \n');
    fprintf(fid, '    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation \n');
    fprintf(fid, '     \\/     M anipulation  | \n');
    fprintf(fid, '------------------------------------------------------------------------------- \n');
    fprintf(fid, 'License \n');
    fprintf(fid, '    This file is part of OpenFOAM. \n');
    fprintf(fid, ' \n');
    fprintf(fid, '    OpenFOAM is free software: you can redistribute it and/or modify it \n');
    fprintf(fid, '    under the terms of the GNU General Public License as published by \n');
    fprintf(fid, '    the Free Software Foundation, either version 3 of the License, or \n');
    fprintf(fid, '    (at your option) any later version. \n');
    fprintf(fid, ' \n');
    fprintf(fid, '    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT \n');
    fprintf(fid, '    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or \n');
    fprintf(fid, '    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License \n');
    fprintf(fid, '    for more details. \n');
    fprintf(fid, ' \n');
    fprintf(fid, '    You should have received a copy of the GNU General Public License \n');
    fprintf(fid, '    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>. \n');
    fprintf(fid, ' \n');
    fprintf(fid, '\\*---------------------------------------------------------------------------*/ \n');
    fprintf(fid, ' \n');
    fprintf(fid, ['inline Foam::scalar Foam::', species_MAP{1}, '::rho(scalar p, scalar T) const \n']);
    fprintf(fid, '{ \n');
    fprintf(fid, '    return rho_.f(p, T); \n');
    fprintf(fid, '} \n');
    fprintf(fid, ' \n');
    fprintf(fid, ' \n');
    fprintf(fid, ['inline Foam::scalar Foam::', species_MAP{1}, '::pv(scalar p, scalar T) const \n']);
    fprintf(fid, '{ \n');
    fprintf(fid, '    return pv_.f(p, T); \n');
    fprintf(fid, '} \n');
    fprintf(fid, ' \n');
    fprintf(fid, ' \n');
    fprintf(fid, ['inline Foam::scalar Foam::', species_MAP{1}, '::hl(scalar p, scalar T) const \n']);
    fprintf(fid, '{ \n');
    fprintf(fid, '    return hl_.f(p, T); \n');
    fprintf(fid, '} \n');
    fprintf(fid, ' \n');
    fprintf(fid, ' \n');
    fprintf(fid, ['inline Foam::scalar Foam::', species_MAP{1}, '::Cp(scalar p, scalar T) const \n']);
    fprintf(fid, '{ \n');
    fprintf(fid, '    return Cp_.f(p, T); \n');
    fprintf(fid, '} \n');
    fprintf(fid, ' \n');
    fprintf(fid, ' \n');
    fprintf(fid, ['inline Foam::scalar Foam::', species_MAP{1}, '::h(scalar p, scalar T) const \n']);
    fprintf(fid, '{ \n');
    fprintf(fid, '    return h_.f(p, T); \n');
    fprintf(fid, '} \n');
    fprintf(fid, ' \n');
    fprintf(fid, ' \n');
    fprintf(fid, ['inline Foam::scalar Foam::', species_MAP{1}, '::Cpg(scalar p, scalar T) const \n']);
    fprintf(fid, '{ \n');
    fprintf(fid, '    return Cpg_.f(p, T); \n');
    fprintf(fid, '} \n');
    fprintf(fid, ' \n');
    fprintf(fid, ' \n');
    fprintf(fid, ['inline Foam::scalar Foam::', species_MAP{1}, '::B(scalar p, scalar T) const \n']);
    fprintf(fid, '{ \n');
    fprintf(fid, '    return B_.f(p, T); \n');
    fprintf(fid, '} \n');
    fprintf(fid, ' \n');
    fprintf(fid, ' \n');
    fprintf(fid, ['inline Foam::scalar Foam::', species_MAP{1}, '::mu(scalar p, scalar T) const \n']);
    fprintf(fid, '{ \n');
    fprintf(fid, '    return mu_.f(p, T); \n');
    fprintf(fid, '} \n');
    fprintf(fid, ' \n');
    fprintf(fid, ' \n');
    fprintf(fid, ['inline Foam::scalar Foam::', species_MAP{1}, '::mug(scalar p, scalar T) const \n']);
    fprintf(fid, '{ \n');
    fprintf(fid, '    return mug_.f(p, T); \n');
    fprintf(fid, '} \n');
    fprintf(fid, ' \n');
    fprintf(fid, ' \n');
    fprintf(fid, ['inline Foam::scalar Foam::', species_MAP{1}, '::kappa(scalar p, scalar T) const \n']);
    fprintf(fid, '{ \n');
    fprintf(fid, '    return kappa_.f(p, T); \n');
    fprintf(fid, '} \n');
    fprintf(fid, ' \n');
    fprintf(fid, ' \n');
    fprintf(fid, ['inline Foam::scalar Foam::', species_MAP{1}, '::kappag(scalar p, scalar T) const \n']);
    fprintf(fid, '{ \n');
    fprintf(fid, '    return kappag_.f(p, T); \n');
    fprintf(fid, '} \n');
    fprintf(fid, ' \n');
    fprintf(fid, ' \n');
    fprintf(fid, ['inline Foam::scalar Foam::', species_MAP{1}, '::sigma(scalar p, scalar T) const \n']);
    fprintf(fid, '{ \n');
    fprintf(fid, '    return sigma_.f(p, T); \n');
    fprintf(fid, '} \n');
    fprintf(fid, ' \n');
    fprintf(fid, ' \n');
    fprintf(fid, ['inline Foam::scalar Foam::', species_MAP{1}, '::D(scalar p, scalar T) const \n']);
    fprintf(fid, '{ \n');
    fprintf(fid, '    return D_.f(p, T); \n');
    fprintf(fid, '} \n');
    fprintf(fid, ' \n');
    fprintf(fid, ' \n');
    fprintf(fid, ['inline Foam::scalar Foam::', species_MAP{1}, '::D(scalar p, scalar T, scalar Wb) const \n']);
    fprintf(fid, '{ \n');
    fprintf(fid, '    return D_.f(p, T, Wb); \n');
    fprintf(fid, '} \n');
    fprintf(fid, ' \n');
    fprintf(fid, ' \n');
    fprintf(fid, '// ************************************************************************* // \n');

    fclose(fid); % Close the I.H file

end

end