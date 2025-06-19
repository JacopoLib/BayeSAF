%     ____                  _____ ___    ______                           %
%    / __ )____ ___  _____ / ___//   |  / ____/                           %
%   / __  / __ `/ / / / _ \\__ \/ /| | / /_                               %
%  / /_/ / /_/ / /_/ /  __/__/ / ___ |/ __/                               %
% /_____/\__,_/\__, /\___/____/_/  |_/_/                                  %
%             /____/                                                      %

% ------------------------------------------------------------------------%
% Contributors / Copyright Notice
% © 2025 Jacopo Liberatori — jacopo.liberatori@centralesupelec.fr  
% Postdoctoral Researcher @ CentraleSupélec, Laboratoire EM2C (CNRS)  
%
% © 2025 Davide Cavalieri — davide.cavalieri@uniroma1.it  
% Postdoctoral Researcher @ Sapienza University of Rome,  
% Department of Mechanical and Aerospace Engineering (DIMA)  
%
% © 2025 Matteo Blandino — matteo.blandino@uniroma1.it  
% Ph.D. Student @ Sapienza University of Rome,  
% Department of Mechanical and Aerospace Engineering (DIMA)
%
% Reference:
% J. Liberatori, D. Cavalieri, M. Blandino, M. Valorani, and P.P. Ciottoli.  
% BayeSAF: Emulation and Design of Sustainable Alternative Fuels via Bayesian  
% Inference and Descriptors-Based Machine Learning. Fuel (under review), 2025.  
% Available at: https://dx.doi.org/10.2139/ssrn.5049145.
% ------------------------------------------------------------------------%

function species_MAP = MAPwriter(families, numComponents, x_MAP, nc_MAP, eta_B_star_MAP, fuel_name)

% DESCRIPTION:
% The MAPwriter function writes a text file and display information about
% the maximum a posteriori (MAP) surrogate.

% Inputs:
% 1) families      : (1 x numComponents) cell array, with the i-th cell containing a string
% that denotes the i-th hydrocarbon family under consideration
% 2) numComponents : number of surrogate mixture components  [-]
% 3) x_MAP: (1 x numComponents) array containing the molar fractions characterizing the MAP surrogate components
% 4) nc_MAP: (1 x numComponents) array containing the numbers of carbon atoms characterizing the MAP surrogate components
% 5) eta_B_star_MAP: (1 x numComponents) array containing the topochemical atom indices characterizing the MAP surrogate components
% 6) fuel_name     : string denoting the real fuel name

% Outputs:
% 1) species_MAP: (1 x numComponents) cell array, with the i-th cell containing a
% string that denotes the name of the i-th species in the MAP surrogate

species_MAP = cell(1, numComponents); % name of MAP species
Wi_MAP = zeros(1, numComponents); % initializing the array containing the molecular weights of MAP components
nc_idx = zeros(1000, numComponents);

for l = 1:numComponents
    if strcmp(families(l), 'nparaffins')
        data = 'nparaffins/nparaffins_MatlabData.csv';
        nparaffins = readtable(data, 'VariableNamingRule', 'preserve');
        nc_idx(1:numel(find(nparaffins.nC == nc_MAP(l))),l) = find(nparaffins.nC == nc_MAP(l));
        eta_B_star_diff = abs(nparaffins.eta_B_star_norm(nc_idx(1:numel(find(nparaffins.nC == nc_MAP(l))),l)) - eta_B_star_MAP(l));
        min_eta_B_star_diff = min(abs(nparaffins.eta_B_star_norm(nc_idx(1:numel(find(nparaffins.nC == nc_MAP(l))),l)) - eta_B_star_MAP(l)));
        row_MAP_idx(l) = find(eta_B_star_diff == min_eta_B_star_diff, 1);

        names = 'nparaffins/nparaffins_Names.csv';
        nparaffins_names = readtable(names, 'VariableNamingRule', 'preserve');

        species_MAP(l) = nparaffins_names.Name(nc_idx(1,l) + row_MAP_idx(l) - 1);
        Wi_MAP(l) = nparaffins.W(nc_idx(1,l) + row_MAP_idx(l) - 1);
    elseif strcmp(families(l), 'isoparaffins')
        data = 'isoparaffins/isoparaffins_MatlabData.csv';
        isoparaffins = readtable(data, 'VariableNamingRule', 'preserve');
        nc_idx(1:numel(find(isoparaffins.nC == nc_MAP(l))),l) = find(isoparaffins.nC == nc_MAP(l));
        eta_B_star_diff = abs(isoparaffins.eta_B_star_norm(nc_idx(1:numel(find(isoparaffins.nC == nc_MAP(l))),l)) - eta_B_star_MAP(l));
        min_eta_B_star_diff = min(abs(isoparaffins.eta_B_star_norm(nc_idx(1:numel(find(isoparaffins.nC == nc_MAP(l))),l)) - eta_B_star_MAP(l)));
        row_MAP_idx(l) = find(eta_B_star_diff == min_eta_B_star_diff, 1);

        names = 'isoparaffins/isoparaffins_Names.csv';
        isoparaffins_names = readtable(names, 'VariableNamingRule', 'preserve');

        species_MAP(l) = isoparaffins_names.Name(nc_idx(1,l) + row_MAP_idx(l) - 1);
        Wi_MAP(l) = isoparaffins.W(nc_idx(1,l) + row_MAP_idx(l) - 1);
    elseif strcmp(families(l), 'cycloparaffins')
        data = 'cycloparaffins/cycloparaffins_MatlabData.csv';
        cycloparaffins = readtable(data, 'VariableNamingRule', 'preserve');
        nc_idx(1:numel(find(cycloparaffins.nC == nc_MAP(l))),l) = find(cycloparaffins.nC == nc_MAP(l));
        eta_B_star_diff = abs(cycloparaffins.eta_B_star_norm(nc_idx(1:numel(find(cycloparaffins.nC == nc_MAP(l))),l)) - eta_B_star_MAP(l));
        min_eta_B_star_diff = min(abs(cycloparaffins.eta_B_star_norm(nc_idx(1:numel(find(cycloparaffins.nC == nc_MAP(l))),l)) - eta_B_star_MAP(l)));
        row_MAP_idx(l) = find(eta_B_star_diff == min_eta_B_star_diff, 1);

        names = 'cycloparaffins/cycloparaffins_Names.csv';
        cycloparaffins_names = readtable(names, 'VariableNamingRule', 'preserve');

        species_MAP(l) = cycloparaffins_names.Name(nc_idx(1,l) + row_MAP_idx(l) - 1);
        Wi_MAP(l) = cycloparaffins.W(nc_idx(1,l) + row_MAP_idx(l) - 1);
    elseif strcmp(families(l), 'dicycloparaffins')
        data = 'dicycloparaffins/dicycloparaffins_MatlabData.csv';
        dicycloparaffins = readtable(data, 'VariableNamingRule', 'preserve');
        nc_idx(1:numel(find(dicycloparaffins.nC == nc_MAP(l))),l) = find(dicycloparaffins.nC == nc_MAP(l));
        eta_B_star_diff = abs(dicycloparaffins.eta_B_star_norm(nc_idx(1:numel(find(dicycloparaffins.nC == nc_MAP(l))),l)) - eta_B_star_MAP(l));
        min_eta_B_star_diff = min(abs(dicycloparaffins.eta_B_star_norm(nc_idx(1:numel(find(dicycloparaffins.nC == nc_MAP(l))),l)) - eta_B_star_MAP(l)));
        row_MAP_idx(l) = find(eta_B_star_diff == min_eta_B_star_diff, 1);

        names = 'dicycloparaffins/dicycloparaffins_Names.csv';
        dicycloparaffins_names = readtable(names, 'VariableNamingRule', 'preserve');

        species_MAP(l) = dicycloparaffins_names.Name(nc_idx(1,l) + row_MAP_idx(l) - 1);
        Wi_MAP(l) = dicycloparaffins.W(nc_idx(1,l) + row_MAP_idx(l) - 1);
    elseif strcmp(families(l), 'alkylbenzenes')
        data = 'alkylbenzenes/alkylbenzenes_MatlabData.csv';
        alkylbenzenes = readtable(data, 'VariableNamingRule', 'preserve');
        nc_idx(1:numel(find(alkylbenzenes.nC == nc_MAP(l))),l) = find(alkylbenzenes.nC == nc_MAP(l));
        eta_B_star_diff = abs(alkylbenzenes.eta_B_star_norm(nc_idx(1:numel(find(alkylbenzenes.nC == nc_MAP(l))),l)) - eta_B_star_MAP(l));
        min_eta_B_star_diff = min(abs(alkylbenzenes.eta_B_star_norm(nc_idx(1:numel(find(alkylbenzenes.nC == nc_MAP(l))),l)) - eta_B_star_MAP(l)));
        row_MAP_idx(l) = find(eta_B_star_diff == min_eta_B_star_diff, 1);

        names = 'alkylbenzenes/alkylbenzenes_Names.csv';
        alkylbenzenes_names = readtable(names, 'VariableNamingRule', 'preserve');

        species_MAP(l) = alkylbenzenes_names.Name(nc_idx(1,l) + row_MAP_idx(l) - 1);
        Wi_MAP(l) = alkylbenzenes.W(nc_idx(1,l) + row_MAP_idx(l) - 1);
    elseif strcmp(families(l), 'alkylnaphtalenes')
        data = 'alkylnaphtalenes/alkylnaphtalenes_MatlabData.csv';
        alkylnaphtalenes = readtable(data, 'VariableNamingRule', 'preserve');
        nc_idx(1:numel(find(alkylnaphtalenes.nC == nc_MAP(l))),l) = find(alkylnaphtalenes.nC == nc_MAP(l));
        eta_B_star_diff = abs(alkylnaphtalenes.eta_B_star_norm(nc_idx(1:numel(find(alkylnaphtalenes.nC == nc_MAP(l))),l)) - eta_B_star_MAP(l));
        min_eta_B_star_diff = min(abs(alkylnaphtalenes.eta_B_star_norm(nc_idx(1:numel(find(alkylnaphtalenes.nC == nc_MAP(l))),l)) - eta_B_star_MAP(l)));
        row_MAP_idx(l) = find(eta_B_star_diff == min_eta_B_star_diff, 1);

        names = 'alkylnaphtalenes/alkylnaphtalenes_Names.csv';
        alkylnaphtalenes_names = readtable(names, 'VariableNamingRule', 'preserve');

        species_MAP(l) = alkylnaphtalenes_names.Name(nc_idx(1,l) + row_MAP_idx(l) - 1);
        Wi_MAP(l) = alkylnaphtalenes.W(nc_idx(1,l) + row_MAP_idx(l) - 1);
    elseif strcmp(families(l), 'cycloaromatics')
        data = 'cycloaromatics/cycloaromatics_MatlabData.csv';
        cycloaromatics = readtable(data, 'VariableNamingRule', 'preserve');
        nc_idx(1:numel(find(cycloaromatics.nC == nc_MAP(l))),l) = find(cycloaromatics.nC == nc_MAP(l));
        eta_B_star_diff = abs(cycloaromatics.eta_B_star_norm(nc_idx(1:numel(find(cycloaromatics.nC == nc_MAP(l))),l)) - eta_B_star_MAP(l));
        min_eta_B_star_diff = min(abs(cycloaromatics.eta_B_star_norm(nc_idx(1:numel(find(cycloaromatics.nC == nc_MAP(l))),l)) - eta_B_star_MAP(l)));
        row_MAP_idx(l) = find(eta_B_star_diff == min_eta_B_star_diff, 1);

        names = 'cycloaromatics/cycloaromatics_Names.csv';
        cycloaromatics_names = readtable(names, 'VariableNamingRule', 'preserve', 'Delimiter', ';');

        species_MAP(l) = cycloaromatics_names.Name(nc_idx(1,l) + row_MAP_idx(l) - 1);
        Wi_MAP(l) = cycloaromatics.W(nc_idx(1,l) + row_MAP_idx(l) - 1);
    end
end

% Calculate mass fractions of MAP components and mean molecular weight of the MAP surrogate mixture
y_MAP = mol2mass(x_MAP, Wi_MAP);
W_MAP = molWeight(x_MAP, Wi_MAP);

% Compute the molecular formula of the MAP surrogate mixture
nCarbon = round(sum(x_MAP.*nc_MAP),2); % number of carbon atoms denoting the molecular formula
nHydrogen = 0.0; % number of hydrogen atoms denoting the molecular formula
for l = 1:numComponents
    if strcmp(families(l), 'nparaffins') || strcmp(families(l), 'isoparaffins')
        nHydrogen = nHydrogen + x_MAP(l)*(2*nc_MAP(l)+2);
    elseif strcmp(families(l), 'cycloparaffins')
        nHydrogen = nHydrogen + x_MAP(l)*2*nc_MAP(l);
    elseif strcmp(families(l), 'dicycloparaffins')
        nHydrogen = nHydrogen + x_MAP(l)*(2*nc_MAP(l)-2);
    elseif strcmp(families(l), 'alkylbenzenes')
        nHydrogen = nHydrogen + x_MAP(l)*(2*nc_MAP(l)-6);
    elseif strcmp(families(l), 'alkylnaphtalenes')
        nHydrogen = nHydrogen + x_MAP(l)*(2*nc_MAP(l)-12);
    elseif strcmp(families(l), 'cycloaromatics')
        nHydrogen = nHydrogen + x_MAP(l)*(2*nc_MAP(l)-8);
    end
end

nHydrogen = round(nHydrogen, 2);

% Hydrogen-to-carbon ratio (H/C)
HC = round(nHydrogen/nCarbon, 3);

% +++ Calculate the derived cetane number (DCN) of the MAP surrogate mixture
rho_i = zeros(1,numComponents);
numerator = zeros(1,numComponents);
DCN_i = zeros(1,numComponents);

for i = 1:numComponents
    if strcmp(families(i), 'nparaffins')
        rho_i(i) = ThermophysicalProperties_SingleLiquid('rho', 300, nparaffins.W(nc_idx(1,i) + row_MAP_idx(i) - 1), nparaffins.Tc(nc_idx(1,i) + row_MAP_idx(i) - 1), nparaffins.Arho(nc_idx(1,i) + row_MAP_idx(i) - 1), nparaffins.Brho(nc_idx(1,i) + row_MAP_idx(i) - 1), nparaffins.Crho(nc_idx(1,i) + row_MAP_idx(i) - 1), nparaffins.Amu(nc_idx(1,i) + row_MAP_idx(i) - 1), nparaffins.Bmu(nc_idx(1,i) + row_MAP_idx(i) - 1), nparaffins.Cmu(nc_idx(1,i) + row_MAP_idx(i) - 1), nparaffins.Dmu(nc_idx(1,i) + row_MAP_idx(i) - 1), nparaffins.Ak(nc_idx(1,i) + row_MAP_idx(i) - 1), nparaffins.Bk(nc_idx(1,i) + row_MAP_idx(i) - 1), nparaffins.Ck(nc_idx(1,i) + row_MAP_idx(i) - 1), nparaffins.Ac(nc_idx(1,i) + row_MAP_idx(i) - 1), nparaffins.Bc(nc_idx(1,i) + row_MAP_idx(i) - 1), nparaffins.Cc(nc_idx(1,i) + row_MAP_idx(i) - 1), nparaffins.Dc(nc_idx(1,i) + row_MAP_idx(i) - 1), nparaffins.Avap(nc_idx(1,i) + row_MAP_idx(i) - 1), nparaffins.Bvap(nc_idx(1,i) + row_MAP_idx(i) - 1), nparaffins.Asat(nc_idx(1,i) + row_MAP_idx(i) - 1), nparaffins.Bsat(nc_idx(1,i) + row_MAP_idx(i) - 1), nparaffins.Csat(nc_idx(1,i) + row_MAP_idx(i) - 1), nparaffins.Dsat(nc_idx(1,i) + row_MAP_idx(i) - 1), nparaffins.Esat(nc_idx(1,i) + row_MAP_idx(i) - 1), nparaffins.Asigma(nc_idx(1,i) + row_MAP_idx(i) - 1), nparaffins.Bsigma(nc_idx(1,i) + row_MAP_idx(i) - 1), 101325);
        numerator(i) = x_MAP(i)*nparaffins.W(nc_idx(1,i) + row_MAP_idx(i) - 1);
        DCN_i(i) = nparaffins.DCN(nc_idx(1,i) + row_MAP_idx(i) - 1);
    elseif strcmp(families(i), 'isoparaffins')
        rho_i(i) = ThermophysicalProperties_SingleLiquid('rho', 300, isoparaffins.W(nc_idx(1,i) + row_MAP_idx(i) - 1), isoparaffins.Tc(nc_idx(1,i) + row_MAP_idx(i) - 1), isoparaffins.Arho(nc_idx(1,i) + row_MAP_idx(i) - 1), isoparaffins.Brho(nc_idx(1,i) + row_MAP_idx(i) - 1), isoparaffins.Crho(nc_idx(1,i) + row_MAP_idx(i) - 1), isoparaffins.Amu(nc_idx(1,i) + row_MAP_idx(i) - 1), isoparaffins.Bmu(nc_idx(1,i) + row_MAP_idx(i) - 1), isoparaffins.Cmu(nc_idx(1,i) + row_MAP_idx(i) - 1), isoparaffins.Dmu(nc_idx(1,i) + row_MAP_idx(i) - 1), isoparaffins.Ak(nc_idx(1,i) + row_MAP_idx(i) - 1), isoparaffins.Bk(nc_idx(1,i) + row_MAP_idx(i) - 1), isoparaffins.Ck(nc_idx(1,i) + row_MAP_idx(i) - 1), isoparaffins.Ac(nc_idx(1,i) + row_MAP_idx(i) - 1), isoparaffins.Bc(nc_idx(1,i) + row_MAP_idx(i) - 1), isoparaffins.Cc(nc_idx(1,i) + row_MAP_idx(i) - 1), isoparaffins.Dc(nc_idx(1,i) + row_MAP_idx(i) - 1), isoparaffins.Avap(nc_idx(1,i) + row_MAP_idx(i) - 1), isoparaffins.Bvap(nc_idx(1,i) + row_MAP_idx(i) - 1), isoparaffins.Asat(nc_idx(1,i) + row_MAP_idx(i) - 1), isoparaffins.Bsat(nc_idx(1,i) + row_MAP_idx(i) - 1), isoparaffins.Csat(nc_idx(1,i) + row_MAP_idx(i) - 1), isoparaffins.Dsat(nc_idx(1,i) + row_MAP_idx(i) - 1), isoparaffins.Esat(nc_idx(1,i) + row_MAP_idx(i) - 1), isoparaffins.Asigma(nc_idx(1,i) + row_MAP_idx(i) - 1), isoparaffins.Bsigma(nc_idx(1,i) + row_MAP_idx(i) - 1), 101325);
        numerator(i) = x_MAP(i)*isoparaffins.W(nc_idx(1,i) + row_MAP_idx(i) - 1);
        DCN_i(i) = isoparaffins.DCN(nc_idx(1,i) + row_MAP_idx(i) - 1);
    elseif strcmp(families(i), 'cycloparaffins')
        rho_i(i) = ThermophysicalProperties_SingleLiquid('rho', 300, cycloparaffins.W(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloparaffins.Tc(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloparaffins.Arho(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloparaffins.Brho(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloparaffins.Crho(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloparaffins.Amu(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloparaffins.Bmu(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloparaffins.Cmu(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloparaffins.Dmu(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloparaffins.Ak(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloparaffins.Bk(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloparaffins.Ck(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloparaffins.Ac(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloparaffins.Bc(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloparaffins.Cc(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloparaffins.Dc(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloparaffins.Avap(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloparaffins.Bvap(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloparaffins.Asat(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloparaffins.Bsat(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloparaffins.Csat(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloparaffins.Dsat(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloparaffins.Esat(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloparaffins.Asigma(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloparaffins.Bsigma(nc_idx(1,i) + row_MAP_idx(i) - 1), 101325);
        numerator(i) = x_MAP(i)*cycloparaffins.W(nc_idx(1,i) + row_MAP_idx(i) - 1);
        DCN_i(i) = cycloparaffins.DCN(nc_idx(1,i) + row_MAP_idx(i) - 1);
    elseif strcmp(families(i), 'dicycloparaffins')
        rho_i(i) = ThermophysicalProperties_SingleLiquid('rho', 300, dicycloparaffins.W(nc_idx(1,i) + row_MAP_idx(i) - 1), dicycloparaffins.Tc(nc_idx(1,i) + row_MAP_idx(i) - 1), dicycloparaffins.Arho(nc_idx(1,i) + row_MAP_idx(i) - 1), dicycloparaffins.Brho(nc_idx(1,i) + row_MAP_idx(i) - 1), dicycloparaffins.Crho(nc_idx(1,i) + row_MAP_idx(i) - 1), dicycloparaffins.Amu(nc_idx(1,i) + row_MAP_idx(i) - 1), dicycloparaffins.Bmu(nc_idx(1,i) + row_MAP_idx(i) - 1), dicycloparaffins.Cmu(nc_idx(1,i) + row_MAP_idx(i) - 1), dicycloparaffins.Dmu(nc_idx(1,i) + row_MAP_idx(i) - 1), dicycloparaffins.Ak(nc_idx(1,i) + row_MAP_idx(i) - 1), dicycloparaffins.Bk(nc_idx(1,i) + row_MAP_idx(i) - 1), dicycloparaffins.Ck(nc_idx(1,i) + row_MAP_idx(i) - 1), dicycloparaffins.Ac(nc_idx(1,i) + row_MAP_idx(i) - 1), dicycloparaffins.Bc(nc_idx(1,i) + row_MAP_idx(i) - 1), dicycloparaffins.Cc(nc_idx(1,i) + row_MAP_idx(i) - 1), dicycloparaffins.Dc(nc_idx(1,i) + row_MAP_idx(i) - 1), dicycloparaffins.Avap(nc_idx(1,i) + row_MAP_idx(i) - 1), dicycloparaffins.Bvap(nc_idx(1,i) + row_MAP_idx(i) - 1), dicycloparaffins.Asat(nc_idx(1,i) + row_MAP_idx(i) - 1), dicycloparaffins.Bsat(nc_idx(1,i) + row_MAP_idx(i) - 1), dicycloparaffins.Csat(nc_idx(1,i) + row_MAP_idx(i) - 1), dicycloparaffins.Dsat(nc_idx(1,i) + row_MAP_idx(i) - 1), dicycloparaffins.Esat(nc_idx(1,i) + row_MAP_idx(i) - 1), dicycloparaffins.Asigma(nc_idx(1,i) + row_MAP_idx(i) - 1), dicycloparaffins.Bsigma(nc_idx(1,i) + row_MAP_idx(i) - 1), 101325);
        numerator(i) = x_MAP(i)*dicycloparaffins.W(nc_idx(1,i) + row_MAP_idx(i) - 1);
        DCN_i(i) = dicycloparaffins.DCN(nc_idx(1,i) + row_MAP_idx(i) - 1);
    elseif strcmp(families(i), 'alkylbenzenes')
        rho_i(i) = ThermophysicalProperties_SingleLiquid('rho', 300, alkylbenzenes.W(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylbenzenes.Tc(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylbenzenes.Arho(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylbenzenes.Brho(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylbenzenes.Crho(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylbenzenes.Amu(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylbenzenes.Bmu(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylbenzenes.Cmu(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylbenzenes.Dmu(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylbenzenes.Ak(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylbenzenes.Bk(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylbenzenes.Ck(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylbenzenes.Ac(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylbenzenes.Bc(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylbenzenes.Cc(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylbenzenes.Dc(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylbenzenes.Avap(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylbenzenes.Bvap(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylbenzenes.Asat(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylbenzenes.Bsat(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylbenzenes.Csat(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylbenzenes.Dsat(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylbenzenes.Esat(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylbenzenes.Asigma(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylbenzenes.Bsigma(nc_idx(1,i) + row_MAP_idx(i) - 1), 101325);
        numerator(i) = x_MAP(i)*alkylbenzenes.W(nc_idx(1,i) + row_MAP_idx(i) - 1);
        DCN_i(i) = alkylbenzenes.DCN(nc_idx(1,i) + row_MAP_idx(i) - 1);
    elseif strcmp(families(i), 'alkylnaphtalenes')
        rho_i(i) = ThermophysicalProperties_SingleLiquid('rho', 300, alkylnaphtalenes.W(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylnaphtalenes.Tc(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylnaphtalenes.Arho(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylnaphtalenes.Brho(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylnaphtalenes.Crho(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylnaphtalenes.Amu(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylnaphtalenes.Bmu(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylnaphtalenes.Cmu(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylnaphtalenes.Dmu(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylnaphtalenes.Ak(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylnaphtalenes.Bk(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylnaphtalenes.Ck(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylnaphtalenes.Ac(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylnaphtalenes.Bc(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylnaphtalenes.Cc(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylnaphtalenes.Dc(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylnaphtalenes.Avap(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylnaphtalenes.Bvap(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylnaphtalenes.Asat(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylnaphtalenes.Bsat(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylnaphtalenes.Csat(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylnaphtalenes.Dsat(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylnaphtalenes.Esat(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylnaphtalenes.Asigma(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylnaphtalenes.Bsigma(nc_idx(1,i) + row_MAP_idx(i) - 1), 101325);
        numerator(i) = x_MAP(i)*alkylnaphtalenes.W(nc_idx(1,i) + row_MAP_idx(i) - 1);
        DCN_i(i) = alkylnaphtalenes.DCN(nc_idx(1,i) + row_MAP_idx(i) - 1);
    elseif strcmp(families(i), 'cycloaromatics')
        rho_i(i) = ThermophysicalProperties_SingleLiquid('rho', 300, cycloaromatics.W(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloaromatics.Tc(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloaromatics.Arho(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloaromatics.Brho(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloaromatics.Crho(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloaromatics.Amu(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloaromatics.Bmu(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloaromatics.Cmu(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloaromatics.Dmu(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloaromatics.Ak(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloaromatics.Bk(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloaromatics.Ck(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloaromatics.Ac(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloaromatics.Bc(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloaromatics.Cc(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloaromatics.Dc(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloaromatics.Avap(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloaromatics.Bvap(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloaromatics.Asat(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloaromatics.Bsat(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloaromatics.Csat(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloaromatics.Dsat(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloaromatics.Esat(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloaromatics.Asigma(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloaromatics.Bsigma(nc_idx(1,i) + row_MAP_idx(i) - 1), 101325);
        numerator(i) = x_MAP(i)*cycloaromatics.W(nc_idx(1,i) + row_MAP_idx(i) - 1);
        DCN_i(i) = cycloaromatics.DCN(nc_idx(1,i) + row_MAP_idx(i) - 1);
    end
end

% Compute the MAP surrogate mixture density
rho_mixture = sum(numerator, 2)./sum(numerator./ rho_i, 2);

% Compute the MAP surrogate mixture volume fractions
Vi_MAP = rho_mixture*y_MAP'./rho_i';

DCN = sum(Vi_MAP.*DCN_i');

% +++ Calculate the flash point in [K] of the MAP surrogate mixture
flash_i = zeros(1,numComponents);
BI_flash_i = zeros(1,numComponents);

for i = 1:numComponents
    if strcmp(families(i), 'nparaffins')
        flash_i(i) = nparaffins.Tf(nc_idx(1,i) + row_MAP_idx(i) - 1);
    elseif strcmp(families(i), 'isoparaffins')
        flash_i(i) = isoparaffins.Tf(nc_idx(1,i) + row_MAP_idx(i) - 1);
    elseif strcmp(families(i), 'cycloparaffins')
        flash_i(i) = cycloparaffins.Tf(nc_idx(1,i) + row_MAP_idx(i) - 1);
    elseif strcmp(families(i), 'dicycloparaffins')
        flash_i(i) = dicycloparaffins.Tf(nc_idx(1,i) + row_MAP_idx(i) - 1);
    elseif strcmp(families(i), 'alkylbenzenes')
        flash_i(i) = alkylbenzenes.Tf(nc_idx(1,i) + row_MAP_idx(i) - 1);
    elseif strcmp(families(i), 'alkylnaphtalenes')
        flash_i(i) = alkylnaphtalenes.Tf(nc_idx(1,i) + row_MAP_idx(i) - 1);
    elseif strcmp(families(i), 'cycloaromatics')
        flash_i(i) = cycloaromatics.Tf(nc_idx(1,i) + row_MAP_idx(i) - 1);
    end
    flash_F = (flash_i(i)-273.15)*9./5+32; % flash point of the i-th MAP surrogate component in [°F]
    BI_flash_i(i) = 51708*exp((log(flash_F)-2.6287)^2/(-0.91725));
end

BI_blend = sum(y_MAP.*BI_flash_i);

flash_point_F = exp(((-0.91725)*log(BI_blend/51708))^0.5+2.6287); % flash point of the MAP surrogate mixture in [°F]
flash_point = (flash_point_F-32)*5./9+273.15; % flash point of the MAP surrogate mixture in [°F]

% +++ Calculate the freezing point in [K] of the MAP surrogate mixture
rho_i = zeros(1,numComponents);
numerator = zeros(1,numComponents);
freezing_i = zeros(1,numComponents);
BI_freezing_i = zeros(1,numComponents);

for i = 1:numComponents
    if strcmp(families(i), 'nparaffins')
        rho_i(i) = ThermophysicalProperties_SingleLiquid('rho', 300, nparaffins.W(nc_idx(1,i) + row_MAP_idx(i) - 1), nparaffins.Tc(nc_idx(1,i) + row_MAP_idx(i) - 1), nparaffins.Arho(nc_idx(1,i) + row_MAP_idx(i) - 1), nparaffins.Brho(nc_idx(1,i) + row_MAP_idx(i) - 1), nparaffins.Crho(nc_idx(1,i) + row_MAP_idx(i) - 1), nparaffins.Amu(nc_idx(1,i) + row_MAP_idx(i) - 1), nparaffins.Bmu(nc_idx(1,i) + row_MAP_idx(i) - 1), nparaffins.Cmu(nc_idx(1,i) + row_MAP_idx(i) - 1), nparaffins.Dmu(nc_idx(1,i) + row_MAP_idx(i) - 1), nparaffins.Ak(nc_idx(1,i) + row_MAP_idx(i) - 1), nparaffins.Bk(nc_idx(1,i) + row_MAP_idx(i) - 1), nparaffins.Ck(nc_idx(1,i) + row_MAP_idx(i) - 1), nparaffins.Ac(nc_idx(1,i) + row_MAP_idx(i) - 1), nparaffins.Bc(nc_idx(1,i) + row_MAP_idx(i) - 1), nparaffins.Cc(nc_idx(1,i) + row_MAP_idx(i) - 1), nparaffins.Dc(nc_idx(1,i) + row_MAP_idx(i) - 1), nparaffins.Avap(nc_idx(1,i) + row_MAP_idx(i) - 1), nparaffins.Bvap(nc_idx(1,i) + row_MAP_idx(i) - 1), nparaffins.Asat(nc_idx(1,i) + row_MAP_idx(i) - 1), nparaffins.Bsat(nc_idx(1,i) + row_MAP_idx(i) - 1), nparaffins.Csat(nc_idx(1,i) + row_MAP_idx(i) - 1), nparaffins.Dsat(nc_idx(1,i) + row_MAP_idx(i) - 1), nparaffins.Esat(nc_idx(1,i) + row_MAP_idx(i) - 1), nparaffins.Asigma(nc_idx(1,i) + row_MAP_idx(i) - 1), nparaffins.Bsigma(nc_idx(1,i) + row_MAP_idx(i) - 1), 101325);
        numerator(i) = x_MAP(i)*nparaffins.W(nc_idx(1,i) + row_MAP_idx(i) - 1);
        freezing_i(i) = nparaffins.Tfz(nc_idx(1,i) + row_MAP_idx(i) - 1);
    elseif strcmp(families(i), 'isoparaffins')
        rho_i(i) = ThermophysicalProperties_SingleLiquid('rho', 300, isoparaffins.W(nc_idx(1,i) + row_MAP_idx(i) - 1), isoparaffins.Tc(nc_idx(1,i) + row_MAP_idx(i) - 1), isoparaffins.Arho(nc_idx(1,i) + row_MAP_idx(i) - 1), isoparaffins.Brho(nc_idx(1,i) + row_MAP_idx(i) - 1), isoparaffins.Crho(nc_idx(1,i) + row_MAP_idx(i) - 1), isoparaffins.Amu(nc_idx(1,i) + row_MAP_idx(i) - 1), isoparaffins.Bmu(nc_idx(1,i) + row_MAP_idx(i) - 1), isoparaffins.Cmu(nc_idx(1,i) + row_MAP_idx(i) - 1), isoparaffins.Dmu(nc_idx(1,i) + row_MAP_idx(i) - 1), isoparaffins.Ak(nc_idx(1,i) + row_MAP_idx(i) - 1), isoparaffins.Bk(nc_idx(1,i) + row_MAP_idx(i) - 1), isoparaffins.Ck(nc_idx(1,i) + row_MAP_idx(i) - 1), isoparaffins.Ac(nc_idx(1,i) + row_MAP_idx(i) - 1), isoparaffins.Bc(nc_idx(1,i) + row_MAP_idx(i) - 1), isoparaffins.Cc(nc_idx(1,i) + row_MAP_idx(i) - 1), isoparaffins.Dc(nc_idx(1,i) + row_MAP_idx(i) - 1), isoparaffins.Avap(nc_idx(1,i) + row_MAP_idx(i) - 1), isoparaffins.Bvap(nc_idx(1,i) + row_MAP_idx(i) - 1), isoparaffins.Asat(nc_idx(1,i) + row_MAP_idx(i) - 1), isoparaffins.Bsat(nc_idx(1,i) + row_MAP_idx(i) - 1), isoparaffins.Csat(nc_idx(1,i) + row_MAP_idx(i) - 1), isoparaffins.Dsat(nc_idx(1,i) + row_MAP_idx(i) - 1), isoparaffins.Esat(nc_idx(1,i) + row_MAP_idx(i) - 1), isoparaffins.Asigma(nc_idx(1,i) + row_MAP_idx(i) - 1), isoparaffins.Bsigma(nc_idx(1,i) + row_MAP_idx(i) - 1), 101325);
        numerator(i) = x_MAP(i)*isoparaffins.W(nc_idx(1,i) + row_MAP_idx(i) - 1);
        freezing_i(i) = isoparaffins.Tfz(nc_idx(1,i) + row_MAP_idx(i) - 1);
    elseif strcmp(families(i), 'cycloparaffins')
        rho_i(i) = ThermophysicalProperties_SingleLiquid('rho', 300, cycloparaffins.W(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloparaffins.Tc(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloparaffins.Arho(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloparaffins.Brho(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloparaffins.Crho(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloparaffins.Amu(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloparaffins.Bmu(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloparaffins.Cmu(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloparaffins.Dmu(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloparaffins.Ak(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloparaffins.Bk(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloparaffins.Ck(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloparaffins.Ac(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloparaffins.Bc(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloparaffins.Cc(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloparaffins.Dc(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloparaffins.Avap(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloparaffins.Bvap(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloparaffins.Asat(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloparaffins.Bsat(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloparaffins.Csat(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloparaffins.Dsat(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloparaffins.Esat(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloparaffins.Asigma(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloparaffins.Bsigma(nc_idx(1,i) + row_MAP_idx(i) - 1), 101325);
        numerator(i) = x_MAP(i)*cycloparaffins.W(nc_idx(1,i) + row_MAP_idx(i) - 1);
        freezing_i(i) = cycloparaffins.Tfz(nc_idx(1,i) + row_MAP_idx(i) - 1);
    elseif strcmp(families(i), 'dicycloparaffins')
        rho_i(i) = ThermophysicalProperties_SingleLiquid('rho', 300, dicycloparaffins.W(nc_idx(1,i) + row_MAP_idx(i) - 1), dicycloparaffins.Tc(nc_idx(1,i) + row_MAP_idx(i) - 1), dicycloparaffins.Arho(nc_idx(1,i) + row_MAP_idx(i) - 1), dicycloparaffins.Brho(nc_idx(1,i) + row_MAP_idx(i) - 1), dicycloparaffins.Crho(nc_idx(1,i) + row_MAP_idx(i) - 1), dicycloparaffins.Amu(nc_idx(1,i) + row_MAP_idx(i) - 1), dicycloparaffins.Bmu(nc_idx(1,i) + row_MAP_idx(i) - 1), dicycloparaffins.Cmu(nc_idx(1,i) + row_MAP_idx(i) - 1), dicycloparaffins.Dmu(nc_idx(1,i) + row_MAP_idx(i) - 1), dicycloparaffins.Ak(nc_idx(1,i) + row_MAP_idx(i) - 1), dicycloparaffins.Bk(nc_idx(1,i) + row_MAP_idx(i) - 1), dicycloparaffins.Ck(nc_idx(1,i) + row_MAP_idx(i) - 1), dicycloparaffins.Ac(nc_idx(1,i) + row_MAP_idx(i) - 1), dicycloparaffins.Bc(nc_idx(1,i) + row_MAP_idx(i) - 1), dicycloparaffins.Cc(nc_idx(1,i) + row_MAP_idx(i) - 1), dicycloparaffins.Dc(nc_idx(1,i) + row_MAP_idx(i) - 1), dicycloparaffins.Avap(nc_idx(1,i) + row_MAP_idx(i) - 1), dicycloparaffins.Bvap(nc_idx(1,i) + row_MAP_idx(i) - 1), dicycloparaffins.Asat(nc_idx(1,i) + row_MAP_idx(i) - 1), dicycloparaffins.Bsat(nc_idx(1,i) + row_MAP_idx(i) - 1), dicycloparaffins.Csat(nc_idx(1,i) + row_MAP_idx(i) - 1), dicycloparaffins.Dsat(nc_idx(1,i) + row_MAP_idx(i) - 1), dicycloparaffins.Esat(nc_idx(1,i) + row_MAP_idx(i) - 1), dicycloparaffins.Asigma(nc_idx(1,i) + row_MAP_idx(i) - 1), dicycloparaffins.Bsigma(nc_idx(1,i) + row_MAP_idx(i) - 1), 101325);
        numerator(i) = x_MAP(i)*dicycloparaffins.W(nc_idx(1,i) + row_MAP_idx(i) - 1);
        freezing_i(i) = dicycloparaffins.Tfz(nc_idx(1,i) + row_MAP_idx(i) - 1);
    elseif strcmp(families(i), 'alkylbenzenes')
        rho_i(i) = ThermophysicalProperties_SingleLiquid('rho', 300, alkylbenzenes.W(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylbenzenes.Tc(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylbenzenes.Arho(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylbenzenes.Brho(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylbenzenes.Crho(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylbenzenes.Amu(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylbenzenes.Bmu(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylbenzenes.Cmu(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylbenzenes.Dmu(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylbenzenes.Ak(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylbenzenes.Bk(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylbenzenes.Ck(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylbenzenes.Ac(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylbenzenes.Bc(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylbenzenes.Cc(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylbenzenes.Dc(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylbenzenes.Avap(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylbenzenes.Bvap(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylbenzenes.Asat(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylbenzenes.Bsat(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylbenzenes.Csat(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylbenzenes.Dsat(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylbenzenes.Esat(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylbenzenes.Asigma(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylbenzenes.Bsigma(nc_idx(1,i) + row_MAP_idx(i) - 1), 101325);
        numerator(i) = x_MAP(i)*alkylbenzenes.W(nc_idx(1,i) + row_MAP_idx(i) - 1);
        freezing_i(i) = alkylbenzenes.Tfz(nc_idx(1,i) + row_MAP_idx(i) - 1);
    elseif strcmp(families(i), 'alkylnaphtalenes')
        rho_i(i) = ThermophysicalProperties_SingleLiquid('rho', 300, alkylnaphtalenes.W(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylnaphtalenes.Tc(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylnaphtalenes.Arho(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylnaphtalenes.Brho(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylnaphtalenes.Crho(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylnaphtalenes.Amu(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylnaphtalenes.Bmu(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylnaphtalenes.Cmu(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylnaphtalenes.Dmu(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylnaphtalenes.Ak(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylnaphtalenes.Bk(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylnaphtalenes.Ck(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylnaphtalenes.Ac(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylnaphtalenes.Bc(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylnaphtalenes.Cc(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylnaphtalenes.Dc(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylnaphtalenes.Avap(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylnaphtalenes.Bvap(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylnaphtalenes.Asat(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylnaphtalenes.Bsat(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylnaphtalenes.Csat(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylnaphtalenes.Dsat(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylnaphtalenes.Esat(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylnaphtalenes.Asigma(nc_idx(1,i) + row_MAP_idx(i) - 1), alkylnaphtalenes.Bsigma(nc_idx(1,i) + row_MAP_idx(i) - 1), 101325);
        numerator(i) = x_MAP(i)*alkylnaphtalenes.W(nc_idx(1,i) + row_MAP_idx(i) - 1);
        freezing_i(i) = alkylnaphtalenes.Tfz(nc_idx(1,i) + row_MAP_idx(i) - 1);
    elseif strcmp(families(i), 'cycloaromatics')
        rho_i(i) = ThermophysicalProperties_SingleLiquid('rho', 300, cycloaromatics.W(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloaromatics.Tc(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloaromatics.Arho(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloaromatics.Brho(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloaromatics.Crho(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloaromatics.Amu(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloaromatics.Bmu(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloaromatics.Cmu(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloaromatics.Dmu(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloaromatics.Ak(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloaromatics.Bk(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloaromatics.Ck(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloaromatics.Ac(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloaromatics.Bc(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloaromatics.Cc(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloaromatics.Dc(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloaromatics.Avap(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloaromatics.Bvap(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloaromatics.Asat(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloaromatics.Bsat(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloaromatics.Csat(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloaromatics.Dsat(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloaromatics.Esat(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloaromatics.Asigma(nc_idx(1,i) + row_MAP_idx(i) - 1), cycloaromatics.Bsigma(nc_idx(1,i) + row_MAP_idx(i) - 1), 101325);
        numerator(i) = x_MAP(i)*cycloaromatics.W(nc_idx(1,i) + row_MAP_idx(i) - 1);
        freezing_i(i) = cycloaromatics.Tfz(nc_idx(1,i) + row_MAP_idx(i) - 1);
    end
    BI_freezing_i(i) = freezing_i(i)^(1/0.05);
end

% Compute the MAP surrogate mixture density
rho_mixture = sum(numerator, 2)./sum(numerator./ rho_i, 2);

% Compute the MAP surrogate mixture volume fractions
Vi_MAP = rho_mixture*y_MAP'./rho_i';

BI_blend = sum(Vi_MAP'.*BI_freezing_i);

freezing_point = BI_blend^0.05; % freezing point of the MAP surrogate mixture in [K]

% +++ Calculate the lower heating value in [MJ/kg] of the MAP surrogate mixture
LHV_i = zeros(1,numComponents);

for i = 1:numComponents
    if strcmp(families(i), 'nparaffins')
        LHV_i(i) = nparaffins.Hc(nc_idx(1,i) + row_MAP_idx(i) - 1);
    elseif strcmp(families(i), 'isoparaffins')
        LHV_i(i) = isoparaffins.Hc(nc_idx(1,i) + row_MAP_idx(i) - 1);
    elseif strcmp(families(i), 'cycloparaffins')
        LHV_i(i) = cycloparaffins.Hc(nc_idx(1,i) + row_MAP_idx(i) - 1);
    elseif strcmp(families(i), 'dicycloparaffins')
        LHV_i(i) = dicycloparaffins.Hc(nc_idx(1,i) + row_MAP_idx(i) - 1);
    elseif strcmp(families(i), 'alkylbenzenes')
        LHV_i(i) = alkylbenzenes.Hc(nc_idx(1,i) + row_MAP_idx(i) - 1);
    elseif strcmp(families(i), 'alkylnaphtalenes')
        LHV_i(i) = alkylnaphtalenes.Hc(nc_idx(1,i) + row_MAP_idx(i) - 1);
    elseif strcmp(families(i), 'cycloaromatics')
        LHV_i(i) = cycloaromatics.Hc(nc_idx(1,i) + row_MAP_idx(i) - 1);
    end
end

LHV = sum(y_MAP.*LHV_i); % LHV of the MAP surrogate mixture in [MJ/kg]

% +++ Calculate the pseudocritical temperature (Tc_m) and pressure (Pc_m) of the MAP surrogate mixture according to Lee-Kesler mixing rule +++ %

Tc_m = 0.0; % initialize critical temperature
Vc_m = 0.0; % initialize critical volume
omega_m = 0.0; % initialize acentric factor

for i = 1:numComponents

    for j = 1:numComponents

        % Compute critical volume and temperature of the i-th and j-th MAP surrogate mixture components
        if strcmp(families(i), 'nparaffins')
            Vc_i = nparaffins.Vc(nc_idx(1,i) + row_MAP_idx(i) - 1);
            Tc_i = nparaffins.Tc(nc_idx(1,i) + row_MAP_idx(i) - 1);
        elseif strcmp(families(i), 'isoparaffins')
            Vc_i = isoparaffins.Vc(nc_idx(1,i) + row_MAP_idx(i) - 1);
            Tc_i = isoparaffins.Tc(nc_idx(1,i) + row_MAP_idx(i) - 1);
        elseif strcmp(families(i), 'cycloparaffins')
            Vc_i = cycloparaffins.Vc(nc_idx(1,i) + row_MAP_idx(i) - 1);
            Tc_i = cycloparaffins.Tc(nc_idx(1,i) + row_MAP_idx(i) - 1);
        elseif strcmp(families(i), 'dicycloparaffins')
            Vc_i = dicycloparaffins.Vc(nc_idx(1,i) + row_MAP_idx(i) - 1);
            Tc_i = dicycloparaffins.Tc(nc_idx(1,i) + row_MAP_idx(i) - 1);
        elseif strcmp(families(i), 'alkylbenzenes')
            Vc_i = alkylbenzenes.Vc(nc_idx(1,i) + row_MAP_idx(i) - 1);
            Tc_i = alkylbenzenes.Tc(nc_idx(1,i) + row_MAP_idx(i) - 1);
        elseif strcmp(families(i), 'alkylnaphtalenes')
            Vc_i = alkylnaphtalenes.Vc(nc_idx(1,i) + row_MAP_idx(i) - 1);
            Tc_i = alkylnaphtalenes.Tc(nc_idx(1,i) + row_MAP_idx(i) - 1);
        elseif strcmp(families(i), 'cycloaromatics')
            Vc_i = cycloaromatics.Vc(nc_idx(1,i) + row_MAP_idx(i) - 1);
            Tc_i = cycloaromatics.Tc(nc_idx(1,i) + row_MAP_idx(i) - 1);
        end

        if strcmp(families(j), 'nparaffins')
            Vc_j = nparaffins.Vc(nc_idx(1,j) + row_MAP_idx(j) - 1);
            Tc_j = nparaffins.Tc(nc_idx(1,j) + row_MAP_idx(j) - 1);
        elseif strcmp(families(j), 'isoparaffins')
            Vc_j = isoparaffins.Vc(nc_idx(1,j) + row_MAP_idx(j) - 1);
            Tc_j = isoparaffins.Tc(nc_idx(1,j) + row_MAP_idx(j) - 1);
        elseif strcmp(families(j), 'cycloparaffins')
            Vc_j = cycloparaffins.Vc(nc_idx(1,j) + row_MAP_idx(j) - 1);
            Tc_j = cycloparaffins.Tc(nc_idx(1,j) + row_MAP_idx(j) - 1);
        elseif strcmp(families(j), 'dicycloparaffins')
            Vc_j = dicycloparaffins.Vc(nc_idx(1,j) + row_MAP_idx(j) - 1);
            Tc_j = dicycloparaffins.Tc(nc_idx(1,j) + row_MAP_idx(j) - 1);
        elseif strcmp(families(j), 'alkylbenzenes')
            Vc_j = alkylbenzenes.Vc(nc_idx(1,j) + row_MAP_idx(j) - 1);
            Tc_j = alkylbenzenes.Tc(nc_idx(1,j) + row_MAP_idx(j) - 1);
        elseif strcmp(families(j), 'alkylnaphtalenes')
            Vc_j = alkylnaphtalenes.Vc(nc_idx(1,j) + row_MAP_idx(j) - 1);
            Tc_j = alkylnaphtalenes.Tc(nc_idx(1,j) + row_MAP_idx(j) - 1);
        elseif strcmp(families(j), 'cycloaromatics')
            Vc_j = cycloaromatics.Vc(nc_idx(1,j) + row_MAP_idx(j) - 1);
            Tc_j = cycloaromatics.Tc(nc_idx(1,j) + row_MAP_idx(j) - 1);
        end

        Vc_ij = 1./8*(Vc_i.^(1./3) + Vc_j.^(1./3)).^3;
        Tc_ij = (Tc_i.*Tc_j).^(1./2);

        Vc_m = Vc_m + x_MAP(i).*x_MAP(j).*Vc_ij;
        Tc_m = Tc_m + x_MAP(i).*x_MAP(j).*Vc_ij.^(1./4).*Tc_ij;

    end

    % Compute acentric factor of the i-th MAP surrogate mixture component
    if strcmp(families(i), 'nparaffins')
        omega_i = nparaffins.omega(nc_idx(1,i) + row_MAP_idx(i) - 1);
    elseif strcmp(families(i), 'isoparaffins')
        omega_i = isoparaffins.omega(nc_idx(1,i) + row_MAP_idx(i) - 1);
    elseif strcmp(families(i), 'cycloparaffins')
        omega_i = cycloparaffins.omega(nc_idx(1,i) + row_MAP_idx(i) - 1);
    elseif strcmp(families(i), 'dicycloparaffins')
        omega_i = dicycloparaffins.omega(nc_idx(1,i) + row_MAP_idx(i) - 1);
    elseif strcmp(families(i), 'alkylbenzenes')
        omega_i = alkylbenzenes.omega(nc_idx(1,i) + row_MAP_idx(i) - 1);
    elseif strcmp(families(i), 'alkylnaphtalenes')
        omega_i = alkylnaphtalenes.omega(nc_idx(1,i) + row_MAP_idx(i) - 1);
    elseif strcmp(families(i), 'cycloaromatics')
        omega_i = cycloaromatics.omega(nc_idx(1,i) + row_MAP_idx(i) - 1);
    end

    omega_m = omega_m + x_MAP(i).*omega_i;

end

Tc_m = Tc_m./(Vc_m.^(1./4)); % pseudocritical temperature of the MAP surrogate mixture [K]

R_univ = 8.314; % universal gas constant [J/mol/K]

Pc_m = (0.2905 - 0.085.*omega_m).*(R_univ.*Tc_m)./Vc_m./1e+05; % pseudocritical pressure of the MAP surrogate mixture [bar]


fid = fopen('MAP.txt', 'w'); % Open the file for writing

% Write the ASCII art header
fprintf(fid, '    ____                  _____ ___    ______                               \n');
fprintf(fid, '   / __ )____ ___  _____ / ___//   |  / ____/                               \n');
fprintf(fid, '  / __  / __ `/ / / / _ \\\\__ \\/ /| | / /_                                \n');
fprintf(fid, ' / /_/ / /_/ / /_/ /  __/__/ / ___ |/ __/                                   \n');
fprintf(fid, '/_____/\\__,_/\\__, /\\___/____/_/  |_/_/                                   \n');
fprintf(fid, '            /____/                                                          \n');
fprintf(fid, '\n'); % Add an extra newline for separation
fprintf(fid, '--------------------------------------------------------------------------  \n');
fprintf(fid, 'Contributors / Copyright Notice                                             \n');
fprintf(fid, '© 2025 Jacopo Liberatori — jacopo.liberatori@centralesupelec.fr             \n');
fprintf(fid, 'Postdoctoral Researcher @ CentraleSupélec, Laboratoire EM2C (CNRS)          \n');
fprintf(fid, '\n'); % Add an extra newline for separation
fprintf(fid, '© 2025 Davide Cavalieri — davide.cavalieri@uniroma1.it                      \n');
fprintf(fid, 'Postdoctoral Researcher @ Sapienza University of Rome,                      \n');
fprintf(fid, 'Department of Mechanical and Aerospace Engineering (DIMA)                   \n');
fprintf(fid, '\n'); % Add an extra newline for separation
fprintf(fid, '© 2025 Matteo Blandino — matteo.blandino@uniroma1.it                        \n');
fprintf(fid, 'Ph.D. Student @ Sapienza University of Rome,                                \n');
fprintf(fid, 'Department of Mechanical and Aerospace Engineering (DIMA)                   \n');
fprintf(fid, '\n'); % Add an extra newline for separation
fprintf(fid, 'Reference:                                                                  \n');
fprintf(fid, 'J. Liberatori, D. Cavalieri, M. Blandino, M. Valorani, and P.P. Ciottoli.   \n');
fprintf(fid, 'BayeSAF: Emulation and Design of Sustainable Alternative Fuels via Bayesian \n');
fprintf(fid, 'Inference and Descriptors-Based Machine Learning. Fuel (under review), 2025.\n');
fprintf(fid, 'Available at: https://dx.doi.org/10.2139/ssrn.5049145.                      \n');
fprintf(fid, '-----------------------------------------------------------------------     \n');
fprintf(fid, '\n'); % Add an extra newline for separation

% Write information to the file
fprintf(fid, ['+++++ MAP surrogate for ', fuel_name, ' +++++\n']);
fprintf(fid, 'Number of MAP components: %s\n', num2str(numComponents));
fprintf(fid, 'MAP hydrocarbon classes: %s\n', strjoin(families, ', '));
fprintf(fid, 'MAP hydrocarbons: %s\n', strjoin(species_MAP, ', '));
fprintf(fid, ['MAP molecular weight: ', num2str(W_MAP), ' g/mol\n']);
fprintf(fid, ['MAP molecular formula: C', num2str(nCarbon), ' H', num2str(nHydrogen), '\n']);
fprintf(fid, 'MAP hydrogen-to-carbon ratio: %s\n', num2str(HC));
fprintf(fid, 'MAP derived cetane number: %s\n', num2str(DCN));
fprintf(fid, ['MAP flash point: ', num2str(round(flash_point,2)), ' K\n']);
fprintf(fid, ['MAP freezing point: ', num2str(round(freezing_point,2)), ' K\n']);
fprintf(fid, ['MAP lower heating value: ', num2str(round(LHV,2)), ' MJ/kg\n']);
fprintf(fid, 'MAP molar fractions: %s\n', num2str(x_MAP));
fprintf(fid, 'MAP mass fractions: %s\n', num2str(y_MAP));
fprintf(fid, 'MAP volume fractions: %s\n', num2str(Vi_MAP'));
fprintf(fid, 'MAP carbon atoms: %s\n', num2str(nc_MAP));
fprintf(fid, 'MAP topochemical atom indices: %s\n', num2str(eta_B_star_MAP));
fprintf(fid, '\n'); % Add an extra newline for separation
fprintf(fid, ['MAP pseudocritical temperature according to Lee-Kesler mixing rule: ', num2str(Tc_m), ' K\n']);
fprintf(fid, ['MAP pseudocritical pressure according to Lee-Kesler mixing rule: ', num2str(Pc_m), ' bar\n']);

fclose(fid); % Close the file

end