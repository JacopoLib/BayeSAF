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

function  species = Hydrocarbons(class, nc_range)

% DESCRIPTION:
% The Hydrocarbons function provides a structure array for a given hydrocarbon family
% and range of number of carbon atoms.

% Auxiliary parameters:
% Nd: number of experimental measurements for the thermophysical property under consideration  [-]
% Np: number of thermophysical properties targeted during the formulation
% of the surrogate mixture  [-]

% Inputs:
% 1) class   : string denoting the hydrocarbon family under consideration
% 2) nc_range: range of the number of carbon atoms

% Outputs:
% 1) species: structure array for the given hydrocarbon family and range of
% number of carbon atoms, containing detailed information about the whole
% set of chemical species (number of carbon atoms, normalized topochemical
% atom index, molecular weight [g/mol], critical temperature [K], set of
% coefficients for the Yaws' polynomials concerning liquid-phase dynamic viscosity,
% liquid-phase density, vapor pressure, liquid-phase thermal conductivity, 
% liquid-phase specific heat capacity, heat of vaporization and liquid-phase surface tension).

data_nparaffins = 'nparaffins_MatlabData.csv';
nparaffins = readtable(data_nparaffins, 'VariableNamingRule', 'preserve');

data_isoparaffins = 'isoparaffins_MatlabData.csv';
isoparaffins = readtable(data_isoparaffins, 'VariableNamingRule', 'preserve');

data_cycloparaffins = 'cycloparaffins_MatlabData.csv';
cycloparaffins = readtable(data_cycloparaffins, 'VariableNamingRule', 'preserve');

data_dicycloparaffins = 'dicycloparaffins_MatlabData.csv';
dicycloparaffins = readtable(data_dicycloparaffins, 'VariableNamingRule', 'preserve');

data_alkylbenzenes = 'alkylbenzenes_MatlabData.csv';
alkylbenzenes = readtable(data_alkylbenzenes, 'VariableNamingRule', 'preserve');

data_alkylnaphtalenes = 'alkylnaphtalenes_MatlabData.csv';
alkylnaphtalenes = readtable(data_alkylnaphtalenes, 'VariableNamingRule', 'preserve');

data_cycloaromatics = 'cycloaromatics_MatlabData.csv';
cycloaromatics = readtable(data_cycloaromatics, 'VariableNamingRule', 'preserve');

struct_size = 0;

for i = 1:numel(nc_range)

    nc = nc_range(i);

    if strcmp(class, 'nparaffins')
        nc_idx = find(nparaffins{:, 1} == nc);
    elseif strcmp(class, 'isoparaffins')
        nc_idx = find(isoparaffins{:, 1} == nc);
    elseif strcmp(class, 'cycloparaffins')
        nc_idx = find(cycloparaffins{:, 1} == nc);
    elseif strcmp(class, 'dicycloparaffins')
        nc_idx = find(dicycloparaffins{:, 1} == nc);
    elseif strcmp(class, 'alkylbenzenes')
        nc_idx = find(alkylbenzenes{:, 1} == nc);
    elseif strcmp(class, 'alkylnaphtalenes')
        nc_idx = find(alkylnaphtalenes{:, 1} == nc);
    elseif strcmp(class, 'cycloaromatics')
        nc_idx = find(cycloaromatics{:, 1} == nc);
    end

    struct_size = struct_size + size(nc_idx,1);

end

species = struct('nC', num2cell(zeros(1, struct_size)), 'eta_B_star_norm', num2cell(zeros(1, struct_size)), 'molWeight', num2cell(zeros(1, struct_size)), 'Tc', num2cell(zeros(1, struct_size)), 'coeffMu', num2cell(zeros(1, struct_size)), 'coeffRho', num2cell(zeros(1, struct_size)), 'coeffPsat', num2cell(zeros(1, struct_size)), 'coeffK', num2cell(zeros(1, struct_size)), 'coeffCl', num2cell(zeros(1, struct_size)), 'coeffHv', num2cell(zeros(1, struct_size)), 'coeffSigma', num2cell(zeros(1, struct_size)), 'DCN', num2cell(zeros(1, struct_size)), 'Tf', num2cell(zeros(1, struct_size)), 'Tfz', num2cell(zeros(1, struct_size)), 'Hc', num2cell(zeros(1, struct_size)));

counter = 1;

for i = 1:numel(nc_range)

    nc = nc_range(i);

    if strcmp(class, 'nparaffins')
        nc_idx = find(nparaffins{:, 1} == nc);
    elseif strcmp(class, 'isoparaffins')
        nc_idx = find(isoparaffins{:, 1} == nc);
    elseif strcmp(class, 'cycloparaffins')
        nc_idx = find(cycloparaffins{:, 1} == nc);
    elseif strcmp(class, 'dicycloparaffins')
        nc_idx = find(dicycloparaffins{:, 1} == nc);
    elseif strcmp(class, 'alkylbenzenes')
        nc_idx = find(alkylbenzenes{:, 1} == nc);
    elseif strcmp(class, 'alkylnaphtalenes')
        nc_idx = find(alkylnaphtalenes{:, 1} == nc);
    elseif strcmp(class, 'cycloaromatics')
        nc_idx = find(cycloaromatics{:, 1} == nc);
    end

    for j = 1:size(nc_idx,1)

        row_idx = nc_idx(j);

        if strcmp(class, 'nparaffins')

            species(counter).nC = nparaffins.nC(row_idx);
            species(counter).eta_B_star_norm = nparaffins.eta_B_star_norm(row_idx);
            species(counter).molWeight = nparaffins.W(row_idx);
            species(counter).Tc = nparaffins.Tc(row_idx);

            species(counter).coeffMu = [nparaffins.Amu(row_idx); nparaffins.Bmu(row_idx); nparaffins.Cmu(row_idx); nparaffins.Dmu(row_idx)];
            species(counter).coeffRho = [nparaffins.Arho(row_idx); nparaffins.Brho(row_idx); nparaffins.Crho(row_idx)];
            species(counter).coeffPsat = [nparaffins.Asat(row_idx); nparaffins.Bsat(row_idx); nparaffins.Csat(row_idx); nparaffins.Dsat(row_idx); nparaffins.Esat(row_idx)];
            species(counter).coeffK = [nparaffins.Ak(row_idx); nparaffins.Bk(row_idx); nparaffins.Ck(row_idx)];
            species(counter).coeffCl = [nparaffins.Ac(row_idx); nparaffins.Bc(row_idx); nparaffins.Cc(row_idx); nparaffins.Dc(row_idx)];
            species(counter).coeffHv = [nparaffins.Avap(row_idx); nparaffins.Bvap(row_idx)];
            species(counter).coeffSigma = [nparaffins.Asigma(row_idx); nparaffins.Bsigma(row_idx)];
            species(counter).DCN = nparaffins.DCN(row_idx);
            species(counter).Tf = nparaffins.Tf(row_idx);
            species(counter).Tfz = nparaffins.Tfz(row_idx);
            species(counter).Hc = nparaffins.Hc(row_idx);

        elseif strcmp(class, 'isoparaffins')

            species(counter).nC = isoparaffins.nC(row_idx);
            species(counter).eta_B_star_norm = isoparaffins.eta_B_star_norm(row_idx);
            species(counter).molWeight = isoparaffins.W(row_idx);
            species(counter).Tc = isoparaffins.Tc(row_idx);

            species(counter).coeffMu = [isoparaffins.Amu(row_idx); isoparaffins.Bmu(row_idx); isoparaffins.Cmu(row_idx); isoparaffins.Dmu(row_idx)];
            species(counter).coeffRho = [isoparaffins.Arho(row_idx); isoparaffins.Brho(row_idx); isoparaffins.Crho(row_idx)];
            species(counter).coeffPsat = [isoparaffins.Asat(row_idx); isoparaffins.Bsat(row_idx); isoparaffins.Csat(row_idx); isoparaffins.Dsat(row_idx); isoparaffins.Esat(row_idx)];
            species(counter).coeffK = [isoparaffins.Ak(row_idx); isoparaffins.Bk(row_idx); isoparaffins.Ck(row_idx)];
            species(counter).coeffCl = [isoparaffins.Ac(row_idx); isoparaffins.Bc(row_idx); isoparaffins.Cc(row_idx); isoparaffins.Dc(row_idx)];
            species(counter).coeffHv = [isoparaffins.Avap(row_idx); isoparaffins.Bvap(row_idx)];
            species(counter).coeffSigma = [isoparaffins.Asigma(row_idx); isoparaffins.Bsigma(row_idx)];
            species(counter).DCN = isoparaffins.DCN(row_idx);
            species(counter).Tf = isoparaffins.Tf(row_idx);
            species(counter).Tfz = isoparaffins.Tfz(row_idx);
            species(counter).Hc = isoparaffins.Hc(row_idx);

        elseif strcmp(class, 'cycloparaffins')

            species(counter).nC = cycloparaffins.nC(row_idx);
            species(counter).eta_B_star_norm = cycloparaffins.eta_B_star_norm(row_idx);
            species(counter).molWeight = cycloparaffins.W(row_idx);
            species(counter).Tc = cycloparaffins.Tc(row_idx);

            species(counter).coeffMu = [cycloparaffins.Amu(row_idx); cycloparaffins.Bmu(row_idx); cycloparaffins.Cmu(row_idx); cycloparaffins.Dmu(row_idx)];
            species(counter).coeffRho = [cycloparaffins.Arho(row_idx); cycloparaffins.Brho(row_idx); cycloparaffins.Crho(row_idx)];
            species(counter).coeffPsat = [cycloparaffins.Asat(row_idx); cycloparaffins.Bsat(row_idx); cycloparaffins.Csat(row_idx); cycloparaffins.Dsat(row_idx); cycloparaffins.Esat(row_idx)];
            species(counter).coeffK = [cycloparaffins.Ak(row_idx); cycloparaffins.Bk(row_idx); cycloparaffins.Ck(row_idx)];
            species(counter).coeffCl = [cycloparaffins.Ac(row_idx); cycloparaffins.Bc(row_idx); cycloparaffins.Cc(row_idx); cycloparaffins.Dc(row_idx)];
            species(counter).coeffHv = [cycloparaffins.Avap(row_idx); cycloparaffins.Bvap(row_idx)];
            species(counter).coeffSigma = [cycloparaffins.Asigma(row_idx); cycloparaffins.Bsigma(row_idx)];
            species(counter).DCN = cycloparaffins.DCN(row_idx);
            species(counter).Tf = cycloparaffins.Tf(row_idx);
            species(counter).Tfz = cycloparaffins.Tfz(row_idx);
            species(counter).Hc = cycloparaffins.Hc(row_idx);

        elseif strcmp(class, 'dicycloparaffins')

            species(counter).nC = dicycloparaffins.nC(row_idx);
            species(counter).eta_B_star_norm = dicycloparaffins.eta_B_star_norm(row_idx);
            species(counter).molWeight = dicycloparaffins.W(row_idx);
            species(counter).Tc = dicycloparaffins.Tc(row_idx);

            species(counter).coeffMu = [dicycloparaffins.Amu(row_idx); dicycloparaffins.Bmu(row_idx); dicycloparaffins.Cmu(row_idx); dicycloparaffins.Dmu(row_idx)];
            species(counter).coeffRho = [dicycloparaffins.Arho(row_idx); dicycloparaffins.Brho(row_idx); dicycloparaffins.Crho(row_idx)];
            species(counter).coeffPsat = [dicycloparaffins.Asat(row_idx); dicycloparaffins.Bsat(row_idx); dicycloparaffins.Csat(row_idx); dicycloparaffins.Dsat(row_idx); dicycloparaffins.Esat(row_idx)];
            species(counter).coeffK = [dicycloparaffins.Ak(row_idx); dicycloparaffins.Bk(row_idx); dicycloparaffins.Ck(row_idx)];
            species(counter).coeffCl = [dicycloparaffins.Ac(row_idx); dicycloparaffins.Bc(row_idx); dicycloparaffins.Cc(row_idx); dicycloparaffins.Dc(row_idx)];
            species(counter).coeffHv = [dicycloparaffins.Avap(row_idx); dicycloparaffins.Bvap(row_idx)];
            species(counter).coeffSigma = [dicycloparaffins.Asigma(row_idx); dicycloparaffins.Bsigma(row_idx)];
            species(counter).DCN = dicycloparaffins.DCN(row_idx);
            species(counter).Tf = dicycloparaffins.Tf(row_idx);
            species(counter).Tfz = dicycloparaffins.Tfz(row_idx);
            species(counter).Hc = dicycloparaffins.Hc(row_idx);

        elseif strcmp(class, 'alkylbenzenes')

            species(counter).nC = alkylbenzenes.nC(row_idx);
            species(counter).eta_B_star_norm = alkylbenzenes.eta_B_star_norm(row_idx);
            species(counter).molWeight = alkylbenzenes.W(row_idx);
            species(counter).Tc = alkylbenzenes.Tc(row_idx);

            species(counter).coeffMu = [alkylbenzenes.Amu(row_idx); alkylbenzenes.Bmu(row_idx); alkylbenzenes.Cmu(row_idx); alkylbenzenes.Dmu(row_idx)];
            species(counter).coeffRho = [alkylbenzenes.Arho(row_idx); alkylbenzenes.Brho(row_idx); alkylbenzenes.Crho(row_idx)];
            species(counter).coeffPsat = [alkylbenzenes.Asat(row_idx); alkylbenzenes.Bsat(row_idx); alkylbenzenes.Csat(row_idx); alkylbenzenes.Dsat(row_idx); alkylbenzenes.Esat(row_idx)];
            species(counter).coeffK = [alkylbenzenes.Ak(row_idx); alkylbenzenes.Bk(row_idx); alkylbenzenes.Ck(row_idx)];
            species(counter).coeffCl = [alkylbenzenes.Ac(row_idx); alkylbenzenes.Bc(row_idx); alkylbenzenes.Cc(row_idx); alkylbenzenes.Dc(row_idx)];
            species(counter).coeffHv = [alkylbenzenes.Avap(row_idx); alkylbenzenes.Bvap(row_idx)];
            species(counter).coeffSigma = [alkylbenzenes.Asigma(row_idx); alkylbenzenes.Bsigma(row_idx)];
            species(counter).DCN = alkylbenzenes.DCN(row_idx);
            species(counter).Tf = alkylbenzenes.Tf(row_idx);
            species(counter).Tfz = alkylbenzenes.Tfz(row_idx);
            species(counter).Hc = alkylbenzenes.Hc(row_idx);

        elseif strcmp(class, 'alkylnaphtalenes')

            species(counter).nC = alkylnaphtalenes.nC(row_idx);
            species(counter).eta_B_star_norm = alkylnaphtalenes.eta_B_star_norm(row_idx);
            species(counter).molWeight = alkylnaphtalenes.W(row_idx);
            species(counter).Tc = alkylnaphtalenes.Tc(row_idx);

            species(counter).coeffMu = [alkylnaphtalenes.Amu(row_idx); alkylnaphtalenes.Bmu(row_idx); alkylnaphtalenes.Cmu(row_idx); alkylnaphtalenes.Dmu(row_idx)];
            species(counter).coeffRho = [alkylnaphtalenes.Arho(row_idx); alkylnaphtalenes.Brho(row_idx); alkylnaphtalenes.Crho(row_idx)];
            species(counter).coeffPsat = [alkylnaphtalenes.Asat(row_idx); alkylnaphtalenes.Bsat(row_idx); alkylnaphtalenes.Csat(row_idx); alkylnaphtalenes.Dsat(row_idx); alkylnaphtalenes.Esat(row_idx)];
            species(counter).coeffK = [alkylnaphtalenes.Ak(row_idx); alkylnaphtalenes.Bk(row_idx); alkylnaphtalenes.Ck(row_idx)];
            species(counter).coeffCl = [alkylnaphtalenes.Ac(row_idx); alkylnaphtalenes.Bc(row_idx); alkylnaphtalenes.Cc(row_idx); alkylnaphtalenes.Dc(row_idx)];
            species(counter).coeffHv = [alkylnaphtalenes.Avap(row_idx); alkylnaphtalenes.Bvap(row_idx)];
            species(counter).coeffSigma = [alkylnaphtalenes.Asigma(row_idx); alkylnaphtalenes.Bsigma(row_idx)];
            species(counter).DCN = alkylnaphtalenes.DCN(row_idx);
            species(counter).Tf = alkylnaphtalenes.Tf(row_idx);
            species(counter).Tfz = alkylnaphtalenes.Tfz(row_idx);
            species(counter).Hc = alkylnaphtalenes.Hc(row_idx);

        elseif strcmp(class, 'cycloaromatics')

            species(counter).nC = cycloaromatics.nC(row_idx);
            species(counter).eta_B_star_norm = cycloaromatics.eta_B_star_norm(row_idx);
            species(counter).molWeight = cycloaromatics.W(row_idx);
            species(counter).Tc = cycloaromatics.Tc(row_idx);

            species(counter).coeffMu = [cycloaromatics.Amu(row_idx); cycloaromatics.Bmu(row_idx); cycloaromatics.Cmu(row_idx); cycloaromatics.Dmu(row_idx)];
            species(counter).coeffRho = [cycloaromatics.Arho(row_idx); cycloaromatics.Brho(row_idx); cycloaromatics.Crho(row_idx)];
            species(counter).coeffPsat = [cycloaromatics.Asat(row_idx); cycloaromatics.Bsat(row_idx); cycloaromatics.Csat(row_idx); cycloaromatics.Dsat(row_idx); cycloaromatics.Esat(row_idx)];
            species(counter).coeffK = [cycloaromatics.Ak(row_idx); cycloaromatics.Bk(row_idx); cycloaromatics.Ck(row_idx)];
            species(counter).coeffCl = [cycloaromatics.Ac(row_idx); cycloaromatics.Bc(row_idx); cycloaromatics.Cc(row_idx); cycloaromatics.Dc(row_idx)];
            species(counter).coeffHv = [cycloaromatics.Avap(row_idx); cycloaromatics.Bvap(row_idx)];
            species(counter).coeffSigma = [cycloaromatics.Asigma(row_idx); cycloaromatics.Bsigma(row_idx)];
            species(counter).DCN = cycloaromatics.DCN(row_idx);
            species(counter).Tf = cycloaromatics.Tf(row_idx);
            species(counter).Tfz = cycloaromatics.Tfz(row_idx);
            species(counter).Hc = cycloaromatics.Hc(row_idx);

        end

        counter = counter + 1;

    end

end

end