%     ____                  _____ ___    ______
%    / __ )____ ___  _____ / ___//   |  / ____/
%   / __  / __ `/ / / / _ \__ \/ /| | / /_
%  / /_/ / /_/ / /_/ /  __/__/ / ___ |/ __/
% /_____/\__,_\__, /\___/____/_/  |_/_/
%             /____/
%
% BayeSAF: Emulation and Design of Sustainable Alternative Fuels
% via Bayesian Inference and Descriptors-Based Machine Learning
%
% Contributors / Copyright Notice
% © 2026 Jacopo Liberatori — jacopo.liberatori@centralesupelec.fr
% Postdoctoral Researcher @ Laboratoire EM2C, CentraleSupélec (CNRS)
%
% © 2026 Davide Cavalieri — davide.cavalieri@uniroma1.it
% Postdoctoral Researcher @ Sapienza University of Rome,
% Department of Mechanical and Aerospace Engineering (DIMA)
%
% © 2026 Matteo Blandino, Ph.D.
%
% Reference:
% J. Liberatori, D. Cavalieri, M. Blandino, M. Valorani, and P.P. Ciottoli.
% BayeSAF: Emulation and Design of Sustainable Alternative Fuels via Bayesian
% Inference and Descriptors-Based Machine Learning. Fuel 419, 138835 (2026).
% Available at: https://doi.org/10.1016/j.fuel.2026.138835.
%
% ------------------------------------------------------------------------
%
% Description:
% Hydrocarbons loads a structure array for a given hydrocarbon family and
% range of carbon-atom numbers from the corresponding CSV database. Each
% element of the returned array stores the species thermochemical and
% transport properties: nC, normalised topochemical atom index, molecular
% weight, critical temperature, Yaws' polynomial coefficients for dynamic
% viscosity, density, vapour pressure, thermal conductivity, specific heat
% capacity, heat of vaporisation, and surface tension, plus derived
% properties (DCN, flash point, freezing point, net heat of combustion).
%
% Inputs:
% 1) class   : string denoting the hydrocarbon family
%              ('nparaffins', 'isoparaffins', 'cycloparaffins',
%               'dicycloparaffins', 'alkylbenzenes',
%               'alkylnaphtalenes', 'cycloaromatics')
% 2) nc_range: vector of carbon-atom counts to include
%
% Outputs:
% 1) species: structure array for the given family and nc_range
% ------------------------------------------------------------------------

function species = Hydrocarbons(class, nc_range)

nparaffins     = readtable('nparaffins_MatlabData.csv',     'VariableNamingRule', 'preserve');
isoparaffins   = readtable('isoparaffins_MatlabData.csv',   'VariableNamingRule', 'preserve');
cycloparaffins = readtable('cycloparaffins_MatlabData.csv', 'VariableNamingRule', 'preserve');
dicycloparaffins  = readtable('dicycloparaffins_MatlabData.csv',  'VariableNamingRule', 'preserve');
alkylbenzenes     = readtable('alkylbenzenes_MatlabData.csv',     'VariableNamingRule', 'preserve');
alkylnaphtalenes  = readtable('alkylnaphtalenes_MatlabData.csv',  'VariableNamingRule', 'preserve');
cycloaromatics    = readtable('cycloaromatics_MatlabData.csv',    'VariableNamingRule', 'preserve');

% ── First pass: count total number of matching species ──────────────────
struct_size = 0;
for i = 1:numel(nc_range)
    nc_idx = get_nc_idx(class, nc_range(i), nparaffins, isoparaffins, ...
        cycloparaffins, dicycloparaffins, alkylbenzenes, alkylnaphtalenes, cycloaromatics);
    struct_size = struct_size + numel(nc_idx);
end

species = struct( ...
    'nC',            num2cell(zeros(1, struct_size)), ...
    'eta_B_star',    num2cell(zeros(1, struct_size)), ...
    'eta_B_star_norm',num2cell(zeros(1, struct_size)), ...
    'molWeight',     num2cell(zeros(1, struct_size)), ...
    'Tc',            num2cell(zeros(1, struct_size)), ...
    'coeffMu',       num2cell(zeros(1, struct_size)), ...
    'coeffRho',      num2cell(zeros(1, struct_size)), ...
    'coeffPsat',     num2cell(zeros(1, struct_size)), ...
    'coeffK',        num2cell(zeros(1, struct_size)), ...
    'coeffCl',       num2cell(zeros(1, struct_size)), ...
    'coeffHv',       num2cell(zeros(1, struct_size)), ...
    'coeffSigma',    num2cell(zeros(1, struct_size)), ...
    'DCN',           num2cell(zeros(1, struct_size)), ...
    'Tf',            num2cell(zeros(1, struct_size)), ...
    'Tfz',           num2cell(zeros(1, struct_size)), ...
    'Hc',            num2cell(zeros(1, struct_size)));

% ── Second pass: fill species fields ────────────────────────────────────
counter = 1;
for i = 1:numel(nc_range)
    nc     = nc_range(i);
    nc_idx = get_nc_idx(class, nc, nparaffins, isoparaffins, ...
        cycloparaffins, dicycloparaffins, alkylbenzenes, alkylnaphtalenes, cycloaromatics);

    tbl = get_table(class, nparaffins, isoparaffins, cycloparaffins, ...
        dicycloparaffins, alkylbenzenes, alkylnaphtalenes, cycloaromatics);

    for j = 1:numel(nc_idx)
        r = nc_idx(j);
        species(counter).nC              = tbl.nC(r);
        species(counter).eta_B_star      = tbl.eta_B_star(r);
        species(counter).eta_B_star_norm = tbl.eta_B_star_norm(r);
        species(counter).molWeight       = tbl.W(r);
        species(counter).Tc              = tbl.Tc(r);
        species(counter).coeffMu         = [tbl.Amu(r); tbl.Bmu(r); tbl.Cmu(r); tbl.Dmu(r)];
        species(counter).coeffRho        = [tbl.Arho(r); tbl.Brho(r); tbl.Crho(r)];
        species(counter).coeffPsat       = [tbl.Asat(r); tbl.Bsat(r); tbl.Csat(r); tbl.Dsat(r); tbl.Esat(r)];
        species(counter).coeffK          = [tbl.Ak(r); tbl.Bk(r); tbl.Ck(r)];
        species(counter).coeffCl         = [tbl.Ac(r); tbl.Bc(r); tbl.Cc(r); tbl.Dc(r)];
        species(counter).coeffHv         = [tbl.Avap(r); tbl.Bvap(r)];
        species(counter).coeffSigma      = [tbl.Asigma(r); tbl.Bsigma(r)];
        species(counter).DCN             = tbl.DCN(r);
        species(counter).Tf              = tbl.Tf(r);
        species(counter).Tfz             = tbl.Tfz(r);
        species(counter).Hc              = tbl.Hc(r);
        counter = counter + 1;
    end
end

end

% ── Helpers ─────────────────────────────────────────────────────────────

function nc_idx = get_nc_idx(class, nc, npar, isopar, cyclopar, dicyclopar, alkben, alknaph, cycloaro)
tbl = get_table(class, npar, isopar, cyclopar, dicyclopar, alkben, alknaph, cycloaro);
nc_idx = find(tbl{:,1} == nc);
end

function tbl = get_table(class, npar, isopar, cyclopar, dicyclopar, alkben, alknaph, cycloaro)
switch class
    case 'nparaffins',        tbl = npar;
    case 'isoparaffins',      tbl = isopar;
    case 'cycloparaffins',    tbl = cyclopar;
    case 'dicycloparaffins',  tbl = dicyclopar;
    case 'alkylbenzenes',     tbl = alkben;
    case 'alkylnaphthalenes', tbl = alknaph;
    case 'cycloaromatics',    tbl = cycloaro;
    otherwise,                error('Unknown hydrocarbon class: %s', class);
end
end
