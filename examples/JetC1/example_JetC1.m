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

% DESCRIPTION:
% The present example illustrates how to develop a two-component surrogate
% mixture emulating the alcohol-to-jet C-1 POSF-11498 fuel.
% Further details are provided at https://dx.doi.org/10.2139/ssrn.5049145. 

clc
clear
close all

addpath('../../thermo-transport')
addpath('../../distillation')
addpath('../../pdfs')
addpath('../../mcmc')
addpath('../../postprocessing')
addpath('../../gsa')
addpath('../../utilities')

addpath('../../database/nparaffins')
addpath('../../database/isoparaffins')
addpath('../../database/cycloparaffins')
addpath('../../database/dicycloparaffins')
addpath('../../database/alkylbenzenes')
addpath('../../database/alkylnaphtalenes')
addpath('../../database/cycloaromatics')

addpath('../../exp/JetC1_POSF11498')

%% --- REAL FUEL DATA  --- %%

% CSV files containing the experimental data about the thermophysical properties of the real fuel
Dataset = {'rho.csv', 'nu.csv', 'distillation.csv', 'sigma.csv', 'specificHeat.csv', 'molWeight.csv', 'HC.csv', 'DCN.csv', 'flash.csv', 'freezing.csv', 'LHV.csv'};
% Plot legend string denoting real fuel experimental data
fuel_name = 'Jet C-1 POSF-11498';
% Bubble temperature of the real fuel [K]
Tbubble_realFuel = 447.45;

%% --- INPUT DATA FOR THE BAYESIAN INFERENCE ANALYSIS --- %%

% Number of components the surrogate mixture will be made of
numComponents = 2;

% Define the hydrocarbon families the candidate surrogate components will be selected from
families = {'isoparaffins', 'isoparaffins'};

% Define the range the number of carbon atoms of the hydrocarbon families can vary within
n_ranges = {12:12, 16:16};

% Define the range the normalized topochemical atom indices can vary within
LowerBound_eta_B_star = [0 0];
UpperBound_eta_B_star = [1 1];

% Define the range the mole fractions of each surrogate mixture component can vary within
LowerBound_molFrac = [0.5 0];
UpperBound_molFrac = [1 0.5];

% Construct a cell array made of structure arrays containing detailed information about the whole
% set of candidate species for each surrogate mixture component 
classes = cell(1, numel(families));
for i = 1:numel(families)
    classes{i} = Hydrocarbons(families{i}, n_ranges{i});
end

%% --- IMPORT EXPERIMENTAL DATA --- %%

% Extract the pressure the distillation curve is computed at [Pa]
pressure_distillation = p_distillation(Dataset);

% Import experimental measurements about the thermophysical properties of the real fuel
[fullData, minT_array, maxT_array, variable_names] = data_import(Dataset, families, n_ranges, Tbubble_realFuel, numComponents, pressure_distillation);

%% --- BUILD FUNCTION HANDLES FOR THE PRIOR, LIKELIHOOD, AND POSTERIOR DISTRIBUTIONS --- %%

% Define the (1 x numComponents) vector of concentration parameters for the Dirichlet distribution of molar fractions
alpha = ones(1, numComponents); % unitary values for a symmetric uniform Dirichlet distribution

% Functional handles for the PDFs
[prior, likelihood, posterior] = pdf_build(fullData, families, classes, variable_names, pressure_distillation, LowerBound_molFrac, UpperBound_molFrac, alpha, n_ranges, LowerBound_eta_B_star, UpperBound_eta_B_star);

%% --- DIFFERENTIAL EVOLUTION MARKOV CHAIN (DE-MC) TO EXPLORE POSTERIOR PDF --- %%

% +++ DE-MC parameters +++ %
maxIterations = 5000; % maximum number of iterations for each chain
t_burnin = 500; % number of iterations to be discarded as burn-in
scaling_factor = 0.025; % scaling factor for the jump rate
noise_X = 1e-06; % noise parameter for the surrogate mixture molar fractions
noise_nc = 0.00; % noise parameter for the surrogate mixture numbers of carbon atoms
noise_eta = 0.00; % noise parameter for the surrogate mixture topochemical atom indices
N_chains = 2*(3*numComponents); % number of chains
outlier_method = 'mad'; % outlier chain detection method: 'iqr' (interquartile range) or 'mad' (median absolute deviation)
R_hat_threshold = 1.1; % threshold value for the R-hat statistic

[chain, posterior_pdf, chain_reshaped, posterior_pdf_reshaped, AR, R_hat, t_convergence]  = DifferentialEvolutionMarkovChain(posterior, classes, LowerBound_molFrac, UpperBound_molFrac, n_ranges, LowerBound_eta_B_star, UpperBound_eta_B_star, maxIterations, t_burnin, scaling_factor, noise_X, noise_nc, noise_eta, N_chains, outlier_method, R_hat_threshold);

%% --- POST-PROCESSING --- %%

% Return molar fractions, numbers of carbon atoms, and topochemical atom indices characterizing the maximum a posteriori surrogate mixture
[x_MAP, nc_MAP, eta_B_star_MAP] = MAPfinder(chain_reshaped, numComponents);

% Write a text file and display information about the MAP surrogate
species_MAP = MAPwriter(families, numComponents, x_MAP, nc_MAP, eta_B_star_MAP, fuel_name);
%%
% Visualize Bayesian inference analysis outcomes
band_percentiles = 'True'; % 90% confidence interval colored with percentiles ('True') or not ('False')
confidence_width = 0.95; % 'confidence_width'% confidence interval
visualizeResults(fullData, families, classes, chain, AR, R_hat, chain_reshaped, t_convergence, t_burnin, n_ranges, numComponents, variable_names, x_MAP, nc_MAP, eta_B_star_MAP, minT_array, maxT_array, pressure_distillation, fuel_name, band_percentiles, confidence_width);

%% --- SAVE VARIABLES FROM WORKSPACE --- %%

save('test_JetC1.mat');
