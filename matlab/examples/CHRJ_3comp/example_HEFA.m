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

% DESCRIPTION:
% The present example illustrates how to develop a three-component surrogate
% mixture emulating the HEFA-SPK POSF-6152 (CHRJ) fuel.

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

addpath('../../../database/nparaffins')
addpath('../../../database/isoparaffins')
addpath('../../../database/cycloparaffins')
addpath('../../../database/dicycloparaffins')
addpath('../../../database/alkylbenzenes')
addpath('../../../database/alkylnaphtalenes')
addpath('../../../database/cycloaromatics')

addpath('../../../exp/CHRJ_POSF6152')

%% --- REAL FUEL DATA  --- %%

% CSV files containing the experimental data about the thermophysical properties of the real fuel
Dataset = {'rho', 'nu', 'distillation', 'deltaT_dist', 'molWeight', 'HC', 'DCN', 'LHV', 'flash', 'freezing'};
% Plot legend string denoting real fuel experimental data
fuel_name = 'CHRJ POSF-6152';
% Bubble temperature of the real fuel [K] (if known, otherwise set it to 1e+18)
Tbubble_realFuel = 445.00;

%% --- INPUT DATA FOR THE BAYESIAN INFERENCE ANALYSIS --- %%

% Number of components the surrogate mixture will be made of
numComponents = 3;

% Define the hydrocarbon families the candidate surrogate components will be selected from
families = {'nparaffins', 'isoparaffins', 'isoparaffins'};

% Define the range the number of carbon atoms of the hydrocarbon families can vary within
n_ranges = {8:15, 9:12, 13:16};

% Define the range the normalized topochemical atom indices can vary within
LowerBound_eta_B_star = [0 0 0];
UpperBound_eta_B_star = [1 1 1];

% Define the range the mole fractions of each surrogate mixture component can vary within
LowerBound_molFrac = [0.1 0.1 0.1];
UpperBound_molFrac = [0.75 0.75 0.75];

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
alpha = 1*ones(1, numComponents); % unitary values for a symmetric uniform Dirichlet distribution

% Functional handles for the PDFs
[prior, likelihood, posterior] = pdf_build_speed(fullData, families, classes, variable_names, pressure_distillation, LowerBound_molFrac, UpperBound_molFrac, alpha, n_ranges, LowerBound_eta_B_star, UpperBound_eta_B_star);

idx_distillation = find(ismember(variable_names, {'deltaT_dist'}));

fullData_nodist = fullData;
fullData_nodist(idx_distillation) = [];

variables_nodist = variable_names;
variables_nodist(idx_distillation) = [];

% Functional handles for the PDFs
[~, ~, posterior_cheap] = pdf_build_speed_cheap(fullData_nodist, families, classes, variables_nodist, pressure_distillation, LowerBound_molFrac, UpperBound_molFrac, alpha, n_ranges, LowerBound_eta_B_star, UpperBound_eta_B_star);

%% --- DIFFERENTIAL EVOLUTION MARKOV CHAIN (DE-MC) TO EXPLORE POSTERIOR PDF --- %%

% +++ DE-MC parameters +++ %
maxIterations = 5000; % maximum number of iterations for each chain
t_burnin = 500; % number of iterations to be discarded as initial burn-in
t_check = 2*t_burnin; % number of iterations to start R-hat convergence check from
scaling_factor_X = 0.1; % scaling factor for the jump rate concerning molar fractions
scaling_factor_nc = 1.0; % scaling factor for the jump rate concerning numbers of carbon atoms
scaling_factor_eta = 1.0; % scaling factor for the jump rate concerning topochemical atom indices
noise_X = 1e-06; % noise parameter for the surrogate mixture molar fractions
noise_nc = 1e-06; % noise parameter for the surrogate mixture numbers of carbon atoms
noise_eta = 1e-06; % noise par puameter for the surrogate mixture topochemical atom indices
N_chains = 2*(3*numComponents-1); % number of chains
outlier_method = 'mad'; % outlier chain detection method: 'iqr' (interquartile range) or 'mad' (median absolute deviation)
R_hat_threshold = 1.2; % threshold value for the R-hat statistic

% +++ Parallel tempering (PT) parameters +++ %
PT_switch = 'True'; % use parallel tempering ('True') or not ('False')
beta_min = 0.0001; % minimum invterse temperature
T_ladder = 'geometric'; % temperature ladder functional form: 'linear' or 'geometric'
swap_freq = 20; % swap frequency

% +++ Proposal parameters +++ %
n_pairs = 1; % maximum number of chain pairs used to propose the new sample
p_gibbs = 0.1; % probability of a Gibbs move
p_snooker = 0.1; % probability of performing a snooker jump rather than a DE move

[chain, posterior_pdf, chain_reshaped, posterior_pdf_reshaped, AR, R_hat, t_convergence]  = DifferentialEvolutionMarkovChain(families, posterior, posterior_cheap, classes, LowerBound_molFrac, UpperBound_molFrac, n_ranges, LowerBound_eta_B_star, UpperBound_eta_B_star, maxIterations, t_burnin, t_check, scaling_factor_X, scaling_factor_nc, scaling_factor_eta, noise_X, noise_nc, noise_eta, N_chains, outlier_method, R_hat_threshold, PT_switch, beta_min, T_ladder, swap_freq, n_pairs, p_gibbs, p_snooker);

%% --- POST-PROCESSING --- %%

% Return molar fractions, numbers of carbon atoms, and topochemical atom indices characterizing the maximum a posteriori surrogate mixture
[x_MAP, nc_MAP, eta_B_star_MAP] = MAPfinder(posterior_pdf_reshaped, chain_reshaped, numComponents);

% Write a text file and display information about the MAP surrogate
species_MAP = MAPwriter(families, numComponents, x_MAP, nc_MAP, eta_B_star_MAP, fuel_name);

% Visualize Bayesian inference analysis outcomes
band_percentiles = 'True'; % 90% confidence interval colored with percentiles ('True') or not ('False')
confidence_width = 0.95; % 'confidence_width'% confidence interval
[temperature_range, property_MAP, volumeFraction_MAP, percentiles_list, mean_model, percentiles_model, model_output_samples, volFrac_interp, sobol_idx] = props_pushforward(fullData, families, classes, chain, AR, R_hat, chain_reshaped, t_convergence, t_burnin, n_ranges, numComponents, variable_names, x_MAP, nc_MAP, eta_B_star_MAP, minT_array, maxT_array, pressure_distillation, fuel_name, band_percentiles, confidence_width);

%% --- SAVE VARIABLES FROM WORKSPACE --- %%

save('CHRJ_3comp.mat');