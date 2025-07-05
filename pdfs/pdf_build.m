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

function  [prior, likelihood, posterior] = pdf_build(fullData, families, classes, variable_names, pressure_distillation, LowerBound_x, UpperBound_x, n_ranges, LowerBound_eta_B_star, UpperBound_eta_B_star)

% DESCRIPTION:
% The pdf_build function calculates the function handles representing the
% prior, likelihood, and posterior distributions.

% Auxiliary parameters:
% Nc: number of surrogate mixture components  [-]
% Nd: number of experimental measurements for the thermophysical property under consideration  [-]
% Np: number of thermophysical properties targeted during the formulation
% of the surrogate mixture  [-]

% Inputs:
% 1) fullData              : (1 x Np) cell array, with the i-th cell being a (Nd x 3) array
% that contains the experimental measurements about the thermophysical property under
% consideration (independent variable, dependent thermophysical property, standard deviation). 
% Note that the i-th cell is a (Nd x 4) array in case the distillation curve is the i-th
% thermophysical property under consideration.
% 2) families              : (1 x Nc) cell array, with the i-th cell containing a string
% that denotes the i-th hydrocarbon family under consideration
% 3) classes               : (1 x Nc) cell array, with the i-th cell containing a structure array
% for the hydrocarbon family and range of number of carbon atoms of the i-th surrogate component
% 4) variable_names        : (1 x Np) cell array, with the i-th cell containing a string
% that denotes the i-th thermophysical property under consideration
% 5) pressure_distillation: pressure the distillation curve is computed at  [Pa]
% 6) LowerBound_x          : (1 x Nc) array, with the i-th element representing the
% lower bound for the range the mole fraction of the i-th surrogate mixture
% component can vary within  [-]
% 7) UpperBound_x          : (1 x Nc) array, with the i-th element representing the
% upper bound for the range the mole fraction of the i-th surrogate mixture
% component can vary within  [-]
% 8) alpha                 : (1 x Nc) vector of concentration parameters for the Dirichlet
% distribution of the molar fractions  [-]
% 9) n_ranges              : (1 x Nc) cell array, with the i-th cell containing the range
% the number of carbon atoms of the i-th surrogate mixture component can vary within
% 10) LowerBound_eta_B_star: (1 x Nc) array, with the i-th element representing the
% lower bound for the range the normalized topochemical atom index for the i-th
% surrogate mixture component can vary within  [-]
% 11) UpperBound_eta_B_star: (1 x Nc) array, with the i-th element representing the
% upper bound for the range the normalized topochemical atom index for the i-th
% surrogate mixture component can vary within  [-]

% Outputs:
% 1) prior: function handle for the prior distribution
% 2) likelihood: function handle for the likelihood distribution
% 3) posterior: function handle for the posterior distribution

prior_informative_x = @(x)dirichlet_pdf(x, alpha, LowerBound_x, UpperBound_x);

prior_informative_eta_B_star = @(eta_B_star)uniform_pdf(eta_B_star, LowerBound_eta_B_star, UpperBound_eta_B_star);

prior_informative_n = @(n) discrete_uniform(n, n_ranges);

prior = @(x, n, eta_B_star) prior_informative_x(x) + prior_informative_n(n) + prior_informative_eta_B_star(eta_B_star);

for i = 1:numel(fullData)

    data = fullData{i};
    if strcmp(variable_names{i},'molWeight') || strcmp(variable_names{i},'HC') || strcmp(variable_names{i},'DCN') || strcmp(variable_names{i},'flash') || strcmp(variable_names{i},'freezing') || strcmp(variable_names{i},'LHV')
        stdData = data(:, 2);
    else
        stdData = data(:, 3);
    end
    variable_name = variable_names{i};

    likelihood = @(x, n, eta_B_star) LikelihoodFunction(families, classes, data, stdData, x, n, eta_B_star, variable_name, pressure_distillation);

    posterior = @(x, n, eta_B_star) likelihood(x, n, eta_B_star) + prior(x, n, eta_B_star);

    prior = @(x, n, eta_B_star) posterior(x, n, eta_B_star);

end

end
