function  [prior, likelihood, posterior] = pdf_build(fullData, classes, variable_names, pressure_distillation, LowerBound_x, UpperBound_x, n_ranges, LowerBound_eta_B_star, UpperBound_eta_B_star)

% DESCRIPTION:
% The pdf_build function calculates the function handles representing the
% prior, likelihood, and posterior distributions.

% Auxiliary parameters:
% Nc: number of surrogate mixture components  [-]
% Nd: number of experimental measurements for the thermophysical property under consideration  [-]
% Np: number of thermophysical properties targeted during the formulation
% of the surrogate mixture  [-]

% Inputs:
% 1) fullData             : (1 x Np) cell array, with the i-th cell being a (Nd x 3) array
% that contains the experimental measurements about the thermophysical property under
% consideration (independent variable, dependent thermophysical property, standard deviation). 
% Note that the i-th cell is a (Nd x 4) array in case the distillation curve is the i-th
% thermophysical property under consideration.
% 2) classes              : (1 x Nc) cell array, with the i-th cell containing a structure array
% for the hydrocarbon family and range of number of carbon atoms of the i-th surrogate component
% 3) variable_names       : (1 x Np) cell array, with the i-th cell containing a string
% that denotes the i-th thermophysical property under consideration
% 4) pressure_distillation: pressure the distillation curve is computed at  [Pa]
% 5) LowerBound_x         : (1 x Nc) array, with the i-th element representing the
% lower bound for the range the mole fraction of the i-th surrogate mixture
% component can vary within  [-]
% 6) UpperBound_x         : (1 x Nc) array, with the i-th element representing the
% upper bound for the range the mole fraction of the i-th surrogate mixture
% component can vary within  [-]
% 7) n_ranges             : (1 x Nc) cell array, with the i-th cell containing the range
% the number of carbon atoms of the i-th surrogate mixture component can vary within
% 8) LowerBound_eta_B_star: (1 x Nc) array, with the i-th element representing the
% lower bound for the range the normalized topochemical atom index for the i-th
% surrogate mixture component can vary within  [-]
% 9) UpperBound_eta_B_star: (1 x Nc) array, with the i-th element representing the
% upper bound for the range the normalized topochemical atom index for the i-th
% surrogate mixture component can vary within  [-]


% Outputs:
% 1) prior: function handle for the prior distribution
% 2) likelihood: function handle for the likelihood distribution
% 3) posterior: function handle for the posterior distribution

% ------------------------------------------------------------------------%
% Contributors/Copyright
% 2024 Jacopo Liberatori, jacopo.liberatori@uniroma1.it
% 2024 Davide Cavalieri,  davide.cavalieri@uniroma1.it
% Department of Mechanical and Aerospace Engineering (DIMA)
% Sapienza University of Rome
% ------------------------------------------------------------------------%

prior_informative_x = @(x)uniform_pdf(x, LowerBound_x, UpperBound_x);

% Set alpha = 1 to get a symmetric Dirichlet distribution
alpha = 1;
prior_constraint_x = @(x)dirichlet_pdf(x, alpha, LowerBound_x, UpperBound_x);

prior_x = @(x) prior_informative_x(x) + prior_constraint_x(x);

prior_informative_eta_B_star = @(eta_B_star)uniform_pdf(eta_B_star, LowerBound_eta_B_star, UpperBound_eta_B_star);

prior_informative_n = @(n) discrete_uniform(n, n_ranges);

prior = @(x, n, eta_B_star) prior_x(x) + prior_informative_n(n) + prior_informative_eta_B_star(eta_B_star);

for i = 1:numel(fullData)

    data = fullData{i};
    stdData = data(:, 3);
    variable_name = variable_names{i};

    likelihood = @(x, n, eta_B_star) LikelihoodFunction(classes, data, stdData, x, n, eta_B_star, variable_name, pressure_distillation);

    posterior = @(x, n, eta_B_star) likelihood(x, n, eta_B_star) + prior(x, n, eta_B_star);

    prior = @(x, n, eta_B_star) posterior(x, n, eta_B_star);

end

end