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
% pdf_build_speed_cheap constructs vectorized function handles for the
% prior, likelihood, and posterior distributions used during the cheap
% MCMC stage. Identical to pdf_build_speed except that the likelihood
% is routed through LikelihoodFunction_speed_cheap, which evaluates only
% the first experimental distillation point.
%
% Inputs:
% 1) fullData               : cell array of experimental data matrices
% 2) families               : (1 x Nc) cell array of family strings
% 3) classes                : (1 x Nc) cell array of struct arrays
% 4) variable_names         : cell array of observable identifiers
% 5) pressure_distillation  : distillation pressure  [Pa]
% 6) LowerBound_x           : (1 x Nc) lower bounds for molar fractions  [-]
% 7) UpperBound_x           : (1 x Nc) upper bounds for molar fractions  [-]
% 8) alpha                  : (1 x Nc) Dirichlet concentration parameters  [-]
% 9) n_ranges               : (1 x Nc) cell array of admissible nC ranges  [-]
% 10) LowerBound_eta_B_star  : (1 x Nc) lower bounds for eta_B_star  [-]
% 11) UpperBound_eta_B_star  : (1 x Nc) upper bounds for eta_B_star  [-]
%
% Outputs:
% 1) prior     : function handle @(X_all, N_all, Eta_all) → log-prior
% 2) likelihood: function handle @(X_all, N_all, Eta_all) → log-likelihood
% 3) posterior : function handle @(X_all, N_all, Eta_all) → log-posterior
% ------------------------------------------------------------------------

function [prior, likelihood, posterior] = pdf_build_speed_cheap(...
    fullData, families, classes, variable_names, pressure_distillation, ...
    LowerBound_x, UpperBound_x, alpha, n_ranges, LowerBound_eta_B_star, UpperBound_eta_B_star)

prior_x       = @(X_all)   arrayfun(@(k) dirichlet_pdf(X_all(k,:),   alpha, LowerBound_x, UpperBound_x), 1:size(X_all,1))';
prior_eta     = @(Eta_all) arrayfun(@(k) uniform_pdf(Eta_all(k,:),   LowerBound_eta_B_star, UpperBound_eta_B_star), 1:size(Eta_all,1))';
prior_n       = @(N_all)   arrayfun(@(k) discrete_uniform(N_all(k,:), n_ranges), 1:size(N_all,1))';

scalar_prior  = @(X_all, N_all, Eta_all) prior_x(X_all) + prior_n(N_all) + prior_eta(Eta_all);

likelihood_handles = cell(1, numel(fullData));
for i = 1:numel(fullData)
    data = fullData{i};
    if ismember(variable_names{i}, {'molWeight','HC','DCN','flash','freezing','LHV'})
        stdData = data(:,2);
    else
        stdData = data(:,3);
    end
    variable_name = variable_names{i};
    likelihood_handles{i} = @(X_all, N_all, Eta_all) LikelihoodFunction_speed_cheap(...
        families, classes, data, stdData, X_all, N_all, Eta_all, variable_name, pressure_distillation);
end

scalar_likelihood = @(X_all, N_all, Eta_all) ...
    sum(cell2mat(cellfun(@(f) f(X_all, N_all, Eta_all), likelihood_handles, 'UniformOutput', false)), 2);

posterior = @(X_all, N_all, Eta_all) scalar_prior(X_all, N_all, Eta_all) + scalar_likelihood(X_all, N_all, Eta_all);
prior     = scalar_prior;
likelihood = scalar_likelihood;

end
