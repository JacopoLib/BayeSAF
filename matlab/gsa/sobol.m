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
% Python equivalent: gsa/sobol.py (sobol_indices)
%
% Description:
% sobol computes grouped first-order Sobol' sensitivity indices based on
% the posterior PDF for three parameter groups — molar fractions, carbon-
% atom numbers, and topochemical atom indices — and their influence on
% non-lumped thermophysical properties (excluding the distillation curve).
% Indices follow the Sobol' and Myshetskaya (2008) estimator.
%
% References:
% Sobol', I.M. and Myshetskaya, E.E., 2008. Monte Carlo estimators for
% small sensitivity indices. Monte Carlo Methods Appl 13(5-6), 455-465.
%
% Inputs:
% 1) classes              : (1 x numComponents) cell array of struct arrays
% 2) numComponents        : number of surrogate mixture components  [-]
% 3) temperature_range    : (1 x 25) temperature grid  [K]
% 4) chain_reshaped       : (N_post x 3*numComponents-1) posterior samples
% 5) model_output_samples : (N_post x 25) model evaluations
% 6) variable_name        : string identifier of the thermophysical property
% 7) pressure_distillation: distillation pressure  [Pa]
%
% Outputs:
% 1) sobol_idx: (3 x 25) grouped first-order Sobol' indices  [-]
% ------------------------------------------------------------------------

function sobol_idx = sobol(classes, numComponents, temperature_range, chain_reshaped, model_output_samples, variable_name, pressure_distillation)

rng('default');

[samples, samples_indices] = unique(round(chain_reshaped, 6), 'rows');
new_output = model_output_samples(samples_indices, :);

if mod(size(samples, 1), 2) ~= 0
    samples(end,:) = [];
end

permutedIndices = randperm(size(samples, 1));
indicesA = permutedIndices(1:numel(permutedIndices)/2);
indicesB = permutedIndices(numel(permutedIndices)/2+1:end);

A   = samples(indicesA, :);
B   = samples(indicesB, :);
y_A = new_output(indicesA, :);
y_B = new_output(indicesB, :);

Vy = var([y_A; y_B]);
c  = mean([y_A; y_B]);
sobol_idx     = zeros(3, numel(temperature_range));
num_temperature = numel(temperature_range);

% +++ Group 1: molar fractions +++ %
A_B = B;
A_B(:, 1:numComponents-1) = A(:, 1:numComponents-1);
index_n_eta_AB = zeros(size(A_B,1), numComponents);
parfor p = 1:size(A_B,1)
    for q = 1:numComponents
        index_n_eta_AB(p,q) = findIndex_eta(classes{q}, A_B(p, numComponents+q-1), A_B(p, 2*numComponents+q-1));
    end
end
y_AB = zeros(size(A_B,1), num_temperature);
parfor p = 1:size(A_B,1)
    for r = 1:num_temperature
        y_AB(p,r) = mix_prop(variable_name, temperature_range(r), A_B(p, 1:numComponents-1), index_n_eta_AB(p,:), pressure_distillation, classes);
    end
end
for s = 1:num_temperature
    sobol_idx(1,s) = abs(mean((y_A(:,s) - c(s)) .* (y_AB(:,s) - y_B(:,s)))) / Vy(s);
end

% +++ Group 2: carbon-atom numbers +++ %
A_B = B;
A_B(:, numComponents:2*numComponents-1) = A(:, numComponents:2*numComponents-1);
parfor p = 1:size(A_B,1)
    for q = 1:numComponents
        index_n_eta_AB(p,q) = findIndex_eta(classes{q}, A_B(p, numComponents+q-1), A_B(p, 2*numComponents+q-1));
    end
end
y_AB = zeros(size(A_B,1), num_temperature);
parfor p = 1:size(A_B,1)
    for r = 1:num_temperature
        y_AB(p,r) = mix_prop(variable_name, temperature_range(r), A_B(p, 1:numComponents-1), index_n_eta_AB(p,:), pressure_distillation, classes);
    end
end
for s = 1:num_temperature
    sobol_idx(2,s) = abs(mean((y_A(:,s) - c(s)) .* (y_AB(:,s) - y_B(:,s)))) / Vy(s);
end

% +++ Group 3: topochemical atom indices +++ %
A_B = B;
A_B(:, 2*numComponents:3*numComponents-1) = A(:, 2*numComponents:3*numComponents-1);
parfor p = 1:size(A_B,1)
    for q = 1:numComponents
        index_n_eta_AB(p,q) = findIndex_eta(classes{q}, A_B(p, numComponents+q-1), A_B(p, 2*numComponents+q-1));
    end
end
y_AB = zeros(size(A_B,1), num_temperature);
parfor p = 1:size(A_B,1)
    for r = 1:num_temperature
        y_AB(p,r) = mix_prop(variable_name, temperature_range(r), A_B(p, 1:numComponents-1), index_n_eta_AB(p,:), pressure_distillation, classes);
    end
end
for s = 1:num_temperature
    sobol_idx(3,s) = abs(mean((y_A(:,s) - c(s)) .* (y_AB(:,s) - y_B(:,s)))) / Vy(s);
end

sobol_idx = sobol_idx ./ sum(sobol_idx, 1);

end

% -------------------------------------------------------------------------
% Local function — callable from parfor (unlike nested functions)
% -------------------------------------------------------------------------
function val = mix_prop(var_name, T, X_row, idx_row, pres, classes)
    raw = ThermophysicalProperties_LiquidMixture_speed(var_name, T, X_row, classes, idx_row, pres);
    if strcmp(var_name, 'mu') || strcmp(var_name, 'nu')
        val = 1e6 * raw;
    elseif strcmp(var_name, 'specificHeat')
        val = 1e-3 * raw;
    else
        val = raw;
    end
end
