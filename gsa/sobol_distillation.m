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

function sobol_idx = sobol_distillation(classes, numComponents, chain_reshaped, model_output_samples, pressure_distillation)

% DESCRIPTION:
% The sobol function calculates grouped first-order Sobol' indices based on the posterior PDF about three sets of
% input parameters and their impact on the distillation curve: i) molar fractions, ii) numbers of carbon atoms, and iii) topochemical atom indices. 
% Note that first-order Sobol' indices are calculated according to the algorithm proposed in Sobol' and Myshetskaya (2008).

% --- References ---    

% Sobol', I.M. and Myshetskaya, E.E., 2008.
% Monte Carlo estimators for small sensitivity indices.
% Monte Carlo Methods Appl 13 (5-6), 455–65.

% Auxiliary parameters:
% N_chains     : number of DE-MC chains
% t_convergence: number of iterations denoting the number of iterations the convergence of the DE-MC run
% is reached at according to the R-hat statistic ( = maxIterations in case convergence is not reached)
% t_burnin     : number of iterations to be discarded as burn-in
% maxIterations: maximum number of iterations for each chain in the DE-MC algorithm

% Inputs:
% 1) classes              : (1 x numComponents) cell array, with the i-th cell containing a structure array
% for the hydrocarbon family and range of number of carbon atoms of the i-th surrogate component
% 2) numComponents        : number of surrogate mixture components  [-]
% 3) chain_reshaped       : ((t_convergence - t_burnin)*N_chains x 3*numComponents-1) array containing the
% samples from the posterior PDF from each chain discarding the burn-in period, concerning
% the molar fractions, numbers of carbon atoms, and topochemical atom indices
% 4) model_output_samples : (t_convergence - t_burnin)*N_chains x 101) array describing the temperature
% evolution along the distillation curve for each sample in chain_reshaped
% 5) pressure_distillation: pressure the distillation curve is computed at  [Pa]

% Outputs:
% 1) sobol_idx: (3 x 101) array describing the evolution of grouped first-order Sobol' indices along
% the distillation curve, i.e., for the recovered volume fraction ranging from 0 to 100 %

rng('default'); % to ensure reproducibility

% Remove duplicate samples
[samples,samples_indices]=unique(round(chain_reshaped,6),'rows');
% Remove model evaluations corresponding to duplicate samples
new_output = model_output_samples(samples_indices,:);

if mod(size(samples,1), 2) == 0
    % Do nothing
else
    samples(end,:) = [];
end

% Randomly permute the row indices
permutedIndices = randperm(size(samples,1));

% Split the permuted indices into two sets
indicesA = permutedIndices(1:numel(permutedIndices)/2);
indicesB = permutedIndices(numel(permutedIndices)/2+1:end);

% Split the sample points into A and B matrices
A = samples(indicesA, :);
B = samples(indicesB, :);
% Model evaluations for A and B
y_A = new_output(indicesA,:);
y_B = new_output(indicesB,:);

% Variance and mean of the overall model evaluations
Vy = var([y_A; y_B]);
c = mean([y_A; y_B]);
% Initialize vector containing first-order Sobol' indices
volFrac_interp = linspace(1, 99, 99);
volFrac_interp = [0.1, volFrac_interp, 99.9];
sobol_idx = zeros(3, numel(volFrac_interp));

% +++ First-order Sobol' indices for molar fractions group +++ %
A_B = B;
y_AB = zeros(size(A_B,1), numel(volFrac_interp));
A_B(:, 1:numComponents-1) = A(:, 1:numComponents-1);
index_n_eta_AB = zeros(size(A_B,1), numComponents);
parfor p = 1:size(A_B,1)
    for q = 1:numComponents
        index_n_eta_AB(p,q) = findIndex_eta(classes{q}, A_B(p, numComponents+q-1), A_B(p, 2*numComponents+q-1));
    end
end
parfor p = 1:size(A_B,1)
    [volumeFraction, temperature] = DistillationCurve(A_B(p, 1:numComponents-1), classes, index_n_eta_AB(p,:), pressure_distillation);
    y_AB(p,:) = interp1(volumeFraction, temperature, volFrac_interp, 'linear');
end
for s = 1:numel(volFrac_interp)
    V_1 = abs(mean((y_A(:,s) - c(s)).*(y_AB(:,s) - y_B(:,s))));
    sobol_idx(1,s) = V_1 / Vy(s);
end

% +++ First-order Sobol' indices for number of carbon atoms group +++ %
A_B = B;
y_AB = zeros(size(A_B,1), numel(volFrac_interp));
A_B(:, numComponents:2*numComponents-1) = A(:, numComponents:2*numComponents-1);
index_n_eta_AB = zeros(size(A_B,1), numComponents);
parfor p = 1:size(A_B,1)
    for q = 1:numComponents
        index_n_eta_AB(p,q) = findIndex_eta(classes{q}, A_B(p, numComponents+q-1), A_B(p, 2*numComponents+q-1));
    end
end
parfor p = 1:size(A_B,1)
    [volumeFraction, temperature] = DistillationCurve(A_B(p, 1:numComponents-1), classes, index_n_eta_AB(p,:), pressure_distillation);
    y_AB(p,:) = interp1(volumeFraction, temperature, volFrac_interp, 'linear');
end
for s = 1:numel(volFrac_interp)
    V_2 = abs(mean((y_A(:,s) - c(s)).*(y_AB(:,s) - y_B(:,s))));
    sobol_idx(2,s) = V_2 / Vy(s);
end

% +++ First-order Sobol' indices for branching indices group +++ %
A_B = B;
y_AB = zeros(size(A_B,1), numel(volFrac_interp));
A_B(:, 2*numComponents:3*numComponents-1) = A(:, 2*numComponents:3*numComponents-1);
index_n_eta_AB = zeros(size(A_B,1), numComponents);
parfor p = 1:size(A_B,1)
    for q = 1:numComponents
        index_n_eta_AB(p,q) = findIndex_eta(classes{q}, A_B(p, numComponents+q-1), A_B(p, 2*numComponents+q-1));
    end
end
parfor p = 1:size(A_B,1)
    [volumeFraction, temperature] = DistillationCurve(A_B(p, 1:numComponents-1), classes, index_n_eta_AB(p,:), pressure_distillation);
    y_AB(p,:) = interp1(volumeFraction, temperature, volFrac_interp, 'linear');
end
for s = 1:numel(volFrac_interp)
    V_3 = abs(mean((y_A(:,s) - c(s)).*(y_AB(:,s) - y_B(:,s))));
    sobol_idx(3,s) = V_3 / Vy(s);
end

end