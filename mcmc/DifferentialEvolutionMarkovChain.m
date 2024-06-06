function [x, p_x, x_nb, p_x_nb, AR, R_hat, t_convergence] = DifferentialEvolutionMarkovChain(Posterior_PDF, families, LowerBound_x, UpperBound_x, n_ranges, LowerBound_eta_B_star, UpperBound_eta_B_star, maxIterations, scaling_factor, noise_X, noise_nc, noise_eta, N_chains, outlier_method)

% DESCRIPTION:
% The DifferentialEvolutionMarkovChain function adopts the differential
% evolution Markov chain (DE-MC) algorithm proposed by ter Braak (2006) to
% explore and sample from the posterior probability density function (PDF).
% Note that the DifferentialEvolutionMarkovChain function is adapted from the DiffeRential
% Evolution Adaptive Metropolis (DREAM) toolbox developed by Vrugt et al. (2008, 2009, 2016).

% Auxiliary parameters:
% Nc: number of surrogate mixture components  [-]

% Inputs:
% 1) Posterior_PDF      : function handle representing the posterior PDF
% 2) families           : (1 x Nc) cell array, with the i-th cell containing a string
% that denotes the i-th hydrocarbon family under consideration
% 3) LowerBound_x       : (1 x Nc) array, with the i-th element representing the
% lower bound for the range the mole fraction of the i-th surrogate mixture
% component can vary within  [-]
% 4) UpperBound_x       : (1 x Nc) array, with the i-th element representing the
% upper bound for the range the mole fraction of the i-th surrogate mixture
% component can vary within  [-]
% 5) n_ranges           : (1 x Nc) cell array, with the i-th cell containing the range
% the number of carbon atoms of the i-th surrogate mixture component can vary within
% 6) LowerBound_eta_B_star : (1 x Nc) array, with the i-th element representing the
% lower bound for the range the normalized topochemical atom index for the i-th
% surrogate mixture component can vary within  [-]
% 7) UpperBound_eta_B_star : (1 x Nc) array, with the i-th element representing the
% upper bound for the range the normalized topochemical atom index for the i-th
% surrogate mixture component can vary within  [-]
% 8) maxIterations      : maximum number of iterations for each chain
% 9) scaling_factor     : scaling factor for the jump rate
% 10) noise_X           : noise parameter for the surrogate mixture
% molar fractions
% 11) noise_nc          : noise parameter for the surrogate mixture
% numbers of carbon atoms
% 12) noise_eta         : noise parameter for the surrogate mixture
% topochemical atom indices
% 13) N_chains          : number of chains
% 14) outlier_method    : string denoting the outlier chain detection
% method, namely, either 'iqr' (interquartile range) or 'mad' (median absolute deviation)

% Outputs:
% 1) x       : (t_convergence x 3*Nc-1 x N_chains) array containing the
% samples from the posterior PDF from each chain along the entire number of
% iterations, concerning the molar fractions, numbers of carbon atoms, and
% topochemical atom indices
% 2) p_x     : (t_convergence x N_chains) array containing the posterior density
% values from each chain along the entire number of iterations
% 3) x_nb    : (t_convergence/2 x 3*Nc-1 x N_chains) array containing the
% samples from the posterior PDF from each chain in the DE-MC algorithm discarding the burn-in period,
% concerning the molar fractions, numbers of carbon atoms, and topochemical atom indices
% 4) p_x_nb  : (t_convergence/2 x N_chains) array containing the posterior density
% values from each chain discarding the burn-in period
% 5) AR      : (t_convergence x N_chains) array containing the acceptance rate from each chain
% along the entire number of iterations
% 6) R_hat   : (t_convergence-0.1*floor(maxIterations)) x 3*Nc-1 array describing the R-hat statistic
% for each parameter to be inferred after 10% of the total number of iterations
% 7) t_convergence: number of iterations denoting the number of iterations the convergence of the DE-MC run
% is reached at according to the R-hat statistic ( = maxIterations in case
% convergence is not reached)

% --- References ---

% ter Braak, C.J.F., 2006. 
% A Markov chain Monte Carlo version of the genetic algorithm
% differential evolution: easy Bayesian computing for real parameter spaces.
% Stat. Comput. 16, 239-249.

% Vrugt, J.A., ter Braak, C.J.F., Clark, M.P., Hyman, J.M., Robinson, B.A., 2008. 
% Treatment of input uncertainty in hydrologic modeling: doing hydrology backward
% with Markov chain Monte Carlo simulation.
% Water Resour. Res. 44, W00B09.

% Vrugt, J.A., ter Braak, C.J.F., Diks, C.G.H., Higdon, D., Robinson, B.A., Hyman, J.M., 2009.
% Accelerating Markov chain Monte Carlo simulation by differential evolution with self-adaptive
% randomized subspace sampling.
% Int. J. Nonlinear Sci. Numer. Simul. 10 (3), 273-290.

% Vrugt, J.A., 2016. 
% Markov chain Monte Carlo simulation using the DREAM software package: Theory, concepts, and MATLAB implementation.
% Environ. Model. Softw. 75, 273-316.

% ------------------------------------------------------------------------%
% Contributors/Copyright
% 2024 Jacopo Liberatori, jacopo.liberatori@uniroma1.it
% 2024 Davide Cavalieri,  davide.cavalieri@uniroma1.it
% Department of Mechanical and Aerospace Engineering (DIMA)
% Sapienza University of Rome
% ------------------------------------------------------------------------%

fid = fopen('MCMC.txt', 'w'); % Open the file for writing

% Write the ASCII art header
fprintf(fid, '    ____                  _____ ___    ______                            \n');
fprintf(fid, '   / __ )____ ___  _____ / ___//   |  / ____/                            \n');
fprintf(fid, '  / __  / __ `/ / / / _ \\\\__ \\/ /| | / /_                             \n');
fprintf(fid, ' / /_/ / /_/ / /_/ /  __/__/ / ___ |/ __/                                \n');
fprintf(fid, '/_____/\\__,_/\\__, /\\___/____/_/  |_/_/                                \n');
fprintf(fid, '            /____/                                                       \n');
fprintf(fid, '\n'); % Add an extra newline for separation
fprintf(fid, '-----------------------------------------------------------------------  \n');
fprintf(fid, 'Contributors/Copyright                                                   \n');
fprintf(fid, '2024 Jacopo Liberatori, jacopo.liberatori@uniroma1.it                    \n');
fprintf(fid, '2024 Davide Cavalieri,  davide.cavalieri@uniroma1.it                     \n');
fprintf(fid, 'Department of Mechanical and Aerospace Engineering (DIMA)                \n');
fprintf(fid, 'Sapienza University of Rome                                              \n');
fprintf(fid, '-----------------------------------------------------------------------  \n');
fprintf(fid, '\n'); % Add an extra newline for separation

% Write header in the text file
fprintf(fid, '+++++ Differential evolution Markov chain (DE-MC) algorithm +++++\n');
fprintf(fid, 'Maximum number of iterations: %s\n', num2str(maxIterations));
fprintf(fid, 'Number of chains: %s\n', num2str(N_chains));

numComponents = numel(families);

gamma = scaling_factor*2.38/sqrt(2*(3*numComponents-1));                                            % Calculate default jump rate
x = nan(maxIterations,3*numComponents-1,N_chains); p_x = nan(maxIterations,N_chains); accept = nan(maxIterations,N_chains); AR = nan(maxIterations,N_chains); R_hat = 1e+18*ones((maxIterations-0.1*floor(maxIterations)), 3*numComponents-1);          % Preallocate memory for chains, density, acceptance rate, and R-hat statistic
% Create initial population by sampling using Latin hypercube sampling (LHS) and fulfilling constraints from the informative priors
X = 1e+18*ones(N_chains,3*numComponents-1);     % (N-1) mole fractions + N numbers of carbon atoms + N topochemical atom indices

p_X = nan(N_chains,1);
lhs = lhsdesign(N_chains, 3*numComponents-1); % generate initial LHS design
chain_init = 0; % number of chains being initialized
for i = 1:N_chains
    while (sum(X(i,1:numComponents-1)) > 1 - LowerBound_x(end)) || (sum(X(i,1:numComponents-1)) < 1 - UpperBound_x(end))
        for j = 1:numComponents-1
            X(i,j) = LowerBound_x(j) + (UpperBound_x(j) - LowerBound_x(j)) * lhs(i,j);
        end
        if (sum(X(i,1:numComponents-1)) <= 1 - LowerBound_x(end)) && (sum(X(i,1:numComponents-1)) >= 1 - UpperBound_x(end))
            chain_init = chain_init + 1;
        else
            lhs = lhsdesign(N_chains, 3*numComponents-1);
        end
    end

    for k = numComponents:2*numComponents-1
        X(i,k) = min(n_ranges{k-numComponents+1}) + floor((max(n_ranges{k-numComponents+1}) - min(n_ranges{k-numComponents+1})) * lhs(i,k));
    end

    for l = 2*numComponents:3*numComponents-1
        X(i,l) = LowerBound_eta_B_star(l-2*numComponents+1) + (UpperBound_eta_B_star(l-2*numComponents+1) - LowerBound_eta_B_star(l-2*numComponents+1)) * lhs(i,l);
        classl = Hydrocarbons(families{l-2*numComponents+1}, n_ranges{l-2*numComponents+1}); % Pick topochemical atom indices from database displaying the closest value to those sampled
        eta_B_star_norm_l_values = [classl.eta_B_star_norm];
        X(i,l) = eta_B_star_norm_l_values(findIndex_eta(classl, X(i,l)));
    end

    p_X(i,1) = Posterior_PDF(X(i,1:numComponents-1),X(i,numComponents:2*numComponents-1),X(i,2*numComponents:3*numComponents-1));     % Compute density initial population
end

x(1,1:3*numComponents-1,1:N_chains) = reshape(X',1,3*numComponents-1,N_chains); p_x(1,1:N_chains) = p_X';   % Store initial states and density

p_Xp = nan(N_chains, 1); % Initialize p_Xp array for storing proposed posterior densities

accept(1,:) = 1;                                                                 % First sample has been accepted
AR(1,:) = 100;

for i = 1:N_chains, R(i,1:N_chains-1) = setdiff(1:N_chains,i); end          % R-matrix: index of chains for differential evolution

% Initialize temporary variables
temp_accept = accept;

counter_convergence = 0;

for t = 2:maxIterations         % Dynamic part: evolution of N chains
    [~, draw] = sort(rand(N_chains-1,N_chains));                            % Permute [1,...,N_chains-1] N_chains times
    g = randsample([gamma 1], 1, 'true', [0.9 0.1]);                    % Select gamma: 90/10 mix [default 1]
    for i = 1:N_chains                                                      % Create proposals and accept/reject
        a = R(i,draw(1,i)); b = R(i,draw(2,i));                             % Extract a and b not equal to i
        Xp(i,1:numComponents-1) = X(i,1:numComponents-1) + g*(X(a,1:numComponents-1)-X(b,1:numComponents-1))...
            + noise_X*randn(1,numComponents-1);
        Xp(i,numComponents:2*numComponents-1) =  X(i,numComponents:2*numComponents-1) + g*(X(a,numComponents:2*numComponents-1)-X(b,numComponents:2*numComponents-1))...
            + noise_nc*randn(1,numComponents);
        Xp(i,numComponents:2*numComponents-1) = round(Xp(i,numComponents:2*numComponents-1));
        Xp(i,2*numComponents:3*numComponents-1) =  X(i,2*numComponents:3*numComponents-1) + g*(X(a,2*numComponents:3*numComponents-1)-X(b,2*numComponents:3*numComponents-1))...
            + noise_eta*randn(1,numComponents);
        for j = 2*numComponents:3*numComponents-1 % Pick topochemical atom indices from database displaying the closest value to those sampled
            classj = Hydrocarbons(families{j-2*numComponents+1}, n_ranges{j-2*numComponents+1});
            eta_B_star_norm_j_values = [classj.eta_B_star_norm];
            Xp(i,j) = eta_B_star_norm_j_values(findIndex_eta(classj, Xp(i,j)));
        end
    end

    parfor i =1:N_chains
        Xp_temp = Xp(i,:); % Create a temporary array for this iteration

        % Boundary handling for molar fractions by folding the parameter space
        for j = 1:numComponents-1
            if Xp_temp(j) < LowerBound_x(j)
                dist = LowerBound_x(j) - Xp_temp(j);
                Xp_temp(j) = max(UpperBound_x(j) - dist, LowerBound_x(j));
            end
            if Xp_temp(j) > UpperBound_x(j)
                dist = Xp_temp(j) - UpperBound_x(j);
                Xp_temp(j) = min(LowerBound_x(j) + dist, UpperBound_x(j));
            end
        end

        % Boundary handling for numbers of carbon atoms by folding the parameter space
        for j = numComponents:2*numComponents-1
            if Xp_temp(j) < min(n_ranges{j-numComponents+1})
                dist = min(n_ranges{j-numComponents+1}) - Xp_temp(j);
                Xp_temp(j) = max(max(n_ranges{j-numComponents+1}) - dist, min(n_ranges{j-numComponents+1}));
            end
            if Xp_temp(j) > max(n_ranges{j-numComponents+1})
                dist = Xp_temp(j) - max(n_ranges{j-numComponents+1});
                Xp_temp(j) = min(min(n_ranges{j-numComponents+1}) + dist, max(n_ranges{j-numComponents+1}));
            end
        end

        % Boundary handling for topochemical atom indices by folding the parameter space
        for j = 2*numComponents:3*numComponents-1
            if Xp_temp(j) < LowerBound_eta_B_star(j-2*numComponents+1)
                dist = LowerBound_eta_B_star(j-2*numComponents+1) - Xp_temp(j);
                Xp_temp(j) = max(UpperBound_eta_B_star(j-2*numComponents+1) - dist, LowerBound_eta_B_star(j-2*numComponents+1));
            end
            if Xp_temp(j) > UpperBound_eta_B_star(j-2*numComponents+1)
                dist = Xp_temp(j) - UpperBound_eta_B_star(j-2*numComponents+1);
                Xp_temp(j) = min(LowerBound_eta_B_star(j-2*numComponents+1) + dist, UpperBound_eta_B_star(j-2*numComponents+1));
            end
        end

        % Posterior evaluates to zero if sum(Xi) > 1 or at least one molar fraction is negative
        if sum(Xp_temp(1:numComponents-1)) > 1 - LowerBound_x(end) || sum(Xp_temp(1:numComponents-1)) < 1 - UpperBound_x(end) || any(Xp_temp < 0)
            p_Xp(i,1) = -1e+300;
        else
            p_Xp(i,1) = Posterior_PDF(Xp_temp(1:numComponents-1),Xp_temp(numComponents:2*numComponents-1),Xp_temp(2*numComponents:3*numComponents-1));      % Calculate density i-th proposal
        end

        % Assign the temporary array back to Xp
        Xp(i, :) = Xp_temp;

        p_acc = min(0,p_Xp(i,1)-p_X(i,1));                                  % Compute acceptance probability
            
        if p_acc > log(rand)                                                     % p_acc larger than U[0,1]?
            X(i,:) = Xp(i,:); p_X(i,1) = p_Xp(i,1);     % True: accept proposal
            temp_accept(t,i) = accept(t-1,i) + 1;                                                % How many samples accepted?
        else
            temp_accept(t,i) = accept(t-1,i);
        end
    end

    % Update 'accept' variable outside the loop
    accept = temp_accept;

    AR(t,:) = 100*(accept(t,:)/(t-1));                                            % Calculate acceptance rate

    x(t,1:3*numComponents-1,1:N_chains) = reshape(X',1,3*numComponents-1,N_chains); p_x(t,1:N_chains) = p_X';       % Append current X and density

    %%% ------- Check for outlier chains at 10%, 20%, 30%, 40%, 50%, 60%, 70%, and 80% of the total number of iterations ------- %%%
    if t == 0.1*floor(maxIterations) || t == 0.2*floor(maxIterations) || t == 0.3*floor(maxIterations) || t == 0.4*floor(maxIterations) || t == 0.5*floor(maxIterations) || t == 0.6*floor(maxIterations) || t == 0.7*floor(maxIterations) || t == 0.8*floor(maxIterations)
        omega(:) = mean(p_x(floor(t/2):t,1:N_chains));
        % Set "best" to be the index of the chain with maximum average
        [~,best] = max(omega);
        omega_sorted = sort(omega);
        %  Identify outlier chains, and replace their later samples with values from the "best" chain
        outlier_num = 0;
        outlier_indices = zeros(1, N_chains);
        if strcmp(outlier_method, 'mad')
            % Compute the median absolute deviation (mad) from the array containing the averages of the latest 50% samples by each chain
            medianPosterior = median(omega_sorted);
            mad = median(abs(omega_sorted-medianPosterior))*1.4826;
            for j = 1 : N_chains
                if omega(j) < medianPosterior - 3*mad
                    outlier_num = outlier_num + 1;
                    outlier_indices(j) = j;
                    X(j,1:3*numComponents-1) = X(best,1:3*numComponents-1);
                    p_X(j,1) = p_X(best,1);
                end
            end
        elseif strcmp(outlier_method, 'iqr')
            % Compute inter-quartile range (iqr) from the array containing the averages of the latest 50% samples by each chain
            q1 = quantile(omega_sorted, 0.25);
            q3 = quantile(omega_sorted, 0.75);
            iqr = q3 - q1;
            for j = 1 : N_chains
                if omega(j) < q1 - iqr
                    outlier_num = outlier_num + 1;
                    outlier_indices(j) = j;
                    X(j,1:3*numComponents-1) = X(best,1:3*numComponents-1);
                    p_X(j,1) = p_X(best,1);
                end
            end
        end
        outlier_indices = outlier_indices(outlier_indices ~= 0);
        if outlier_num == 1
            fprintf(fid, ['Chain # ', num2str(outlier_indices), ' was detected as an outlier chain after ', num2str(t), ' iterations. Its latest samples were replaced by samples from chain # ', num2str(best), ', characterized by the largest mean of the posterior densities computed over the second half of samples so far.\n']);
        elseif outlier_num > 1
            fprintf(fid, ['Chains # ', num2str(outlier_indices), ' were detected as outlier chains after ', num2str(t), ' iterations. Their latest samples were replaced by samples from chain # ', num2str(best), ', characterized by the largest mean of the posterior densities computed over the second half of samples so far.\n']);
        end
    end

    %%% ------- Check convergence through R-statistic at least after 10% of the total number of iterations (if convergence has not been reached yet) ------- %%%
    if t >= 0.1*floor(maxIterations)

        W = zeros(3*numComponents-1,1);
        B = zeros(3*numComponents-1,1);
        sigma = zeros(3*numComponents-1,1);

        for j = 1:3*numComponents-1
            x_j_r = x(floor(t/2):t, j, :);
            x_j_r_bar = zeros(N_chains,1);
            for r = 1:N_chains
                x_j_r_bar(r) = 2/(t-2)*sum(x_j_r(:, r));
            end
            inner_sum_W = 0;
            for r = 1:N_chains
                for i = 1:size(x_j_r,1)
                    inner_sum_W = inner_sum_W + (x_j_r(i, r) - x_j_r_bar(r))^2;
                end
            end
            W(j) = (2 / (N_chains * (t - 2))) * inner_sum_W;
            x_j_bar_bar = 1/N_chains*sum(x_j_r_bar);
            inner_sum_B = 0;
            for r = 1:N_chains
                inner_sum_B = inner_sum_B + (x_j_r_bar(r) - x_j_bar_bar)^2;
            end
            B(j) = t/(2*(N_chains-1))*inner_sum_B;
            sigma(j) = (t-2)/t*W(j)+2/t*B(j);
            R_hat(t-0.1*floor(maxIterations)+1,j) = sqrt((N_chains+1)/N_chains*sigma(j)/W(j)-(t-2)/N_chains/t);
        end

        if all(R_hat(t-0.1*floor(maxIterations)+1,:) <= 1.2) % ---- Convergence is reached ---- %
            t_convergence = t;
            R_hat = R_hat(1:t_convergence-0.1*floor(maxIterations)+1,:);
            AR = AR(1:t_convergence,:);
            fprintf(fid, ['R-statistic below the critical threshold of 1.2 for all parameters of the target distribution: convergence reached after ', num2str(t_convergence), ' iterations within each chain. The second half of these samples will be used for posterior analysis, while the first half is discarded as burn-in.\n']);
            counter_convergence = 1;
            break
        end
        
    end
end         % End dynamic part

if counter_convergence == 0     % ---- Convergence not reached ---- %
    t_convergence = maxIterations;
    fprintf(fid, ['Convergence not reached! Posterior analysis will be performed considering the second half from the entirety of ', num2str(maxIterations), ' MCMC samples, with the first half being discarded as burn-in.\n']);
end

t_burnin = floor(t_convergence/2); % number of iterations to be discarded as burn-in

%%% ------- Discard burn-in samples ------- %%%
for l = 1:N_chains
    x_nb((l-1)*t_burnin+1:l*(t_convergence-t_burnin), :) = x(t_burnin+1:t_convergence, :, l);
    p_x_nb((l-1)*t_burnin+1:l*(t_convergence-t_burnin), :, 1) = p_x(t_burnin+1:t_convergence, l);
end

fclose(fid); % Close the file