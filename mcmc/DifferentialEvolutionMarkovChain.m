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

function [x, p_x, x_nb, p_x_nb, AR, R_hat, t_convergence] = DifferentialEvolutionMarkovChain(Posterior_PDF, classes, LowerBound_x, UpperBound_x, n_ranges, LowerBound_eta_B_star, UpperBound_eta_B_star, maxIterations, t_burnin, scaling_factor, noise_X, noise_nc, noise_eta, N_chains, outlier_method, R_hat_threshold)

% DESCRIPTION:
% The DifferentialEvolutionMarkovChain function adopts the differential
% evolution Markov chain (DE-MC) algorithm proposed by ter Braak (2006) to
% explore and sample from the posterior probability density function (PDF).
% Note that the DifferentialEvolutionMarkovChain function is adapted from the DiffeRential
% Evolution Adaptive Metropolis (DREAM) toolbox developed by Vrugt et al. (2008, 2009, 2016).

% Auxiliary parameters:
% Nc: number of surrogate mixture components  [-]

% Inputs:
% 1) Posterior_PDF        : function handle representing the posterior PDF
% 2) classes              : (1 x Nc) cell array, with the i-th cell containing a structure
% array for the hydrocarbon family and range of number of carbon atoms of the i-th surrogate component
% 3) LowerBound_x         : (1 x Nc) array, with the i-th element representing the
% lower bound for the range the mole fraction of the i-th surrogate mixture
% component can vary within  [-]
% 4) UpperBound_x         : (1 x Nc) array, with the i-th element representing the
% upper bound for the range the mole fraction of the i-th surrogate mixture
% component can vary within  [-]
% 5) n_ranges             : (1 x Nc) cell array, with the i-th cell containing the range
% the number of carbon atoms of the i-th surrogate mixture component can vary within
% 6) LowerBound_eta_B_star: (1 x Nc) array, with the i-th element representing the
% lower bound for the range the normalized topochemical atom index for the i-th
% surrogate mixture component can vary within  [-]
% 7) UpperBound_eta_B_star: (1 x Nc) array, with the i-th element representing the
% upper bound for the range the normalized topochemical atom index for the i-th
% surrogate mixture component can vary within  [-]
% 8) maxIterations        : maximum number of iterations for each chain
% 9) t_burnin             : number of iterations to be discarded as burn-in
% 10) scaling_factor      : scaling factor for the jump rate
% 11) noise_X             : noise parameter for the surrogate mixture
% molar fractions
% 12) noise_nc            : noise parameter for the surrogate mixture
% numbers of carbon atoms
% 13) noise_eta           : noise parameter for the surrogate mixture
% topochemical atom indices
% 14) N_chains            : number of chains
% 15) outlier_method      : string denoting the outlier chain detection
% method, namely, either 'iqr' (interquartile range) or 'mad' (median absolute deviation)
% 16) R_hat_threshold     : threshold value for the R-hat statistic indicating convergence of the DE-MC algorithm

% Outputs:
% 1) x       : (t_convergence x 3*Nc-1 x N_chains) array containing the
% samples from the posterior PDF from each chain along the entire number of
% iterations, concerning the molar fractions, numbers of carbon atoms, and
% topochemical atom indices
% 2) p_x     : (t_convergence x N_chains) array containing the posterior density
% values from each chain along the entire number of iterations
% 3) x_nb    : (((t_convergence - t_burnin) x N_chains) x 3*Nc-1) array containing the
% samples from the posterior PDF from each chain in the DE-MC algorithm discarding the burn-in period,
% concerning the molar fractions, numbers of carbon atoms, and topochemical atom indices
% 4) p_x_nb  : ((t_convergence - t_burnin) x 1) array containing the posterior density
% values from each chain discarding the burn-in period
% 5) AR      : (t_convergence x N_chains) array containing the acceptance rate from each chain
% along the entire number of iterations
% 6) R_hat   : floor((t_convergence-t_burnin)/2) x 3*Nc-1 array describing the R-hat statistic
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

fid = fopen('MCMC.txt', 'w'); % Open the file for writing

% Write the ASCII art header
fprintf(fid, '    ____                  _____ ___    ______                               \n');
fprintf(fid, '   / __ )____ ___  _____ / ___//   |  / ____/                               \n');
fprintf(fid, '  / __  / __ `/ / / / _ \\\\__ \\/ /| | / /_                                \n');
fprintf(fid, ' / /_/ / /_/ / /_/ /  __/__/ / ___ |/ __/                                   \n');
fprintf(fid, '/_____/\\__,_/\\__, /\\___/____/_/  |_/_/                                   \n');
fprintf(fid, '            /____/                                                          \n');
fprintf(fid, '\n'); % Add an extra newline for separation
fprintf(fid, '--------------------------------------------------------------------------  \n');
fprintf(fid, 'Contributors / Copyright Notice                                             \n');
fprintf(fid, '© 2025 Jacopo Liberatori — jacopo.liberatori@centralesupelec.fr             \n');
fprintf(fid, 'Postdoctoral Researcher @ CentraleSupélec, Laboratoire EM2C (CNRS)          \n');
fprintf(fid, '\n'); % Add an extra newline for separation
fprintf(fid, '© 2025 Davide Cavalieri — davide.cavalieri@uniroma1.it                      \n');
fprintf(fid, 'Postdoctoral Researcher @ Sapienza University of Rome,                      \n');
fprintf(fid, 'Department of Mechanical and Aerospace Engineering (DIMA)                   \n');
fprintf(fid, '\n'); % Add an extra newline for separation
fprintf(fid, '© 2025 Matteo Blandino — matteo.blandino@uniroma1.it                        \n');
fprintf(fid, 'Ph.D. Student @ Sapienza University of Rome,                                \n');
fprintf(fid, 'Department of Mechanical and Aerospace Engineering (DIMA)                   \n');
fprintf(fid, '\n'); % Add an extra newline for separation
fprintf(fid, 'Reference:                                                                  \n');
fprintf(fid, 'J. Liberatori, D. Cavalieri, M. Blandino, M. Valorani, and P.P. Ciottoli.   \n');
fprintf(fid, 'BayeSAF: Emulation and Design of Sustainable Alternative Fuels via Bayesian \n');
fprintf(fid, 'Inference and Descriptors-Based Machine Learning. Fuel (under review), 2025.\n');
fprintf(fid, 'Available at: https://dx.doi.org/10.2139/ssrn.5049145.                      \n');
fprintf(fid, '-----------------------------------------------------------------------     \n');
fprintf(fid, '\n'); % Add an extra newline for separation

% Write header in the text file
fprintf(fid, '+++++ Differential evolution Markov chain (DE-MC) algorithm +++++\n');
fprintf(fid, 'Maximum number of iterations: %s\n', num2str(maxIterations));
fprintf(fid, 'Number of chains: %s\n', num2str(N_chains));

rng('default'); % to ensure reproducibility

numComponents = numel(classes);

gamma = scaling_factor*2.38/sqrt(2*(3*numComponents-1));                                            % Calculate default jump rate
x = nan(maxIterations,3*numComponents-1,N_chains); eta_B_star_chain = nan(maxIterations,numComponents,N_chains); p_x = nan(maxIterations,N_chains); accept = nan(maxIterations,N_chains); AR = nan(maxIterations,N_chains); R_hat = 1e+18*ones(floor((maxIterations-2*t_burnin)/2), 3*numComponents-1);          % Preallocate memory for chains, topochemical atom indices, density, acceptance rate, and R-hat statistic
% Create initial population by sampling using Latin hypercube sampling (LHS) and fulfilling constraints from the informative priors
X = 1e+18*ones(N_chains,3*numComponents-1);     % (N-1) mole fractions + N numbers of carbon atoms + N topochemical atom indices

p_X = nan(N_chains,1);
% Constrained Latin hypercube sampling (LHS) to initialize DE-MC chains
validSamples = [];

while size(validSamples, 1) < N_chains
    
    % Oversample
    lhs = lhsdesign(10*N_chains, 3*numComponents - 1);
    
    % === 1. Extract LHS parts ===
    lhs_x = lhs(:, 1:numComponents-1); % mole fractions
    lhs_rest = lhs(:, numComponents:end); % number of carbon atoms and normalized topochemical atom indices

    % === 2. Rescale mole fractions ===
    scaled_x = lhs_x .* (UpperBound_x(1:end-1) - LowerBound_x(1:end-1)) + LowerBound_x(1:end-1);

    % === 3. Apply constraint on sum ===
    xSum = sum(scaled_x, 2);
    constraint = (xSum < 1 - LowerBound_x(end)) & (xSum > 1 - UpperBound_x(end));

    % === 4. Retain only valid samples ===
    valid_x = scaled_x(constraint, :);
    valid_rest = lhs_rest(constraint, :);

    % === 5. Recombine and store ===
    validSamples = [validSamples; [valid_x, valid_rest]];
end

% Trim to the number of DE-MC chains
lhs_constrained = validSamples(1:N_chains, :);

t = 1; % starting iteration

for i = 1:N_chains

    for j = 1:numComponents-1
        X(i,j) = lhs_constrained(i,j);
    end

    for k = numComponents:2*numComponents-1
        nC_k = min(n_ranges{k-numComponents+1}) + round((max(n_ranges{k-numComponents+1}) - min(n_ranges{k-numComponents+1})) * lhs_constrained(i,k));
        [~, idx_neighb] = min(abs(n_ranges{k-numComponents+1} - nC_k));
        X(i,k) = n_ranges{k-numComponents+1}(idx_neighb);
    end

    for l = 2*numComponents:3*numComponents-1
        eta_l = LowerBound_eta_B_star(l-2*numComponents+1) + (UpperBound_eta_B_star(l-2*numComponents+1) - LowerBound_eta_B_star(l-2*numComponents+1)) * lhs_constrained(i,l);
        classl = classes{l-2*numComponents+1}; % Pick topochemical atom indices from database displaying the closest value to those sampled
        eta_B_star_norm_list = [classl.eta_B_star_norm];
        eta_B_star_norm_values(i,l-2*numComponents+1) = eta_B_star_norm_list(findIndex_eta(classl, X(i,l-numComponents), eta_l));
        X(i,l) = find(sort([classl([classl.nC] == X(i,l-numComponents)).eta_B_star_norm])==eta_B_star_norm_values(i,l-2*numComponents+1))/numel([classl([classl.nC] == X(i,l-numComponents)).eta_B_star_norm]);
    end

    p_X(i,t) = Posterior_PDF(X(i,1:numComponents-1),X(i,numComponents:2*numComponents-1),eta_B_star_norm_values(i,:));     % Compute density initial population

end

x(t,1:3*numComponents-1,1:N_chains) = reshape(X',1,3*numComponents-1,N_chains); eta_B_star_chain(t,1:numComponents,1:N_chains) = eta_B_star_norm_values'; p_x(t,1:N_chains) = p_X';   % Store initial states, topochemical atom indices and density

p_Xp = nan(N_chains, 1); % Initialize p_Xp array for storing proposed posterior densities

accept(t,:) = 1;                                                                 % First sample has been accepted
AR(t,:) = 100;

for i = 1:N_chains, R(i,1:N_chains-1) = setdiff(1:N_chains,i); end          % R-matrix: index of chains for differential evolution

% Initialize temporary variables
temp_accept = accept;

convergence = false; % switch

counter_check = 1;

% Create a waitbar
f = waitbar(0,'1','Name','Differential evolution Markov Chain (DE-MC)','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(f,'canceling',0);
hTitle = findall(f, 'Type', 'Axes');
hText = findall(f, 'Type', 'Text');
set(get(hTitle, 'Title'), 'FontName', 'Times New Roman', 'FontSize', 16);
set(hText, 'FontName', 'Times New Roman', 'FontSize', 16);
set(f, 'Color', [0.8, 0.8, 0.8]);
hPatch = findall(f, 'Type', 'Patch');
set(hPatch, 'FaceColor', [0, 0.5, 0]);
set(hText, 'Color', [0.3, 0.3, 0.3]);
figure('visible', 'off');

while ~convergence && t <= maxIterations         % Dynamic part: evolution of N chains

    t = t + 1;

    if getappdata(f,'canceling')
        break
    end
    waitbar(t/maxIterations,f,sprintf(strcat('Number of iterations: ', num2str(t))))

    [~, draw] = sort(rand(N_chains-1,N_chains));                            % Permute [1,...,N_chains-1] N_chains times
    g = randsample([gamma 1], 1, 'true', [0.9 0.1]);                        % Select gamma: 90/10 mix [default 1]
    Xp = zeros(N_chains,3*numComponents-1);
    for i = 1:N_chains                                                      % Create proposals and accept/reject
        a = R(i,draw(1,i)); b = R(i,draw(2,i));                             % Extract a and b not equal to i
        Xp(i,1:numComponents-1) = X(i,1:numComponents-1) + g*(X(a,1:numComponents-1)-X(b,1:numComponents-1))...
            + noise_X*randn(1,numComponents-1);
        for j = numComponents:2*numComponents-1
            Xp(i,j) =  X(i,j) + g*(X(a,j)-X(b,j)) + noise_nc*randn;
            Xp(i,numComponents:2*numComponents-1) = round(Xp(i,numComponents:2*numComponents-1));
        end
        for k = 2*numComponents:3*numComponents-1
            Xp(i,k) =  X(i,k) + g*(X(a,k)-X(b,k)) + noise_eta*randn;
            Xp(i,2*numComponents:3*numComponents-1) = Xp(i,2*numComponents:3*numComponents-1);
        end
        
        % Boundary handling for molar fractions by folding the parameter space
        for j = 1:numComponents-1
            if Xp(i,j) < LowerBound_x(j)
                dist = LowerBound_x(j) - Xp(i,j);
                Xp(i,j) = max(UpperBound_x(j) - dist, LowerBound_x(j));
            end
            if Xp(i,j) > UpperBound_x(j)
                dist = Xp(i,j) - UpperBound_x(j);
                Xp(i,j) = min(LowerBound_x(j) + dist, UpperBound_x(j));
            end
        end

        % Boundary handling for numbers of carbon atoms by folding the parameter space
        for j = numComponents:2*numComponents-1
            if Xp(i,j) < min(n_ranges{j-numComponents+1})
                dist = min(n_ranges{j-numComponents+1}) - Xp(i,j);
                Xp(i,j) = max(max(n_ranges{j-numComponents+1}) - dist, min(n_ranges{j-numComponents+1}));
            end
            if Xp(i,j) > max(n_ranges{j-numComponents+1})
                dist = Xp(i,j) - max(n_ranges{j-numComponents+1});
                Xp(i,j) = min(min(n_ranges{j-numComponents+1}) + dist, max(n_ranges{j-numComponents+1}));
            end
        end

        % Boundary handling for topochemical atom indices by folding the parameter space and pick topochemical atom indices from database
        for j = 2*numComponents:3*numComponents-1
            if Xp(i,j) < 0
                dist = 0 - Xp(i,j);
                Xp(i,j) = max(1 - abs(dist), 0);
            end
            if Xp(i,j) > 1
                dist = Xp(i,j) - 1;
                Xp(i,j) = min(1 + dist, 1);
            end
            classj = classes{j-2*numComponents+1};
            eta_B_star_norm_list = sort([classj([classj.nC] == Xp(i,j-numComponents)).eta_B_star_norm]);
            eta_index = max(round(numel(eta_B_star_norm_list)*Xp(i,j)),1);
            eta_B_star_norm_values(i,j-2*numComponents+1) = eta_B_star_norm_list(eta_index);
        end

    end

    u = rand(1,N_chains);

    parfor i =1:N_chains
        Xp_temp = Xp(i,:); % Create a temporary array for this iteration

        % Posterior evaluates to -Inf if sum(Xi) > 1 or at least one molar fraction is negative
        if sum(Xp_temp(1:numComponents-1)) > 1 - LowerBound_x(end) || sum(Xp_temp(1:numComponents-1)) < 1 - UpperBound_x(end)
            p_Xp(i,1) = -Inf;
        else
            p_Xp(i,1) = Posterior_PDF(Xp_temp(1:numComponents-1),Xp_temp(numComponents:2*numComponents-1),eta_B_star_norm_values(i,:));      % Calculate density i-th proposal
        end

        % Assign the temporary array back to Xp
        Xp(i, :) = Xp_temp;

        % Metropolis-Hastings acceptance criterion
        if p_Xp(i,1) > p_X(i,1)
            X(i,:) = Xp(i,:); p_X(i,1) = p_Xp(i,1);     % True: accept proposal
            temp_accept(t,i) = accept(t-1,i) + 1;       % How many samples accepted?
        else
            if u(i) < exp(p_Xp(i,1) - p_X(i,1))
                X(i,:) = Xp(i,:); p_X(i,1) = p_Xp(i,1);     % True: accept proposal
                temp_accept(t,i) = accept(t-1,i) + 1;       % How many samples accepted?
            else
                temp_accept(t,i) = accept(t-1,i);
            end
        end
    end

    % Update 'accept' variable outside the loop
    accept = temp_accept;

    AR(t,:) = 100*(accept(t,:)/(t-1));                                            % Calculate acceptance rate

    x(t,1:3*numComponents-1,1:N_chains) = reshape(X',1,3*numComponents-1,N_chains); eta_B_star_chain(t,1:numComponents,1:N_chains) = eta_B_star_norm_values'; p_x(t,1:N_chains) = p_X';       % Append current X, topochemical atom indices and density

    %%% ------- Check for outlier chains from 2*t_burnin on ------- %%%
    if t >= 2*t_burnin && mod(t,t_burnin) == 0
        omega(:) = mean(p_x(t_burnin:t,1:N_chains));
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
                if omega(j) < q1 - 1.5*iqr
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

    %%% ------- Check convergence through R-statistic from 2*t_burnin on ------- %%%
    if t >= 2*t_burnin && mod(t,2) == 0

        W = zeros(3*numComponents-1,1);
        B = zeros(3*numComponents-1,1);
        sigma = zeros(3*numComponents-1,1);

        for j = 1:3*numComponents-1
            x_j_r = x(t_burnin:t, j, :);
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
            R_hat(counter_check,j) = sqrt((N_chains+1)/N_chains*sigma(j)/W(j)-(t-2)/N_chains/t);
        end

        if all(R_hat(counter_check,:) <= R_hat_threshold) % ---- Convergence is reached ---- %
            t_convergence = t;
            R_hat = R_hat(1:counter_check,:);
            AR = AR(1:t_convergence,:);
            fprintf(fid, ['R-statistic below the critical threshold of ', num2str(R_hat_threshold), ' for all parameters of the target distribution: convergence reached after ', num2str(t_convergence), ' iterations within each chain. These samples will be used for posterior analysis after discarding the burn-in period.\n']);
            convergence = true;
        end

        counter_check = counter_check + 1;
        
    end

end         % End dynamic part

% Close the waitbar
delete(f);

if ~convergence     % ---- Convergence not reached ---- %
    t_convergence = maxIterations;
    fprintf(fid, ['Convergence not reached! Posterior analysis will be performed considering the entirety of ', num2str(maxIterations), ' MCMC samples after discarding the burn-in period.\n']);
end

x = x(1:t_convergence,:,:);
eta_B_star_chain = eta_B_star_chain(1:t_convergence,:,:);
x(:, 2*numComponents:3*numComponents-1, :) = eta_B_star_chain;
p_x = p_x(1:t_convergence,:);

%%% ------- Discard burn-in samples ------- %%%
for l = 1:N_chains
    x_nb((l-1)*(t_convergence-t_burnin)+1:l*(t_convergence-t_burnin), :) = x(t_burnin+1:end, :, l);
    p_x_nb((l-1)*(t_convergence-t_burnin)+1:l*(t_convergence-t_burnin), :, 1) = p_x(t_burnin+1:end, l);
end

fclose(fid); % Close the file
