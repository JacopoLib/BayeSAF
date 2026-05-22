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
% The DifferentialEvolutionMarkovChain function adopts the differential
% evolution Markov chain (DE-MC) algorithm proposed by ter Braak (2006) to
% explore and sample from the posterior probability density function (PDF).
% Note that the DifferentialEvolutionMarkovChain function is adapted from the DiffeRential
% Evolution Adaptive Metropolis (DREAM) toolbox developed by Vrugt et al. (2008, 2009, 2016).

% Auxiliary parameters:
% Nc: number of surrogate mixture components  [-]

% Inputs:
% 1) families             : (1 x Nc) cell array, with the i-th cell containing a string
% that denotes the i-th hydrocarbon family under consideration
% 2) Posterior_PDF        : function handle representing the posterior PDF
% 3) classes              : (1 x Nc) cell array, with the i-th cell containing a structure
% array for the hydrocarbon family and range of number of carbon atoms of the i-th surrogate component
% 4) LowerBound_x         : (1 x Nc) array, with the i-th element representing the
% lower bound for the range the mole fraction of the i-th surrogate mixture
% component can vary within  [-]
% 5) UpperBound_x         : (1 x Nc) array, with the i-th element representing the
% upper bound for the range the mole fraction of the i-th surrogate mixture
% component can vary within  [-]
% 6) n_ranges             : (1 x Nc) cell array, with the i-th cell containing the range
% the number of carbon atoms of the i-th surrogate mixture component can vary within
% 7) maxIterations        : maximum number of iterations for each chain
% 8) t_burnin             : number of iterations to be discarded as burn-in
% chains after the latest burn-in period is detected (if at least one
% outlier chain is detected, all samples collected so far are discarded as
% burn-in)
% 9) t_check              : number of iterations to start R-hat convergence check from
% 10) scaling_factor_X    : scaling factor for the jump rate concerning molar fractions
% 11) scaling_factor_nc   : scaling factor for the jump rate concerning numbers of carbon atoms
% 12) scaling_factor_eta  : scaling factor for the jump rate concerning topochemical atom indices
% 13) noise_X             : noise parameter for the surrogate mixture
% molar fractions
% 14) noise_nc            : noise parameter for the surrogate mixture
% numbers of carbon atoms
% 15) noise_eta           : noise parameter for the surrogate mixture
% topochemical atom indices
% 16) N_chains            : number of chains
% 17) outlier_method      : string denoting the outlier chain detection
% method, namely, either 'iqr' (interquartile range) or 'mad' (median absolute deviation)
% 18) R_hat_threshold     : threshold value for the R-hat statistic indicating convergence of the DE-MC algorithm
% 19) PT_switch           : use parallel tempering ('True') or not ('False')
% 20) beta_min            : minimum inverse temperature for parallel tempering
% 21) T_ladder            : temperature ladder functional form (1. 'linear', 2. 'geometric')
% 22) swap_freq           : swap frequency for parallel tempering
% 23) n_pairs             : maximum number of chain pairs used to propose the new sample
% 24) p_gibbs             : probability of a Gibbs move
% 25) p_snooker           : probability of performing a snooker jump rather than a DE move

% Outputs:
% 1) x            : (t_convergence x 3*Nc-1 x N_chains) array containing the
% samples from the posterior PDF from each chain along the entire number of
% iterations, concerning the molar fractions, numbers of carbon atoms, and
% topochemical atom indices
% 2) p_x          : (t_convergence x N_chains) array containing the posterior density
% values from each chain along the entire number of iterations
% 3) x_nb         : (((t_convergence - t_burnin) x N_chains) x 3*Nc-1) array containing the
% samples from the posterior PDF from each chain in the DE-MC algorithm discarding the burn-in period,
% concerning the molar fractions, numbers of carbon atoms, and topochemical atom indices
% 4) p_x_nb       : ((t_convergence - t_burnin) x 1) array containing the posterior density
% values from each chain discarding the burn-in period
% 5) AR           : (t_convergence x N_chains) array containing the acceptance rate from each chain
% along the entire number of iterations
% 6) R_hat        : floor((t_convergence - t_burnin)/2) x 3*Nc-1 array describing
% the R-hat statistic for each parameter to be inferred
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
% ------------------------------------------------------------------------

function [x, p_x, x_nb, p_x_nb, AR, R_hat, t_convergence] = DifferentialEvolutionMarkovChain_singlePT_delayed(families, Posterior_PDF, Posterior_PDF_cheap, classes, LowerBound_x, UpperBound_x, n_ranges, LowerBound_eta_B_star, UpperBound_eta_B_star, maxIterations, t_burnin, t_check, scaling_factor_X, scaling_factor_nc, scaling_factor_eta, noise_X, noise_nc, noise_eta, N_chains, outlier_method, R_hat_threshold, PT_switch, beta_min, T_ladder, swap_freq, n_pairs, p_gibbs, p_snooker, varargin)
% Optional benchmark arguments (passed as varargin):
%   varargin{1} = X_init_bench  : (N_chains x 3*numComponents-1) normalised init states
%   varargin{2} = p_X_init_bench: (N_chains x 1) log-posterior at init states

% === Deterministic Seeding Pattern (Official Implementation) ===
% Ensures reproducibility across runs and parallel workers.

% Choose a fixed base seed for this DE-MC experiment
baseSeed = 1234;

% Set the global RNG state deterministically
rng(baseSeed, 'twister');

% If a parallel pool exists, assign deterministic but distinct streams
pool = gcp('nocreate');
if ~isempty(pool)
    % Create an independent stream for each worker using the baseSeed offset
    spmd
        workerSeed = baseSeed + spmdIndex;
        rng(workerSeed, 'twister');
    end
end
% ================================================================

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
fprintf(fid, '© 2026 Jacopo Liberatori — jacopo.liberatori@centralesupelec.fr             \n');
fprintf(fid, 'Postdoctoral Researcher @ Laboratoire EM2C, CentraleSupélec (CNRS)          \n');
fprintf(fid, '\n'); % Add an extra newline for separation
fprintf(fid, '© 2026 Davide Cavalieri — davide.cavalieri@uniroma1.it                      \n');
fprintf(fid, 'Postdoctoral Researcher @ Sapienza University of Rome,                      \n');
fprintf(fid, 'Department of Mechanical and Aerospace Engineering (DIMA)                   \n');
fprintf(fid, '\n'); % Add an extra newline for separation
fprintf(fid, '© 2026 Matteo Blandino, Ph.D.                                               \n');
fprintf(fid, '\n'); % Add an extra newline for separation
fprintf(fid, 'Reference:                                                                  \n');
fprintf(fid, 'J. Liberatori, D. Cavalieri, M. Blandino, M. Valorani, and P.P. Ciottoli.   \n');
fprintf(fid, 'BayeSAF: Emulation and Design of Sustainable Alternative Fuels via Bayesian \n');
fprintf(fid, 'Inference and Descriptors-Based Machine Learning. Fuel 419, 138835 (2026).  \n');
fprintf(fid, 'Available at: https://doi.org/10.1016/j.fuel.2026.138835.                    \n');
fprintf(fid, '-----------------------------------------------------------------------     \n');
fprintf(fid, '\n'); % Add an extra newline for separation

% Write header in the text file
fprintf(fid, '+++++ Differential evolution Markov chain (DE-MC) algorithm +++++\n');
fprintf(fid, 'Maximum number of iterations: %s\n', num2str(maxIterations));
fprintf(fid, 'Number of chains: %s\n', num2str(N_chains));

numComponents = numel(classes);

counter_swap = 0;

x = nan(maxIterations,3*numComponents-1,N_chains); molFrac_chain = nan(maxIterations,numComponents-1,N_chains); nC_chain = nan(maxIterations,numComponents,N_chains); eta_B_star_chain = nan(maxIterations,numComponents,N_chains); p_x = nan(maxIterations,N_chains); accept = nan(maxIterations,N_chains); AR = nan(maxIterations,N_chains); R_hat = 1e+18*ones(floor((maxIterations - t_burnin)/2), 3*numComponents-1); ESS = 1e+18*ones(floor((maxIterations - t_burnin)/2), 3*numComponents-1);          % Preallocate memory for chains, topochemical atom indices, density, acceptance rate, and R-hat statistic
% Create initial population by sampling using Latin hypercube sampling (LHS) and fulfilling constraints from the informative priors
X = 1e+18*ones(N_chains,3*numComponents-1);     % (N-1) mole fractions + N numbers of carbon atoms + N topochemical atom indices

p_X = nan(N_chains,1);

% === Initialization ===
t = 1; % starting iteration — defined here so it is available in both branches

% --- Benchmark shortcut: use externally provided initial states if supplied ---
use_bench_init = numel(varargin) >= 2 && ~isempty(varargin{1});
if use_bench_init
    X_init_bench  = varargin{1};   % N_chains x n_params (normalised)
    p_X_init_bench = varargin{2}(:); % N_chains x 1

    X = X_init_bench;
    p_X = p_X_init_bench;

    % Decode physical parameters for chain metadata storage
    molFrac_values        = zeros(N_chains, numComponents - 1);
    nC_values             = zeros(N_chains, numComponents);
    eta_B_star_norm_values = zeros(N_chains, numComponents);

    for i_b = 1:N_chains
        % Molar fractions
        for j_b = 1:numComponents-1
            molFrac_values(i_b,j_b) = X_init_bench(i_b,j_b) * ...
                (UpperBound_x(j_b) - LowerBound_x(j_b)) + LowerBound_x(j_b);
        end
        % Carbon numbers
        for k_b = 1:numComponents
            n_range_k = n_ranges{k_b};
            kk = numel(n_range_k);
            u = X_init_bench(i_b, numComponents - 1 + k_b);
            idx_nc = min(floor(u * kk) + 1, kk);
            nC_values(i_b, k_b) = n_range_k(idx_nc);
        end
        % Topochemical indices
        for l_b = 1:numComponents
            classl = classes{l_b};
            eta_list  = [classl.eta_B_star];
            eta_sorted = sort(eta_list);
            eta_norm  = [classl.eta_B_star_norm];
            u = X_init_bench(i_b, 2*numComponents - 1 + l_b);
            idx_eta = round(1 + u * (numel(eta_sorted) - 1));
            idx_eta = max(1, min(idx_eta, numel(eta_sorted)));
            eta_target = eta_sorted(idx_eta);
            nC_match   = find([classl.nC] == nC_values(i_b, l_b));
            [~, idx_nb] = min(abs(eta_list(nC_match) - eta_target));
            eta_B_star_norm_values(i_b, l_b) = eta_norm(nC_match(idx_nb));
        end
    end
else
% --- Standard random initialisation (unchanged) ---
Nneeded = 10 * N_chains;
validSamples = [];

K = numComponents; % total components
K_rest = 3 * K - 1 - (K - 1); % same shape as your lhs_rest in original code

batchSize = max(1000, 5 * Nneeded);

while size(validSamples, 1) < Nneeded

    % --- 1) Sample mole fractions uniformly on the simplex (Dirichlet(1))
    Y = -log(rand(batchSize, K));     % Exponential(1)
    DirSamples = Y ./ sum(Y, 2);      % Normalize → each row sums to 1

    % --- 2) Check per-component lower/upper bounds
    lb = LowerBound_x(:)';
    ub = UpperBound_x(:)';

    ok_mask = all((DirSamples >= lb) & (DirSamples <= ub), 2);
    val_dir = DirSamples(ok_mask, :);

    if isempty(val_dir)
        continue
    end

    % --- 3) Compute normalized mole fractions (for first K-1 components)
    % Normalized to [0,1] within their individual bounds
    val_norm = (val_dir(:, 1:K-1) - lb(1:K-1)) ./ (ub(1:K-1) - lb(1:K-1));

    % --- 4) Generate random values for remaining parameters in [0,1]
    n_valid = size(val_norm, 1);
    lhs_rest_cont = rand(n_valid, K_rest); % uniform in [0,1]

    % --- 5) Combine into one array
    validSamples = [validSamples; [val_norm, lhs_rest_cont]];

end

p_X_lhs = nan(size(validSamples,1),1);
X_lhs = 1e+18*ones(size(validSamples,1),3*numComponents-1);
molFrac_values_lhs = 1e+18*ones(size(validSamples,1),numComponents-1);
eta_B_star_norm_values_lhs = 1e+18*ones(size(validSamples,1),numComponents);
nC_values_lhs = 1e+18*ones(size(validSamples,1),numComponents);

for i = 1:size(validSamples,1)

    for j = 1:numComponents-1
        X_lhs(i,j) = validSamples(i,j);
        molFrac_values_lhs(i,j) = validSamples(i,j).* (UpperBound_x(j) - LowerBound_x(j)) + LowerBound_x(j);
    end

    for k = numComponents:2*numComponents-1

        X_lhs(i,k) = validSamples(i,k);

        % === 1. Retrieve allowed carbon numbers ===
        n_range_k = n_ranges{k - numComponents + 1};
        kk = numel(n_range_k);

        % === 2. Extract normalized parameter (in [0,1]) ===
        u = validSamples(i,k);

        % === 3. Map normalized u to discrete index using equal-width bins ===
        idx_nc = min(floor(u * kk) + 1, kk);   % idx_nc ∈ [1, k]

        % === 4. Assign discrete carbon number ===
        nC_values_lhs(i, k - numComponents + 1) = n_range_k(idx_nc);

    end

    for l = 2*numComponents:3*numComponents-1
        classl = classes{l-2*numComponents+1}; % Pick topochemical atom indices from database displaying the closest value to those sampled
        eta_B_star_list = [classl.eta_B_star];
        eta_B_star_list_sorted = sort(eta_B_star_list);
        eta_B_star_norm_list = [classl.eta_B_star_norm];
        idx_eta_l = round(1 + validSamples(i,l) * (numel(eta_B_star_list_sorted) - 1));
        eta_l = eta_B_star_list_sorted(idx_eta_l);

        nC_indices = find([classl.nC] == nC_values_lhs(i,l-2*numComponents+1));
        eta_B_star_list_l = eta_B_star_list(nC_indices);
        eta_B_star_norm_list_l = eta_B_star_norm_list(nC_indices);

        [~, idx_neighb] = min(abs(eta_B_star_list_l - eta_l));

        eta_B_star_norm_values_lhs(i,l-2*numComponents+1) = eta_B_star_norm_list_l(idx_neighb);
        X_lhs(i,l) = validSamples(i,l);
    end

end

p_X_lhs = Posterior_PDF(molFrac_values_lhs,nC_values_lhs,eta_B_star_norm_values_lhs);     % Compute density initial population

phys_samples = [molFrac_values_lhs, nC_values_lhs, eta_B_star_norm_values_lhs];

N = size(phys_samples,1);
selectedIdx = zeros(N_chains,1);

% Step 1: randomly pick the first row
selectedIdx(1) = randi(N);

% Step 2: iteratively select rows maximizing min distance to current set
for k = 2:N_chains
    remaining = setdiff(1:N, selectedIdx(1:k-1));
    maxMinDist = -Inf;
    bestIdx = remaining(1);
    for i = remaining
        dists = pdist2(phys_samples(i,:), phys_samples(selectedIdx(1:k-1),:));
        minDist = min(dists); % min distance to already selected rows
        if minDist > maxMinDist
            maxMinDist = minDist;
            bestIdx = i;
        end
    end
    selectedIdx(k) = bestIdx;
end

p_X = p_X_lhs(selectedIdx,:);
X = X_lhs(selectedIdx,:);
molFrac_values = molFrac_values_lhs(selectedIdx,:);
nC_values = nC_values_lhs(selectedIdx,:);
eta_B_star_norm_values = eta_B_star_norm_values_lhs(selectedIdx,:);

perm = randperm(N_chains);

X = X(perm, :);
p_X = p_X(perm);
molFrac_values = molFrac_values(perm, :);
nC_values = nC_values(perm, :);
eta_B_star_norm_values = eta_B_star_norm_values(perm, :);

end % if use_bench_init / else

x(t,1:3*numComponents-1,1:N_chains) = reshape(X',1,3*numComponents-1,N_chains);
molFrac_chain(t,1:numComponents-1,1:N_chains) = molFrac_values';
nC_chain(t,1:numComponents,1:N_chains) = nC_values';
eta_B_star_chain(t,1:numComponents,1:N_chains) = eta_B_star_norm_values';
p_x(t,1:N_chains) = p_X';   % Store initial states, topochemical atom indices and density

p_Xp = nan(N_chains, 1); % Initialize p_Xp array for storing proposed posterior densities

accept(t,:) = 1;                                                                 % First sample has been accepted
AR(t,:) = 100;

for i = 1:N_chains, R(i,1:N_chains-1) = setdiff(1:N_chains,i); end          % R-matrix: index of chains for differential evolution

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

% Use parallel tempering or not
if strcmp(PT_switch, 'True')
    % Define inverse temperatures (β) for parallel tempering
    if strcmp(T_ladder, 'linear')
        beta = linspace(1, beta_min, N_chains);   % coldest first, hottest last
    elseif strcmp(T_ladder, 'geometric')
        beta = beta_min.^((0:N_chains-1)/(N_chains-1));
        beta = sort(beta,'descend');              % coldest first, hottest last
    end
else
    beta = ones(1,N_chains);
    swap_freq = Inf;
end

% Preallocate swap counters (for adjacent pairs)
nPairs = floor(N_chains/2);
swap_accept_count = zeros(floor(maxIterations/swap_freq), nPairs);
swap_attempt_count = zeros(floor(maxIterations/swap_freq), nPairs);

x_archive = zeros(maxIterations,3*numComponents-1,N_chains);

while ~convergence && t <= maxIterations         % Dynamic part: evolution of N chains

    [~, idx_best] = maxk(p_X,1);

    param_orig = [molFrac_values, nC_values, eta_B_star_norm_values];

    t = t + 1;

    if t > t_burnin
        beta = ones(1,N_chains);
        swap_freq = Inf;
    end

    if getappdata(f,'canceling')
        break
    end
    waitbar(t/maxIterations,f,sprintf(strcat('Number of iterations: ', num2str(t))))

    lambda = unifrnd(-0.1, 0.1, N_chains, 1); % draw N_chains lambda values
    [~, draw] = sort(rand(N_chains-1,N_chains));
    Xp = X;

    J = zeros(1,N_chains);

    S1f = ones(1,N_chains);
    log_pi1_x = zeros(1,N_chains);
    log_pi1_y = zeros(1,N_chains);

    for i = 1:N_chains                                                      % Create proposals and accept/reject

        r_gibbs = rand;

        Xp_temp = Xp(i,:);
        molFrac_values_temp = molFrac_values(i,:);
        nC_values_temp = nC_values(i,:);
        eta_B_star_norm_values_temp = eta_B_star_norm_values(i,:);

        if r_gibbs > p_gibbs
            block = 'group';
            group_idx = randi(3); %randi(3); % pick one subset of compositional parameters
            if group_idx == 1
                block_group = 'fractions';
            elseif group_idx == 2
                block_group = 'nC';
            elseif group_idx == 3
                block_group = 'eta';
            end
        else
            % block = 'all';
            block = 'group';
            block_group = 'fractions';
            % block_group = 'none';
        end

        % Pick random number to determine whether to perform a DE move or a snooker jump
        r_move = rand;

        if t > 10*(3*numComponents-1) && t > t_burnin

            % --------------------- SNOOKER JUMP --------------------- %
            idx_chain = randperm(N_chains);
            idx_chain(idx_chain==i) = [];
            r1 = idx_chain(1); r2 = idx_chain(2); r3 = idx_chain(3);

            rowsToRemove = all(x_archive == 0, [2 3]);
            x_archive_filt = x_archive(~rowsToRemove, :, :);
            idx_sample = round(unifrnd(1,size(x_archive_filt,1),3,1));

            % Define reference point xR and direction vectors
            xR = x_archive_filt(idx_sample(1), :, r1);
            z_snooker = xR - X(i,:);
            v = x_archive_filt(idx_sample(2),:,r2) - x_archive_filt(idx_sample(3),:,r3);
            alpha = (z_snooker * v') / (z_snooker * z_snooker' + 1e-8);
            z_proj = alpha * z_snooker;
            g_snooker = unifrnd(1.2, 2.2);

        end

        % --------------------- DE MOVE --------------------- %
        D = randsample(1:n_pairs, 1, 'true');
        a = R(i,draw(1:D,i));
        b = R(i,draw(D+1:2*D,i));

        switch block

            case 'group'


                % Calculate default jump rates
                gamma_X = scaling_factor_X*2.38/sqrt(2*D*(numComponents-1));
                gamma_nc = scaling_factor_nc*2.38/sqrt(2*D*numComponents);
                gamma_eta = scaling_factor_eta*2.38/sqrt(2*D*numComponents);
                g_X = randsample([gamma_X 1], 1, 'true', [0.9 0.1]); % Select gamma: 90/10 mix [default 1]
                mask = rand(size(gamma_nc)) < 0.9;   % true with 90% prob
                g_nc = gamma_nc;                     % start with original values
                g_nc(~mask) = 1;                     % replace 10% of them with 1
                g_eta = randsample([gamma_eta 1], 1, 'true', [0.9 0.1]); % Select gamma: 90/10 mix [default 1]

                switch block_group

                    case 'fractions'

                        if r_move < p_snooker && t > 10*(3*numComponents-1) && t > t_burnin

                            Xp_temp(1:numComponents-1) = X(i,1:numComponents-1) + g_snooker * z_proj(1:numComponents-1) + noise_X * randn(1,numComponents-1);

                            % --- Compute Snooker Jacobian correction ---
                            d = numComponents - 1;
                            zp = xR(1:numComponents-1) - Xp_temp(1:numComponents-1);
                            epsJ = 1e-12;
                            J(i) = (d - 1) * ( log(norm(z_snooker(1:numComponents-1)) + epsJ) - ...
                                log(norm(zp) + epsJ) );

                        else

                            for j = 1:numComponents-1
                                Xp_temp(j) = X(i,j) ...
                                    + (1 - lambda(i)) * (g_X * sum(X(a,j)-X(b,j), 1)) ...
                                    + noise_X*randn;
                            end

                        end

                    case 'nC'

                        for j = numComponents:2*numComponents-1
                            Xp_temp(j) = X(i,j) ...
                                + (1 - lambda(i)) * (g_nc * sum(X(a,j)-X(b,j), 1)) ...
                                + noise_nc*randn;
                        end

                    case 'eta'

                        for j = 2*numComponents:3*numComponents-1
                            Xp_temp(j) = X(i,j) ...
                                + (1 - lambda(i)) * (g_eta * sum(X(a,j)-X(b,j), 1)) ...
                                + noise_eta*randn;
                        end

                    case 'none'

                end

        end

        for j = 1:numComponents-1
            if Xp_temp(j) < 0
                dist = 0 - Xp_temp(j);
                Xp_temp(j) = max(1 - abs(dist), 0);
            end
            if Xp_temp(j) > 1
                dist = Xp_temp(j) - 1;
                Xp_temp(j) = min(0 + dist, 1);
            end
            molFrac_values_temp(j) = Xp_temp(j).* (UpperBound_x(j) - LowerBound_x(j)) + LowerBound_x(j);
        end

        if sum(molFrac_values_temp) <= 1 - LowerBound_x(end) && sum(molFrac_values_temp) >= 1 - UpperBound_x(end) && r_gibbs <= p_gibbs

            % === CONFIGURATION ===
            if t <= t_burnin %- 1e+18 %rand < 0.9 - 1e+18
                numToUpdate = 1;
                maxSamples  = Inf*100;
                maxComb     = Inf*10*numToUpdate;
            else
                numToUpdate = numComponents;
                maxSamples  = Inf*100;
                maxComb     = 10*numToUpdate + (N_chains-1);
            end
            N_batches   = 1;
            NBH         = Inf;

            % === Select random distinct components to update ===
            components = randperm(numComponents, numToUpdate);

            % === Containers for pair lists and metadata ===
            pair_lists_fwd      = cell(1, numToUpdate);
            n_ranges_sel        = cell(1, numToUpdate);
            classes_sel         = cell(1, numToUpdate);
            eta_sorted_all      = cell(1, numToUpdate);
            eta_list_all        = cell(1, numToUpdate);
            cur_idx_each        = zeros(1, numToUpdate);

            % -------------------------------------------------------------------------
            % === Build candidate (nC, eta) pairs for each selected component ===
            for cc = 1:numToUpdate
                j_comp  = components(cc);
                jn = j_comp + (numComponents - 1);
                je = jn + numComponents;

                % === Retrieve per-component info ===
                n_ranges_sel{cc} = n_ranges{j_comp};
                classes_sel{cc}  = classes{je - 2*numComponents + 1};

                classj        = classes_sel{cc};
                eta_norm_list = [classj.eta_B_star_norm];

                % === Identify neighborhood for nC ===
                n_range_j = n_ranges_sel{cc};
                nCurr = nC_values_temp(j_comp);
                [~, idx_curr] = min(abs(n_range_j - nCurr));
                idx_low  = max(1, idx_curr - NBH);
                idx_high = min(numel(n_range_j), idx_curr + NBH);
                idx_neighb = idx_low:idx_high;

                % === Step 1: include all nC in neighborhood ===
                nC_candidates = n_range_j(idx_neighb);

                % === Step 2: randomly sample eta for each nC ===
                pair_list = [];
                for zz = 1:numel(nC_candidates)
                    nC_j = nC_candidates(zz);
                    nC_idx = find([classj.nC] == nC_j);
                    eta_norm_j = eta_norm_list(nC_idx);
                    if isempty(eta_norm_j), continue; end
                    n_eta = numel(eta_norm_j);
                    take_eta = min(maxSamples, n_eta);
                    take_eta = max(1, take_eta);
                    eta_idcs = randperm(n_eta, take_eta);
                    eta_samples = eta_norm_j(eta_idcs);
                    pair_list = [pair_list; [repmat(nC_j, take_eta, 1), eta_samples(:)]];
                end

                % === Step 3: ensure current (nC, eta) is included ===
                cur_pair = [nC_values_temp(j_comp), eta_B_star_norm_values_temp(j_comp)];
                [tf, loc] = ismember(cur_pair, pair_list, 'rows');
                if ~tf
                    pair_list = [pair_list; cur_pair];
                    loc = size(pair_list, 1);
                end
                cur_idx_each(cc) = loc;

                if t > t_burnin %- 1e+18
                    % === Step 3b: include discrete states from all other chains ===
                    other_idx = setdiff(1:N_chains, i); % all chains except current
                    donor_pairs = [nC_values(other_idx, j_comp), eta_B_star_norm_values(other_idx, j_comp)];
                    donor_pairs = unique(donor_pairs, 'rows'); % remove duplicates
                    pair_list = unique([pair_list; donor_pairs], 'rows'); % append & deduplicate
                end

                % === Step 4: store component-specific data ===
                pair_lists_fwd{cc} = pair_list;
                eta_list_all{cc}   = [classj.eta_B_star];
                eta_sorted_all{cc} = sort([classj.eta_B_star]);
            end

            % -------------------------------------------------------------------------
            % === Cartesian product (with sub-sampling if needed) ===
            numCombEach = cellfun(@(x) size(x,1), pair_lists_fwd);
            rangeCells   = arrayfun(@(n) 1:n, numCombEach, 'UniformOutput', false);
            idx_grids    = cell(1, numToUpdate);

            nCombFull = prod(numCombEach);
            targetPerDim = max(1, ceil(maxComb^(1/numToUpdate)));

            if nCombFull > maxComb
                % --- Stratified sub-sample per dimension ---
                idx_sub = cell(1, numToUpdate);
                for cc = 1:numToUpdate
                    P   = pair_lists_fwd{cc};
                    n   = size(P,1);
                    inc = cur_idx_each(cc);
                    targ_cc = max(1, ceil(maxComb^(1/numToUpdate)));

                    if n <= targ_cc
                        idx_sub{cc} = 1:n;
                        continue
                    end

                    % Group by nC
                    nC_vals = P(:,1);
                    [nC_unique, ~, grp] = unique(nC_vals);
                    G = numel(nC_unique);
                    targ_cc = min(n, max(targ_cc, G));

                    rows_by_grp = cell(G,1);
                    for g = 1:G, rows_by_grp{g} = find(grp == g); end

                    selected = zeros(0,1);
                    for g = 1:G
                        rg = rows_by_grp{g};
                        selected(end+1,1) = rg(randi(numel(rg))); %#ok<AGROW>
                    end

                    % Ensure current included
                    if ~ismember(inc, selected)
                        g_inc = grp(inc);
                        sel_mask = (grp(selected) == g_inc);
                        if any(sel_mask)
                            selected(find(sel_mask,1,'first')) = inc;
                        else
                            selected(end+1,1) = inc; %#ok<AGROW>
                        end
                    end

                    % Fill remaining quota round-robin
                    rem_by_grp = cell(G,1);
                    for g = 1:G
                        rem = setdiff(rows_by_grp{g}, selected, 'stable');
                        rem_by_grp{g} = rem(randperm(numel(rem)));
                    end
                    g_ptr = ones(G,1);
                    while numel(selected) < targ_cc
                        for g = 1:G
                            if numel(selected) >= targ_cc, break; end
                            if g_ptr(g) <= numel(rem_by_grp{g})
                                selected(end+1,1) = rem_by_grp{g}(g_ptr(g)); %#ok<AGROW>
                                g_ptr(g) = g_ptr(g) + 1;
                            end
                        end
                        if all(arrayfun(@(g) g_ptr(g) > numel(rem_by_grp{g}), 1:G)), break; end
                    end
                    idx_sub{cc} = selected(:).';
                end
                [idx_grids{:}] = ndgrid(idx_sub{:});
            else
                [idx_grids{:}] = ndgrid(rangeCells{:});
            end

            idx_pairs = cellfun(@(x) x(:), idx_grids, 'UniformOutput', false);
            idx_pairs = [idx_pairs{:}];
            nComb = size(idx_pairs, 1);

            % -------------------------------------------------------------------------
            % === Preallocate posterior input matrices ===
            molFrac_all     = repmat(molFrac_values_temp, nComb, 1);
            nC_values_all   = repmat(nC_values_temp, nComb, 1);
            eta_values_all  = repmat(eta_B_star_norm_values_temp, nComb, 1);

            % === Fill candidate combinations ===
            nC_all_mat  = zeros(nComb, numToUpdate);
            eta_all_mat = zeros(nComb, numToUpdate);
            for cc = 1:numToUpdate
                P = pair_lists_fwd{cc};
                nC_all_mat(:,cc)  = P(idx_pairs(:,cc), 1);
                eta_all_mat(:,cc) = P(idx_pairs(:,cc), 2);
            end
            nC_values_all(:, components)  = nC_all_mat;
            eta_values_all(:, components) = eta_all_mat;

            % -------------------------------------------------------------------------
            % === Identify current joint state and append if missing ===
            if numToUpdate == 1
                cur_idx_in_grid = cur_idx_each(1);
            else
                [tf_row, cur_idx_in_grid] = ismember(cur_idx_each, idx_pairs, 'rows');
                if ~tf_row
                    cur_idx_in_grid = nComb + 1;
                    idx_pairs = [idx_pairs; cur_idx_each];
                    nComb = nComb + 1;
                    molFrac_all(end+1,:) = molFrac_values_temp;
                    nC_values_all(end+1,:) = nC_values_temp;
                    eta_values_all(end+1,:) = eta_B_star_norm_values_temp;
                end
            end
            cur_idx_in_grid = max(1, min(cur_idx_in_grid, size(molFrac_all,1)));

            if t <= t_burnin %- 1e+18

                % --- Random sub-sampling that always keeps current state ---
            if nComb > maxComb
                idx_keep = randperm(nComb, maxComb - 1);

                % Guarantee inclusion of the current state
                if ~ismember(cur_idx_in_grid, idx_keep)
                    idx_keep(end+1) = cur_idx_in_grid;
                end

                % Defensive clamp — remove any accidental > nComb indices
                idx_keep = unique(idx_keep(idx_keep <= size(molFrac_all,1)));
            else
                idx_keep = 1:nComb;
            end

            else


            % -------------------------------------------------------------------------
            % === Append all other chains’ current discrete states (multi-chain guarantee) ===
            other_idx = setdiff(1:N_chains, i);
            for kk = other_idx
                molFrac_all(end+1,:)   = molFrac_values_temp;
                nC_values_all(end+1,:) = nC_values(kk,:);
                eta_values_all(end+1,:) = eta_B_star_norm_values(kk,:);
            end

            % Update counts after appending
            nComb = size(molFrac_all,1);
            idx_all_states = (nComb - numel(other_idx)) : nComb; % indices of appended donor states

            % -------------------------------------------------------------------------
            % === Random sub-sampling that always keeps current and donor states ===
            if nComb > maxComb
                idx_keep = randperm(nComb, maxComb - (1 + numel(idx_all_states)));
                % Guarantee inclusion of current + donors
                keep_set = unique([idx_keep, cur_idx_in_grid, idx_all_states]);
                keep_set = keep_set(keep_set <= nComb);
                idx_keep = keep_set;
            else
                idx_keep = 1:nComb;
            end

            end

            % -------------------------------------------------------------------------
            % === Subset candidates safely ===
            mola = molFrac_all(idx_keep,:);
            nca  = nC_values_all(idx_keep,:);
            etaa = eta_values_all(idx_keep,:);
            logp1_fwd = Posterior_PDF_cheap(mola,nca,etaa);
            nComb = numel(idx_keep);

            % -------------------------------------------------------------------------
            % === Stable normalization and sampling ===
            maxlogp1_fwd = max(logp1_fwd);
            w1_fwd = exp(logp1_fwd - maxlogp1_fwd);
            S1f(i) = sum(w1_fwd);
            p1_fwd_norm = w1_fwd / S1f(i);
            idx_sel = randsample(1:nComb, 1, true, p1_fwd_norm);
            log_pi1_y(i) = logp1_fwd(idx_sel);

            % -------------------------------------------------------------------------
            % === Extract selected (nC, eta) for all updated components ===
            for cc = 1:numToUpdate
                j_comp = components(cc);
                nC_sel  = nca(idx_sel, j_comp);
                eta_sel = etaa(idx_sel, j_comp);
                nC_values_temp(j_comp)              = nC_sel;
                eta_B_star_norm_values_temp(j_comp) = eta_sel;

                % normalized X updates
                n_range_j = n_ranges_sel{cc};
                jn = j_comp + (numComponents - 1);
                Xp_temp(jn) = (nC_sel - min(n_range_j)) / (max(n_range_j) - min(n_range_j));

                classj = classes_sel{cc};
                eta_list_sorted = eta_sorted_all{cc};
                eta_list_j = eta_list_all{cc}([classj.nC] == nC_sel);
                eta_unnorm = min(eta_list_j) + eta_sel * (max(eta_list_j) - min(eta_list_j));
                [~, idx_neighb] = min(abs(eta_list_sorted - eta_unnorm));
                je = jn + numComponents;
                Xp_temp(je) = (idx_neighb - 1) / (numel(eta_list_sorted) - 1);
            end

            log_pi1_x(i) = Posterior_PDF_cheap(molFrac_values(i,:), nC_values(i,:), eta_B_star_norm_values(i,:));

        elseif strcmp(block_group, 'nC') || strcmp(block_group, 'eta') % t <= t_burnin || r_gibbs > p_gibbs

            % Boundary handling for numbers of carbon atoms by folding the parameter space
            for j = numComponents:2*numComponents-1

                % if Xp_temp(j) < 0
                %     dist = 0 - Xp_temp(j);
                %     Xp_temp(j) = max(1 - abs(dist), 0);
                % end
                % if Xp_temp(j) > 1
                %     dist = Xp_temp(j) - 1;
                %     Xp_temp(j) = min(0 + dist, 1);
                % end

                % xj = Xp_temp(j);
                % xj = abs(mod(xj, 2));
                % xj(xj > 1) = 2 - xj(xj > 1);
                % Xp_temp(j) = xj;

                % --- Reflect normalized parameter into [0,1] ---
                xj = 1 - abs(1 - mod(Xp_temp(j), 2));  % periodic triangle-wave reflection
                Xp_temp(j) = xj;

                % === 1. Retrieve allowed carbon numbers ===
                n_range_j = n_ranges{j - numComponents + 1};
                kk = numel(n_range_j);

                % === 2. Extract normalized parameter (in [0,1]) ===
                u = Xp_temp(j);

                % === 3. Map normalized u to discrete index using equal-width bins ===
                idx_nc = min(floor(u * kk) + 1, kk);   % idx_nc ∈ [1, kk]

                % === 4. Assign discrete carbon number ===
                nC_values_temp(j - numComponents + 1) = n_range_j(idx_nc);

            end

            % Boundary handling for topochemical atom indices by folding the parameter space and pick topochemical atom indices from database
            for j = 2*numComponents:3*numComponents-1

                % if Xp_temp(j) < 0
                %     dist = 0 - Xp_temp(j);
                %     Xp_temp(j) = max(1 - abs(dist), 0);
                % end
                % if Xp_temp(j) > 1
                %     dist = Xp_temp(j) - 1;
                %     Xp_temp(j) = min(0 + dist, 1);
                % end
                % xj = Xp_temp(j);
                % xj = abs(mod(xj, 2));
                % xj(xj > 1) = 2 - xj(xj > 1);
                % Xp_temp(j) = xj;

                % --- Reflect normalized parameter into [0,1] ---
                xj = 1 - abs(1 - mod(Xp_temp(j), 2));  % same folding for normalized η
                Xp_temp(j) = xj;

                classj = classes{j-2*numComponents+1}; % Pick topochemical atom indices from database displaying the closest value to those sampled
                eta_B_star_list = [classj.eta_B_star];
                eta_B_star_list_sorted = sort(eta_B_star_list);
                eta_B_star_norm_list = [classj.eta_B_star_norm];
                idx_eta_j = round(1 + Xp_temp(j) * (numel(eta_B_star_list_sorted) - 1));
                eta_j = eta_B_star_list_sorted(idx_eta_j);

                nC_indices = find([classj.nC] == nC_values_temp(j-2*numComponents+1));
                eta_B_star_list_j = eta_B_star_list(nC_indices);
                eta_B_star_norm_list_j = eta_B_star_norm_list(nC_indices);

                [~, idx_neighb] = min(abs(eta_B_star_list_j - eta_j));

                eta_B_star_norm_values_temp(j-2*numComponents+1) = eta_B_star_norm_list_j(idx_neighb);

            end

        end

        Xp(i,:) = Xp_temp;
        molFrac_values(i,:) = molFrac_values_temp;
        nC_values(i,:) = nC_values_temp;
        eta_B_star_norm_values(i,:) = eta_B_star_norm_values_temp;

    end

    param_p = [molFrac_values, nC_values, eta_B_star_norm_values];

    % -------- Delayed-Acceptance: Stage 1 (cheap, vectorized) --------
    % Uses your precomputed cheap logs:
    %   log_pi1_x(i) = Posterior_PDF_cheap(current_i)
    %   log_pi1_y(i) = Posterior_PDF_cheap(proposed_i)
    % J(i) is snooker Jacobian if applicable, else 0.

    % if r_gibbs <= p_gibbs

    log_alpha1_vec = beta(:) .* (log_pi1_y(:) - log_pi1_x(:)) + J(:);
    u1 = log(rand(N_chains,1));
    pass1 = u1 < log_alpha1_vec;   % chains that pass Stage 1
    fail1 = ~pass1;

    % default: carry forward previous accept counter for all chains
    accept(t,:) = accept(t-1,:);

    if any(pass1)
        % -------- Delayed-Acceptance: Stage 2 (full, batched only on pass1) --------
        idx_pass = find(pass1);

        % Build the batched proposed states for those chains (use your *_temp + Xp)
        mola_pass = molFrac_values(idx_pass, :);
        nca_pass  = nC_values(idx_pass, :);
        etaa_pass = eta_B_star_norm_values(idx_pass, :);

        % One vectorized full-posterior call for Stage-2 candidates
        p_y_full_pass = Posterior_PDF(mola_pass, nca_pass, etaa_pass);  % column vector

        if any(p_y_full_pass == -Inf)
            aa = 2;
        end

        % DA correction: log α2 = β * ( [ℓ(y)-ℓ(x)] - [ℓc(y)-ℓc(x)] )
        delta_full  = p_y_full_pass - p_X(idx_pass);                         % full diff
        delta_cheap = (log_pi1_y(idx_pass) - log_pi1_x(idx_pass));           % cheap diff
        log_alpha2  = beta(idx_pass)' .* (delta_full(:) - delta_cheap(:)); % no J here

        u2 = log(rand(numel(idx_pass),1));
        accept2 = u2 < log_alpha2;

        if any(accept2)
            idx_acc = idx_pass(accept2);

            % Commit accepted proposals: normalized, physical, and full posterior
            X(idx_acc, :) = Xp(idx_acc, :);
            p_X(idx_acc, 1) = p_y_full_pass(accept2);

            % Update accept counters only for finally accepted chains
            accept(t, idx_acc) = accept(t-1, idx_acc) + 1;
        end
    end

    AR(t,:) = 100*(accept(t,:)/(t-1));                                            % Calculate acceptance rate

    % ----- Replica Exchange -----

    % Attempt swaps between adjacent temperature levels every 'swap_freq' iterations
    swap_idx = floor(t / swap_freq); % which swap interval we're in

    if mod(t, swap_freq) == 0
        if mod(swap_idx, 2) == 0
            pair_start = 1; % (1,2), (3,4), ...
        else
            pair_start = 2; % (2,3), (4,5), ...
        end
        for k = pair_start:2:(N_chains-1)
            i = k; j = k + 1;
            Delta = (beta(i) - beta(j)) * (p_X(j) - p_X(i));
            accept_swap = (Delta >= 0) || (log(rand) < Delta);
            swap_attempt_count(swap_idx, ceil(k/2)) = ...
                swap_attempt_count(swap_idx, ceil(k/2)) + 1;
            if accept_swap
                X([i j], :) = X([j i], :);
                p_X([i j]) = p_X([j i]);
                swap_accept_count(swap_idx, ceil(k/2)) = ...
                    swap_accept_count(swap_idx, ceil(k/2)) + 1;
            end
        end
    end

    if t == round(0.25*t_burnin) || t == round(0.5*t_burnin) || t == round(0.75*t_burnin) || t == t_burnin

        %%% ------- Check for outlier chains within the initial burn-in period ------- %%%

        omega(:) = mean(p_x(round(t/2):t,1:N_chains));
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
            fprintf(fid, ['Chain # ', num2str(outlier_indices), ' was detected as an outlier chain after ', num2str(t), ' iterations during the initial burn-in period. Its latest samples were replaced by samples from chain # ', num2str(best), ', characterized by the largest mean of the posterior densities.\n']);
        elseif outlier_num > 1
            fprintf(fid, ['Chains # ', num2str(outlier_indices), ' were detected as outlier chains after ', num2str(t), ' iterations during the initial burn-in period. Their latest samples were replaced by samples from chain # ', num2str(best), ', characterized by the largest mean of the posterior densities.\n']);
        end

    end

    x(t,1:3*numComponents-1,1:N_chains) = reshape(X',1,3*numComponents-1,N_chains);

    for i = 1:N_chains

        % +++ Molar fractions +++
        for j = 1:numComponents-1
            molFrac_values(i,j) = X(i,j).* (UpperBound_x(j) - LowerBound_x(j)) + LowerBound_x(j);
        end

        % +++ Numbers of carbon atoms +++
        for j = numComponents:2*numComponents-1

            % === 1. Retrieve allowed carbon numbers ===
            n_range_j = n_ranges{j - numComponents + 1};
            k = numel(n_range_j);

            % === 2. Extract normalized parameter (in [0,1]) ===
            u = X(i, j);

            % === 3. Map normalized u to discrete index using equal-width bins ===
            idx_nc = min(floor(u * k) + 1, k);   % idx_nc ∈ [1, k]

            % === 4. Assign discrete carbon number ===
            nC_values(i, j - numComponents + 1) = n_range_j(idx_nc);

        end

        % +++ Topochemical atom indices +++
        for j = 2*numComponents:3*numComponents-1
            classj = classes{j-2*numComponents+1}; % Pick topochemical atom indices from database displaying the closest value to those sampled
            eta_B_star_list = [classj.eta_B_star];
            eta_B_star_list_sorted = sort(eta_B_star_list);
            eta_B_star_norm_list = [classj.eta_B_star_norm];
            idx_eta_j = round(1 + X(i,j) * (numel(eta_B_star_list_sorted) - 1));
            eta_j = eta_B_star_list_sorted(max(1, min(idx_eta_j, numel(eta_B_star_list_sorted))));

            nC_indices = find([classj.nC] == nC_values(i,j-2*numComponents+1));
            eta_B_star_list_j = eta_B_star_list(nC_indices);
            eta_B_star_norm_list_j = eta_B_star_norm_list(nC_indices);

            [~, idx_neighb] = min(abs(eta_B_star_list_j - eta_j));

            eta_B_star_norm_values(i,j-2*numComponents+1) = eta_B_star_norm_list_j(idx_neighb);
        end

    end

    molFrac_chain(t,1:numComponents-1,1:N_chains) = molFrac_values';
    nC_chain(t,1:numComponents,1:N_chains) = nC_values';
    eta_B_star_chain(t,1:numComponents,1:N_chains) = eta_B_star_norm_values';

    p_x(t,1:N_chains) = p_X';       % Append current X, topochemical atom indices and density

    if mod(t,10) == 0
        x_archive(t,1:3*numComponents-1,1:N_chains) = x(t,:,:);
    end

    %%% ------- Check convergence through R-statistic after t_check iterations ------- %%%

    if t > t_check && mod(t,2) == 0

        T = t - (t_burnin + round((t - t_burnin)/2)); % number of samples in each chain
        idx_start = t_burnin + round((t - t_burnin)/2) + 1;

        % Preallocation
        nParams = 3*numComponents - 1;
        W = zeros(nParams,1);
        B = zeros(nParams,1);
        sigma = zeros(nParams,1);

        % Build x_rhat directly
        x_rhat = cat(2, molFrac_chain, nC_chain, eta_B_star_chain);
        % x_rhat = cat(2, molFrac_chain, nC_chain, x(:,2*numComponents:3*numComponents-1,:));

        % Extract portion of interest only once
        x_sel = x_rhat(idx_start:t, :, :); % size = [T x nParams x N_chains]
        T_eff = size(x_sel,1);

        % Compute chain means (x_j_r_bar)
        x_j_r_bar = squeeze(2/(T_eff-2) * sum(x_sel,1)); % [nParams x N_chains]

        % Compute within-chain variance (W)
        x_bar_expanded = permute(x_j_r_bar, [3 1 2]); % [1 x nParams x N_chains]
        diffs = x_sel - x_bar_expanded;                % broadcast subtraction
        inner_sum_W = squeeze(sum(diffs.^2, [1 3]));   % sum over T and chains
        W = (2 / (N_chains * (T_eff - 2))) * inner_sum_W(:);

        % Compute between-chain variance (B)
        x_j_bar_bar = mean(x_j_r_bar, 2);              % [nParams x 1]
        diffs_B = x_j_r_bar - x_j_bar_bar;             % [nParams x N_chains]
        inner_sum_B = sum(diffs_B.^2, 2);
        B = (T_eff / (2*(N_chains - 1))) * inner_sum_B(:);

        % Compute sigma and R_hat
        sigma = ((T_eff-2)/T_eff) .* W + (2/T_eff) .* B;
        R_hat(counter_check, :) = sqrt( ((N_chains+1)/N_chains) .* sigma ./ W - (T_eff-2)/(N_chains*T_eff) );

        if all(R_hat(counter_check,:) <= R_hat_threshold)  % ---- Convergence is reached ---- %
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
molFrac_chain = molFrac_chain(1:t_convergence,:,:);
x(:, 1:numComponents-1, :) = molFrac_chain;
nC_chain = nC_chain(1:t_convergence,:,:);
x(:, numComponents:2*numComponents-1, :) = nC_chain;
eta_B_star_chain = eta_B_star_chain(1:t_convergence,:,:);
x(:, 2*numComponents:3*numComponents-1, :) = eta_B_star_chain;
p_x = p_x(1:t_convergence,:);

%%% ------- Discard burn-in samples ------- %%%
for l = 1:N_chains
    x_nb((l-1)*(t_convergence-t_burnin)+1:l*(t_convergence-t_burnin), :) = x(t_burnin+1:end, :, l);
    p_x_nb((l-1)*(t_convergence-t_burnin)+1:l*(t_convergence-t_burnin), :, 1) = p_x(t_burnin+1:end, l);
end

fclose(fid); % Close the file