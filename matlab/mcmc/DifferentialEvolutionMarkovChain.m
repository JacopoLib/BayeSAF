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
% It supports parallel tempering, delayed acceptance, snooker jumps, blocked updates,
% an exact multiple-try Gibbs move on carbon atoms and isomers with a uniform prior
% over isomers, and an adaptive burn-in with outlier-chain resets. Convergence is
% monitored via the classic Gelman-Rubin R-hat statistic and the effective sample size.

% Auxiliary parameters:
% Nc: number of surrogate mixture components  [-]

% Inputs:
% 1) families             : (1 x Nc) cell array, with the i-th cell containing a string
% that denotes the i-th hydrocarbon family under consideration
% 2) Posterior_PDF        : function handle representing the posterior PDF
% 3) Posterior_PDF_cheap  : function handle representing a cheap surrogate of the
% posterior PDF (used for delayed acceptance)
% 4) classes              : (1 x Nc) cell array, with the i-th cell containing a structure
% array for the hydrocarbon family and range of number of carbon atoms of the i-th surrogate component
% 5) LowerBound_x         : (1 x Nc) array, with the i-th element representing the
% lower bound for the range the mole fraction of the i-th surrogate mixture
% component can vary within  [-]
% 6) UpperBound_x         : (1 x Nc) array, with the i-th element representing the
% upper bound for the range the mole fraction of the i-th surrogate mixture
% component can vary within  [-]
% 7) n_ranges             : (1 x Nc) cell array, with the i-th cell containing the range
% the number of carbon atoms of the i-th surrogate mixture component can vary within
% 8) LowerBound_eta_B_star: (1 x Nc) array with lower bounds for the normalized topochemical atom indices
% 9) UpperBound_eta_B_star: (1 x Nc) array with upper bounds for the normalized topochemical atom indices
% 10) maxIterations       : maximum number of iterations for each chain
% 11) t_burnin            : number of iterations of the parallel-tempering phase of the burn-in
% 12) scaling_factor_X    : scaling factor for the jump rate concerning molar fractions
% 13) scaling_factor_nc   : scaling factor for the jump rate concerning numbers of carbon atoms
% 14) scaling_factor_eta  : scaling factor for the jump rate concerning topochemical atom indices
% 15) noise_X             : noise parameter for the surrogate mixture
% molar fractions
% 16) noise_nc            : noise parameter for the surrogate mixture
% numbers of carbon atoms
% 17) noise_eta           : noise parameter for the surrogate mixture
% topochemical atom indices
% 18) N_chains            : number of chains
% 19) beta_min            : minimum inverse temperature during the burn-in
% 20) T_ladder            : temperature ladder functional form (1. 'linear', 2. 'geometric')
% 21) swap_freq           : swap frequency for parallel tempering during the burn-in
% 22) n_pairs             : maximum number of chain pairs used to propose the new sample
% 23) p_gibbs             : probability of a Gibbs move
% 24) p_snooker           : probability of performing a snooker jump rather than a DE move
% 25) gibbs_tries         : number of tries of the multiple-try Gibbs move
% 26) n_cold              : number of chains at beta = 1 after the burn-in
% 27) rhat_stop           : R-hat convergence threshold
% 28) rhat_min_samples    : minimum number of iterations after the burn-in
% 29) ess_min             : minimum effective sample size of the molar fractions

% Outputs:
% 1) x            : (t_convergence x 3*Nc-1 x n_cold) array containing the
% samples from the posterior PDF from each chain at beta = 1 along the entire number of
% iterations, concerning the molar fractions, numbers of carbon atoms, and
% topochemical atom indices
% 2) p_x          : (t_convergence x n_cold) array containing the posterior density
% values from each chain at beta = 1 along the entire number of iterations
% 3) x_nb         : (((t_convergence - t_burnin) x n_cold) x 3*Nc-1) array containing the
% samples from the posterior PDF from each chain at beta = 1 discarding the burn-in period,
% concerning the molar fractions, numbers of carbon atoms, and topochemical atom indices
% 4) p_x_nb       : (((t_convergence - t_burnin) x n_cold) x 1) array containing the posterior density
% values from each chain at beta = 1 discarding the burn-in period
% 5) AR           : (t_convergence x n_cold) array containing the acceptance rate from each chain
% at beta = 1 along the entire number of iterations
% 6) R_hat        : floor((t_convergence - t_burnin)/2) x 3*Nc-1 array describing
% the classic R-hat statistic for each parameter, every 2 iterations after the burn-in
% 7) t_convergence: number of iterations the convergence of the DE-MC run
% is reached at ( = maxIterations in case convergence is not reached)
% 8) info         : structure array with the burn-in length (t_burnin), the number of
% outlier-chain resets (n_resets), the swap acceptance rates after the burn-in (swap_rates),
% the inverse temperatures after the burn-in (betas), and the minimum effective sample
% size of the molar fractions at the stop (ESS_min)

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

% Liu, J.S., Liang, F., Wong, W.H., 2000.
% The multiple-try method and local optimization in Metropolis sampling.
% J. Am. Stat. Assoc. 95, 121-134.

% Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., Bürkner, P.-C., 2021.
% Rank-normalization, folding, and localization: an improved R-hat for
% assessing convergence of MCMC.
% Bayesian Anal. 16, 667-718.
% ------------------------------------------------------------------------

function [x, p_x, x_nb, p_x_nb, AR, R_hat, t_convergence, info] = DifferentialEvolutionMarkovChain(families, Posterior_PDF, Posterior_PDF_cheap, classes, LowerBound_x, UpperBound_x, n_ranges, LowerBound_eta_B_star, UpperBound_eta_B_star, maxIterations, t_burnin, scaling_factor_X, scaling_factor_nc, scaling_factor_eta, noise_X, noise_nc, noise_eta, N_chains, beta_min, T_ladder, swap_freq, n_pairs, p_gibbs, p_snooker, gibbs_tries, n_cold, rhat_stop, rhat_min_samples, ess_min)

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

% Parallel tempering after the burn-in
beta_min_sampling = 0.05; % minimum inverse temperature of the hot chains
swap_every_sampling = 1; % swap frequency

% Burn-in
outlier_every = max(50, floor(t_burnin/4)); % iterations between two outlier-chain checks
burnin_cap = 0.75; % maximum burn-in length, as a fraction of maxIterations

% Convergence
rhat_checks = 3; % number of consecutive R-hat checks below rhat_stop
rhat_spacing = 500; % iterations between two R-hat checks

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
nParams = 3*numComponents - 1;
cols_x = 1:numComponents-1;
cols_nc = numComponents:2*numComponents-1;
cols_eta = 2*numComponents:3*numComponents-1;

% Decoding tables of the topochemical atom indices and uniform prior over isomers
isomers = isomer_tables(classes, n_ranges);
Posterior_PDF = @(X_all, N_all, Eta_all) Posterior_PDF(X_all, N_all, Eta_all) + isomer_prior(N_all, Eta_all, isomers, n_ranges);
Posterior_PDF_cheap = @(X_all, N_all, Eta_all) Posterior_PDF_cheap(X_all, N_all, Eta_all) + isomer_prior(N_all, Eta_all, isomers, n_ranges);
posterior_rows = @(f, rows) f(rows(:,cols_x), rows(:,cols_nc), rows(:,cols_eta));

% Topochemical atom indices that are constant within a family (e.g., n-paraffins) are skipped by the R-hat test
rhat_skip = false(1, nParams);
for k = 1:numComponents
    rhat_skip(cols_eta(k)) = numel(unique([classes{k}.eta_B_star_norm])) == 1;
end

x = nan(maxIterations,nParams,N_chains); p_x = nan(maxIterations,N_chains); accept = nan(maxIterations,N_chains); AR = nan(maxIterations,N_chains); R_hat = 1e+18*ones(floor(maxIterations/2), nParams); rhat_history = nan(maxIterations,1);          % Preallocate memory for chains, density, acceptance rate, and R-hat statistic
x_archive = zeros(maxIterations,nParams,N_chains);

% === Initialization ===
% Create initial population by sampling mole fractions uniformly on the simplex and fulfilling constraints from the informative priors
t = 1;
Nneeded = 10 * N_chains;
validSamples = [];

K = numComponents; % total components
K_rest = 2 * K; % numbers of carbon atoms and topochemical atom indices

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
    val_norm = (val_dir(:, 1:K-1) - lb(1:K-1)) ./ (ub(1:K-1) - lb(1:K-1));

    % --- 4) Generate random values for remaining parameters in [0,1]
    n_valid = size(val_norm, 1);
    lhs_rest_cont = rand(n_valid, K_rest);

    % --- 5) Combine into one array
    validSamples = [validSamples; [val_norm, lhs_rest_cont]];

end

X_lhs = validSamples;
phys_lhs = decode_rows(X_lhs, isomers, n_ranges, LowerBound_x, UpperBound_x);
p_X_lhs = posterior_rows(Posterior_PDF, phys_lhs);     % Compute density initial population

N = size(phys_lhs,1);
selectedIdx = zeros(N_chains,1);

% Step 1: randomly pick the first row
selectedIdx(1) = randi(N);

% Step 2: iteratively select rows maximizing min distance to current set
for k = 2:N_chains
    remaining = setdiff(1:N, selectedIdx(1:k-1));
    dists = min(pdist2(phys_lhs(remaining,:), phys_lhs(selectedIdx(1:k-1),:)), [], 2);
    [~, idx_max] = max(dists);
    selectedIdx(k) = remaining(idx_max);
end

selectedIdx = selectedIdx(randperm(N_chains));
X = X_lhs(selectedIdx,:);
phys = phys_lhs(selectedIdx,:);
p_X = p_X_lhs(selectedIdx);

x(t,:,:) = reshape(phys',1,nParams,N_chains);
p_x(t,:) = p_X';   % Store initial states and density

accept(t,:) = 1;                                                                 % First sample has been accepted
AR(t,:) = 100;

for i = 1:N_chains, R(i,1:N_chains-1) = setdiff(1:N_chains,i); end          % R-matrix: index of chains for differential evolution

% Parallel tempering after the burn-in: n_cold chains at beta = 1 and a ladder of hot chains
n_hot = N_chains - n_cold;
if n_cold < 4 || n_hot < 4
    error('At least 4 cold and 4 hot chains are needed.');
end
beta_sampling = [ones(1,n_cold), beta_min_sampling.^((1:n_hot)/n_hot)];
group_of = [ones(1,n_cold), 2*ones(1,n_hot)];
groups = {1:n_cold, n_cold+1:N_chains};
swap_try = zeros(1,n_hot);
swap_acc = zeros(1,n_hot);

% Define inverse temperatures (β) for parallel tempering during the burn-in
if strcmp(T_ladder, 'linear')
    beta_burnin = linspace(1, beta_min, N_chains);   % coldest first, hottest last
elseif strcmp(T_ladder, 'geometric')
    beta_burnin = beta_min.^((0:N_chains-1)/(N_chains-1));
    beta_burnin = sort(beta_burnin,'descend');       % coldest first, hottest last
end

convergence = false; % switch
t_convergence = maxIterations;
counter_check = 0;
ESS_stop = NaN;

% Adaptive burn-in
t_burn = t_burnin;
burnin_done = false;
n_resets = 0;
margin_best = 5*sqrt(nParams/2);

% Create a waitbar
show_waitbar = usejava('desktop');
if show_waitbar
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
end

while ~convergence && t < maxIterations         % Dynamic part: evolution of N chains

    t = t + 1;

    in_sampling = burnin_done;
    if in_sampling
        beta = beta_sampling;
    elseif t > t_burnin
        beta = ones(1,N_chains);
    else
        beta = beta_burnin;
    end

    if show_waitbar
        if getappdata(f,'canceling')
            break
        end
        waitbar(t/maxIterations,f,sprintf(strcat('Number of iterations: ', num2str(t))))
    end

    lambda = unifrnd(-0.1, 0.1, N_chains, 1); % draw N_chains lambda values
    [~, draw] = sort(rand(N_chains-1,N_chains));
    Xp = X;
    phys_p = phys;

    J = zeros(N_chains,1);
    log_pi1_x = zeros(N_chains,1);
    log_pi1_y = zeros(N_chains,1);
    log_alpha1_gibbs = nan(N_chains,1);

    for i = 1:N_chains                                                      % Create proposals

        Xp_temp = Xp(i,:);
        phys_temp = phys(i,:);

        if in_sampling
            partners = groups{group_of(i)};
            partners(partners == i) = [];
        else
            partners = R(i,draw(:,i));
        end

        if rand <= p_gibbs

            % --------------------- GIBBS MOVE --------------------- %
            % Multiple-try Metropolis on (nC, eta) of one component, after a symmetric DE jump of the molar fractions
            k = randi(numComponents);
            jn = cols_nc(k);
            je = cols_eta(k);
            log_pc_x = posterior_rows(Posterior_PDF_cheap, phys(i,:));

            if in_sampling
                partners = partners(randperm(numel(partners)));
            end
            a = partners(1);
            b = partners(2);
            g_X = scaling_factor_X*2.38/sqrt(2*(numComponents-1));
            if rand >= 0.9
                g_X = 1;
            end
            Xp_temp(cols_x) = X(i,cols_x) + (1 - lambda(i)) * g_X * (X(a,cols_x) - X(b,cols_x)) + noise_X*randn(1,numComponents-1);
            [Xp_temp, phys_temp] = fold_fractions(Xp_temp, phys_temp, cols_x, LowerBound_x, UpperBound_x);

            if sum(phys_temp(cols_x)) > 1 - LowerBound_x(end) || sum(phys_temp(cols_x)) < 1 - UpperBound_x(end)
                log_alpha1_gibbs(i) = -Inf;
                Xp(i,:) = Xp_temp;
                phys_p(i,:) = phys_temp;
                continue
            end

            % Tries at the proposed molar fractions
            u_tries = rand(gibbs_tries, 2);
            rows_tries = repmat(phys_temp, gibbs_tries, 1);
            for r = 1:gibbs_tries
                [rows_tries(r,jn), rows_tries(r,je)] = decode_component(u_tries(r,1), u_tries(r,2), isomers{k}, n_ranges{k});
            end
            log_pc_tries = posterior_rows(Posterior_PDF_cheap, rows_tries);
            log_w = beta(i) * log_pc_tries;

            if any(isfinite(log_w))
                w = exp(log_w - max(log_w(isfinite(log_w))));
                w(~isfinite(log_w)) = 0;
                sel = find(rand*sum(w) < cumsum(w), 1);

                % Reference points at the current molar fractions, including the current state
                u_ref = rand(gibbs_tries-1, 2);
                rows_ref = repmat(phys(i,:), gibbs_tries-1, 1);
                for r = 1:gibbs_tries-1
                    [rows_ref(r,jn), rows_ref(r,je)] = decode_component(u_ref(r,1), u_ref(r,2), isomers{k}, n_ranges{k});
                end
                log_w_ref = beta(i) * [posterior_rows(Posterior_PDF_cheap, rows_ref); log_pc_x];
                log_w_ref = log_w_ref(isfinite(log_w_ref));

                log_alpha1_gibbs(i) = log_sum_exp(log_w(isfinite(log_w))) - log_sum_exp(log_w_ref);
                phys_temp([jn je]) = rows_tries(sel,[jn je]);
                Xp_temp([jn je]) = u_tries(sel,:);
                log_pi1_y(i) = log_pc_tries(sel);
                log_pi1_x(i) = log_pc_x;
            else
                log_alpha1_gibbs(i) = -Inf;
            end

            Xp(i,:) = Xp_temp;
            phys_p(i,:) = phys_temp;
            continue

        end

        group_idx = randi(3); % pick one subset of compositional parameters
        if group_idx == 1
            block_group = 'fractions';
        elseif group_idx == 2
            block_group = 'nC';
        elseif group_idx == 3
            block_group = 'eta';
        end

        % --------------------- DE MOVE --------------------- %
        D = randi(n_pairs);
        if in_sampling
            partners = partners(randperm(numel(partners)));
        end
        a = partners(1:D);
        b = partners(D+1:2*D);

        % Calculate default jump rates
        gamma_X = scaling_factor_X*2.38/sqrt(2*D*(numComponents-1));
        gamma_nc = scaling_factor_nc*2.38/sqrt(2*D*numComponents);
        gamma_eta = scaling_factor_eta*2.38/sqrt(2*D*numComponents);
        g_X = randsample([gamma_X 1], 1, 'true', [0.9 0.1]); % Select gamma: 90/10 mix [default 1]
        mask = rand(size(gamma_nc)) < 0.9;   % true with 90% prob
        g_nc = gamma_nc;                     % start with original values
        g_nc(~mask) = 1;                     % replace 10% of them with 1
        g_eta = randsample([gamma_eta 1], 1, 'true', [0.9 0.1]); % Select gamma: 90/10 mix [default 1]

        % Pick random number to determine whether to perform a DE move or a snooker jump
        r_move = rand;
        use_snooker = r_move < p_snooker && t > 10*nParams && t > t_burnin && strcmp(block_group, 'fractions');

        if use_snooker

            % --------------------- SNOOKER JUMP --------------------- %
            rowsToRemove = all(x_archive == 0, [2 3]);
            x_archive_filt = x_archive(~rowsToRemove, :, :);

            if size(x_archive_filt,1) >= 3
                idx_chain = partners(randperm(numel(partners)));
                r1 = idx_chain(1); r2 = idx_chain(2); r3 = idx_chain(3);
                idx_sample = randi(size(x_archive_filt,1), 3, 1);

                % Define reference point xR and direction vectors
                xR = x_archive_filt(idx_sample(1), :, r1);
                z_snooker = xR - X(i,:);
                v = x_archive_filt(idx_sample(2),:,r2) - x_archive_filt(idx_sample(3),:,r3);
                alpha = (z_snooker * v') / (z_snooker * z_snooker' + 1e-8);
                z_proj = alpha * z_snooker;
                g_snooker = unifrnd(1.2, 2.2);
            else
                use_snooker = false;
            end

        end

        switch block_group

            case 'fractions'

                if use_snooker

                    Xp_temp(cols_x) = X(i,cols_x) + g_snooker * z_proj(cols_x) + noise_X * randn(1,numComponents-1);

                    % --- Compute Snooker Jacobian correction ---
                    d = numComponents - 1;
                    zp = xR(cols_x) - Xp_temp(cols_x);
                    epsJ = 1e-12;
                    J(i) = (d - 1) * ( log(norm(z_snooker(cols_x)) + epsJ) - ...
                        log(norm(zp) + epsJ) );

                else

                    for j = cols_x
                        Xp_temp(j) = X(i,j) ...
                            + (1 - lambda(i)) * (g_X * sum(X(a,j)-X(b,j), 1)) ...
                            + noise_X*randn;
                    end

                end

            case 'nC'

                for j = cols_nc
                    Xp_temp(j) = X(i,j) ...
                        + (1 - lambda(i)) * (g_nc * sum(X(a,j)-X(b,j), 1)) ...
                        + noise_nc*randn;
                end

            case 'eta'

                for j = cols_eta
                    Xp_temp(j) = X(i,j) ...
                        + (1 - lambda(i)) * (g_eta * sum(X(a,j)-X(b,j), 1)) ...
                        + noise_eta*randn;
                end

        end

        [Xp_temp, phys_temp] = fold_fractions(Xp_temp, phys_temp, cols_x, LowerBound_x, UpperBound_x);

        if ~strcmp(block_group, 'fractions')
            % Boundary handling for numbers of carbon atoms and topochemical atom indices by folding the parameter space
            for k = 1:numComponents
                Xp_temp(cols_nc(k)) = 1 - abs(1 - mod(Xp_temp(cols_nc(k)), 2));
                Xp_temp(cols_eta(k)) = 1 - abs(1 - mod(Xp_temp(cols_eta(k)), 2));
                [phys_temp(cols_nc(k)), phys_temp(cols_eta(k))] = decode_component(Xp_temp(cols_nc(k)), Xp_temp(cols_eta(k)), isomers{k}, n_ranges{k});
            end
        end

        Xp(i,:) = Xp_temp;
        phys_p(i,:) = phys_temp;

    end

    % -------- Delayed-Acceptance: Stage 1 (cheap) --------
    log_alpha1_vec = beta(:) .* (log_pi1_y - log_pi1_x) + J;
    idx_gibbs = ~isnan(log_alpha1_gibbs);
    log_alpha1_vec(idx_gibbs) = log_alpha1_gibbs(idx_gibbs);
    u1 = log(rand(N_chains,1));
    pass1 = u1 < log_alpha1_vec;   % chains that pass Stage 1

    % default: carry forward previous accept counter for all chains
    accept(t,:) = accept(t-1,:);

    if any(pass1)
        % -------- Delayed-Acceptance: Stage 2 (full, batched only on pass1) --------
        idx_pass = find(pass1);

        % One vectorized full-posterior call for Stage-2 candidates
        p_y_full_pass = posterior_rows(Posterior_PDF, phys_p(idx_pass,:));

        % DA correction: log α2 = β * ( [ℓ(y)-ℓ(x)] - [ℓc(y)-ℓc(x)] )
        delta_full  = p_y_full_pass - p_X(idx_pass);
        delta_cheap = log_pi1_y(idx_pass) - log_pi1_x(idx_pass);
        log_alpha2  = beta(idx_pass)' .* (delta_full - delta_cheap);

        u2 = log(rand(numel(idx_pass),1));
        accept2 = u2 < log_alpha2;

        if any(accept2)
            idx_acc = idx_pass(accept2);

            % Commit accepted proposals: normalized, physical, and full posterior
            X(idx_acc, :) = Xp(idx_acc, :);
            phys(idx_acc, :) = phys_p(idx_acc, :);
            p_X(idx_acc) = p_y_full_pass(accept2);

            % Update accept counters only for finally accepted chains
            accept(t, idx_acc) = accept(t-1, idx_acc) + 1;
        end
    end

    AR(t,:) = 100*(accept(t,:)/(t-1));                                            % Calculate acceptance rate

    % ----- Replica Exchange during the burn-in -----
    if ~in_sampling && t <= t_burnin && mod(t, swap_freq) == 0
        if mod(floor(t / swap_freq), 2) == 0
            pair_start = 1; % (1,2), (3,4), ...
        else
            pair_start = 2; % (2,3), (4,5), ...
        end
        for k = pair_start:2:(N_chains-1)
            i = k; j = k + 1;
            Delta = (beta(i) - beta(j)) * (p_X(j) - p_X(i));
            if (Delta >= 0) || (log(rand) < Delta)
                X([i j], :) = X([j i], :);
                phys([i j], :) = phys([j i], :);
                p_X([i j]) = p_X([j i]);
            end
        end
    end

    % ----- Replica Exchange after the burn-in: one cold chain drawn at random against the hot ladder -----
    if in_sampling && mod(t, swap_every_sampling) == 0
        for lev = 1 + mod(floor(t / swap_every_sampling), 2):2:n_hot
            if lev == 1
                i = randi(n_cold);
            else
                i = n_cold + lev - 1;
            end
            j = n_cold + lev;
            Delta = (beta(i) - beta(j)) * (p_X(j) - p_X(i));
            swap_try(lev) = swap_try(lev) + 1;
            if (Delta >= 0) || (log(rand) < Delta)
                X([i j], :) = X([j i], :);
                phys([i j], :) = phys([j i], :);
                p_X([i j]) = p_X([j i]);
                swap_acc(lev) = swap_acc(lev) + 1;
            end
        end
    end

    %%% ------- Adaptive burn-in: reset outlier chains to good chains drawn at random ------- %%%
    if ~burnin_done && t > t_burnin && mod(t - t_burnin, outlier_every) == 0

        omega = mean(p_x(t-outlier_every:t-1,:), 1);
        threshold = max(omega) - margin_best;
        outliers = find(omega < threshold);
        good = find(omega >= threshold);

        if ~isempty(outliers)
            donors = good(randi(numel(good), 1, numel(outliers)));
            X(outliers,:) = X(donors,:);
            phys(outliers,:) = phys(donors,:);
            p_X(outliers) = p_X(donors);
            n_resets = n_resets + numel(outliers);
            fprintf(fid, ['Chains # ', num2str(outliers), ' below the others (mean log-posterior < ', num2str(threshold, '%.1f'), ') after ', num2str(t), ' iterations. Their states were replaced by those of chains # ', num2str(donors), ', drawn at random among the good chains.\n']);
        elseif t >= 2*t_burnin
            burnin_done = true;
            t_burn = t;
            fprintf(fid, ['No outlier chains after ', num2str(t), ' iterations: the burn-in period ends here.\n']);
        else
            fprintf(fid, ['No outlier chains after ', num2str(t), ' iterations: the burn-in period continues to at least ', num2str(2*t_burnin), ' iterations.\n']);
        end

        if ~burnin_done && t >= burnin_cap*maxIterations
            burnin_done = true;
            t_burn = t;
            fprintf(fid, ['WARNING: maximum burn-in length reached after ', num2str(t), ' iterations with outlier chains still present.\n']);
        end

    end

    x(t,:,:) = reshape(phys',1,nParams,N_chains);
    p_x(t,:) = p_X';       % Append current states and density

    if mod(t,10) == 0
        x_archive(t,:,:) = reshape(X',1,nParams,N_chains);
    end

    %%% ------- Classic R-hat statistic of the chains at beta = 1 every 2 iterations after the burn-in ------- %%%
    if burnin_done && t > t_burn && mod(t - t_burn, 2) == 0
        counter_check = counter_check + 1;
        R_hat(counter_check,:) = rhat_classic(x(:,:,1:n_cold), t_burn, t);
        R_hat(counter_check,rhat_skip) = NaN;
    end

    %%% ------- Check convergence through the R-hat statistic and the effective sample size ------- %%%
    if mod(t,100) == 0 && burnin_done && t - t_burn >= 20

        rc = rhat_classic(x(:,:,1:n_cold), t_burn, t);
        rc(rhat_skip) = 1;
        rc(~isfinite(rc)) = 1e+18;
        rhat_history(t) = max(rc);

        t_checks = t - (0:rhat_checks-1)*rhat_spacing;
        stop = t - t_burn >= rhat_min_samples && all(t_checks >= 1) && all(rhat_history(max(t_checks,1)) <= rhat_stop);

        if stop
            ESS_stop = Inf;
            for j = cols_x
                ESS_stop = min(ESS_stop, ess_bulk(squeeze(x(t_burn:t,j,1:n_cold))));
            end
            stop = ESS_stop >= ess_min;
        end

        if ~show_waitbar
            fprintf('DE-MC iteration %d/%d  AR = %.1f%%  R-hat max = %.3f\n', t, maxIterations, mean(AR(t,:)), max(rc));
        end

        if stop
            t_convergence = t;
            fprintf(fid, ['R-statistic below the critical threshold of ', num2str(rhat_stop), ' after ', num2str(fliplr(t_checks)), ' iterations and minimum effective sample size of the molar fractions equal to ', num2str(ESS_stop, '%.0f'), ': convergence reached after ', num2str(t_convergence), ' iterations within each chain. These samples will be used for posterior analysis after discarding the burn-in period.\n']);
            convergence = true;
        end

    elseif ~show_waitbar && mod(t,100) == 0
        fprintf('DE-MC iteration %d/%d  AR = %.1f%%\n', t, maxIterations, mean(AR(t,:)));
    end

end         % End dynamic part

% Close the waitbar
if show_waitbar
    delete(f);
end

if ~convergence     % ---- Convergence not reached ---- %
    t_convergence = t;
    fprintf(fid, ['Convergence not reached! Posterior analysis will be performed considering the entirety of ', num2str(t_convergence), ' MCMC samples after discarding the burn-in period.\n']);
end

if ~burnin_done
    t_burn = t_convergence - 1;
end

% Only the chains at beta = 1 are returned
x = x(1:t_convergence,:,1:n_cold);
p_x = p_x(1:t_convergence,1:n_cold);
AR = AR(1:t_convergence,1:n_cold);
R_hat = R_hat(1:counter_check,:);

%%% ------- Discard burn-in samples ------- %%%
x_nb = zeros((t_convergence-t_burn)*n_cold, nParams);
p_x_nb = zeros((t_convergence-t_burn)*n_cold, 1);
for l = 1:n_cold
    x_nb((l-1)*(t_convergence-t_burn)+1:l*(t_convergence-t_burn), :) = x(t_burn+1:end, :, l);
    p_x_nb((l-1)*(t_convergence-t_burn)+1:l*(t_convergence-t_burn)) = p_x(t_burn+1:end, l);
end

info.t_burnin = t_burn;
info.n_resets = n_resets;
info.swap_rates = swap_acc ./ max(swap_try, 1);
info.betas = beta_sampling;
info.ESS_min = ESS_stop;

fclose(fid); % Close the file

end

% ------------------------------------------------------------------------

function isomers = isomer_tables(classes, n_ranges)
% For each component and number of carbon atoms: the topochemical atom index decoded from each
% rank of the class-wide sorted list of eta_B_star, and the log-prior correction giving a
% uniform prior over the isomers in place of the measure of their decoder cells
isomers = cell(1, numel(classes));
for k = 1:numel(classes)
    eta_list = [classes{k}.eta_B_star];
    eta_sorted = sort(eta_list);
    eta_norm_list = [classes{k}.eta_B_star_norm];
    nC_list = [classes{k}.nC];
    N = numel(eta_sorted);
    if N > 1
        lo = max(0, ((0:N-1) - 0.5)/(N-1));
        hi = min(1, ((0:N-1) + 0.5)/(N-1));
    else
        lo = 0;
        hi = 1;
    end
    n_range_k = n_ranges{k};
    isomers{k}.decode = zeros(N, numel(n_range_k));
    isomers{k}.eta = cell(1, numel(n_range_k));
    isomers{k}.log_prior = cell(1, numel(n_range_k));
    for j = 1:numel(n_range_k)
        nC_indices = find(nC_list == n_range_k(j));
        for r = 1:N
            [~, idx_neighb] = min(abs(eta_list(nC_indices) - eta_sorted(r)));
            isomers{k}.decode(r,j) = eta_norm_list(nC_indices(idx_neighb));
        end
        eta_j = unique(isomers{k}.decode(:,j))';
        cell_measure = zeros(size(eta_j));
        for e = 1:numel(eta_j)
            cell_measure(e) = sum(hi(isomers{k}.decode(:,j) == eta_j(e)) - lo(isomers{k}.decode(:,j) == eta_j(e)));
        end
        isomers{k}.eta{j} = eta_j;
        isomers{k}.log_prior{j} = -log(nnz(cell_measure > 0)) - log(cell_measure);
    end
end
end

function log_prior = isomer_prior(N_all, Eta_all, isomers, n_ranges)
log_prior = zeros(size(N_all,1), 1);
for r = 1:size(N_all,1)
    for k = 1:size(N_all,2)
        j = find(n_ranges{k} == N_all(r,k), 1);
        e = find(isomers{k}.eta{j} == Eta_all(r,k), 1);
        if isempty(e)
            log_prior(r) = -Inf;
        else
            log_prior(r) = log_prior(r) + isomers{k}.log_prior{j}(e);
        end
    end
end
end

function [nC, eta_B_star_norm] = decode_component(u_nc, u_eta, isomers_k, n_range_k)
% Map normalized coordinates in [0,1] to the number of carbon atoms (equal-width bins) and to the
% topochemical atom index of the isomer closest to the sampled rank of the class-wide eta_B_star list
kk = numel(n_range_k);
idx_nc = min(floor(u_nc * kk) + 1, kk);
nC = n_range_k(idx_nc);
N = size(isomers_k.decode, 1);
idx_eta = max(1, min(round(1 + u_eta * (N - 1)), N));
eta_B_star_norm = isomers_k.decode(idx_eta, idx_nc);
end

function phys = decode_rows(X, isomers, n_ranges, LowerBound_x, UpperBound_x)
numComponents = numel(isomers);
phys = zeros(size(X));
for i = 1:size(X,1)
    for j = 1:numComponents-1
        phys(i,j) = X(i,j) * (UpperBound_x(j) - LowerBound_x(j)) + LowerBound_x(j);
    end
    for k = 1:numComponents
        [phys(i,numComponents-1+k), phys(i,2*numComponents-1+k)] = decode_component(X(i,numComponents-1+k), X(i,2*numComponents-1+k), isomers{k}, n_ranges{k});
    end
end
end

function [Xp_temp, phys_temp] = fold_fractions(Xp_temp, phys_temp, cols_x, LowerBound_x, UpperBound_x)
for j = cols_x
    if Xp_temp(j) < 0
        Xp_temp(j) = max(1 - abs(Xp_temp(j)), 0);
    elseif Xp_temp(j) > 1
        Xp_temp(j) = min(abs(Xp_temp(j) - 1), 1);
    end
    phys_temp(j) = Xp_temp(j) * (UpperBound_x(j) - LowerBound_x(j)) + LowerBound_x(j);
end
end

function s = log_sum_exp(v)
m = max(v);
s = m + log(sum(exp(v - m)));
end

function R_hat = rhat_classic(x, t_burn, t)
% Classic Gelman-Rubin R-hat over the second half of the iterations after the burn-in
x = x(t_burn + floor((t - t_burn)/2):t, :, :);
n = size(x,1);
W = mean(var(x, 0, 1), 3);
B = n * var(mean(x, 1), 0, 3);
R_hat = sqrt(((n-1)/n * W + B/n) ./ W);
end

function ESS = ess_bulk(draws)
% Bulk effective sample size of draws (iterations x chains), Vehtari et al. (2021)
if all(draws(:) == draws(1))
    ESS = Inf;
    return
end
z = norminv((reshape(tiedrank(draws(:)), size(draws)) - 0.375) / (numel(draws) + 0.25));
n = floor(size(z,1)/2);
z = [z(1:n,:), z(n+1:2*n,:)];
m = size(z,2);
zc = z - mean(z,1);
F = fft(zc, 2*n, 1);
acov = real(ifft(F .* conj(F), [], 1));
acov = acov(1:n,:) / n;
W = mean(acov(1,:) * n / (n-1));
B = n * var(mean(z,1));
var_plus = (n-1)/n * W + B/n;
rho = 1 - (W - mean(acov,2)) / var_plus;
rho(1) = 1;
P = rho(1:2:n-1) + rho(2:2:n);
k_neg = find(P < 0, 1);
if ~isempty(k_neg)
    P = P(1:k_neg-1);
end
if isempty(P)
    P = 1;
else
    P = cummin(P);
end
tau = -1 + 2*sum(P);
ESS = m*n / max(tau, 1/log10(m*n));
end
