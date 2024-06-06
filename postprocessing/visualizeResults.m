function visualizeResults(fullData, families, chain, AR, R_hat, chain_reshaped, t_convergence, maxIterations, n_ranges, numComponents, variable_names, x_MAP, eta_B_star_MAP, minT_array, maxT_array, pressure_distillation, fuel_name)

% DESCRIPTION:
% The visualizeResults function provides tools to analyze and visualize
% results from the Bayesian inference analysis.

% Auxiliary parameters:
% Nd           : number of experimental measurements for the thermophysical property under consideration  [-]
% Np           : number of thermophysical properties targeted during the formulation
% of the surrogate mixture  [-]
% N_chains     : number of chains

% Inputs:
% 1) fullData              : (1 x Np) cell array, with the i-th cell being a (Nd x 3) array
% that contains the experimental measurements about the thermophysical property under
% consideration (independent variable, dependent thermophysical property, standard deviation). 
% Note that the i-th cell is a (Nd x 4) array in case the distillation curve is the i-th
% thermophysical property under consideration.
% 2) families              : (1 x numComponents) cell array, with the i-th cell containing a string
% that denotes the i-th hydrocarbon family under consideration
% 3) chain                 : (t_convergence x 3*numComponents-1 x N_chains) array containing the
% samples from the posterior PDF from each chain along the entire number of
% iterations, concerning the molar fractions, numbers of carbon atoms, and
% topochemical atom indices
% 4) AR                    : (t_convergence x N_chains) array containing the acceptance rate from each chain
% along the entire number of iterations
% 5) R_hat                 : (t_convergence-0.1*floor(maxIterations)) x 3*numComponents-1 array containing the 
% acceptance rate from each chain along the entire number of iterations
% 6) chain_reshaped        : (t_convergence/2 x 3*numComponents-1 x N_chains) array containing the
% samples from the posterior PDF from each chain discarding the burn-in period, concerning
% the molar fractions, numbers of carbon atoms, and topochemical atom indices
% 7) t_convergence         : number of iterations denoting the number of iterations the convergence of the DE-MC run
% is reached at according to the R-hat statistic ( = maxIterations in case
% convergence is not reached)
% 8) maxIterations         : maximum number of iterations for each chain
% 9) n_ranges              : (1 x numComponents) cell array, with the i-th cell containing the range
% the number of carbon atoms of the i-th surrogate mixture component can vary within
% 10) numComponents        : number of surrogate mixture components  [-]
% 11) variable_names       : (1 x Np) cell array, with the i-th cell containing a string
% that denotes the i-th thermophysical property under consideration
% 12) x_MAP                : (1 x numComponents) array containing the molar fractions of the MAP surrogate components
% 13) eta_B_star_MAP       : topochemical atom indices characterizing the MAP surrogate components
% 14) minT_array           : (1 x Np) array, with the i-th element denoting the minimum
% temperature in the dataset of the i-th thermophysical property
% 15) maxT_array           : (1 x Np) array, with the i-th element denoting the maximum
% temperature in the dataset of the i-th thermophysical property
% 16) pressure_distillation: pressure the distillation curve is computed at  [Pa]
% 17) fuel_name            : plot legend string denoting real fuel experimental data

% ------------------------------------------------------------------------%
% Contributors/Copyright
% 2024 Jacopo Liberatori, jacopo.liberatori@uniroma1.it
% 2024 Davide Cavalieri,  davide.cavalieri@uniroma1.it
% Department of Mechanical and Aerospace Engineering (DIMA)
% Sapienza University of Rome
% ------------------------------------------------------------------------%

% Plot settings
fontsize = 16;
fontname = 'Times';
line = 1.3;
pix = 500;
blue  = [0.0039 0.451 0.741];
green = [0.45, 0.80, 0.69];
red   = [0.6470 0.09019 0.1843];
gray = [0.7 0.7 0.7];
num_levels = 256;

% DE-MC parameters
N_chains = size(chain, 3);

% ---------- Plot the marginal and joint PDFs of the molar fractions ---------- %
for m = 1:numComponents-1
    chain_xm = chain_reshaped(1:(floor(t_convergence/2))*N_chains,m);
    % Use kernel density estimation to estimate the marginal PDF of the m-th component's molar fraction
    [f_x, x_values] = ksdensity(chain_xm);

    % Plot the marginal PDF
    figure('visible', 'off'); box on;
    plot(x_values, f_x, 'LineWidth', 2, 'Color',blue);
    xlabelName = strcat('$X_', num2str(m), '$');
    ylabelName = 'Marginal Probability';
    xlabelHandle = xlabel(xlabelName, 'Interpreter', 'latex');
    ylabelHandle = ylabel(ylabelName);
    set(xlabelHandle, 'FontSize', fontsize, 'FontName', fontname);
    set(ylabelHandle, 'FontSize', fontsize, 'FontName', fontname);
    set(0,'defaultfigurecolor',[1 1 1]);
    axesHandle = gca;
    set(axesHandle, 'FontSize', fontsize, 'FontName', fontname, 'Linewidth', line);
    set(gcf, 'Position', [100, 100, pix, pix/1.3801]);
    % Set the x-axis limit dynamically
    xlim([max(0, min(x_values)), min(1, max(x_values))]);
    % Specify the folder path and filename
    folderPathMarginal = './Figures/PosteriorPDF/MarginalProbability/';
    % Create the folder if it doesn't exist
    if ~isfolder(folderPathMarginal)
        mkdir(folderPathMarginal);
    end
    filenameMarginalXm = strcat(folderPathMarginal, 'Marginal_X', num2str(m));
    print(filenameMarginalXm, '-dpng', '-r600');

    for p = m+1:numComponents-1
        chain_xp = chain_reshaped(1:(floor(t_convergence/2))*N_chains,p);
        % Use 2D kernel density estimation to estimate the joint PDF of the m-th and p-th components' molar fractions
        [xm_joint, xp_joint] = meshgrid(linspace(min(chain_xm), max(chain_xm), 100), linspace(min(chain_xp), max(chain_xp), 100));
        positions = [xm_joint(:), xp_joint(:)];
        f_joint = ksdensity([chain_xm, chain_xp], positions);
        % Reshape the estimated PDF
        f_joint = reshape(f_joint, size(xm_joint));

        % Plot the joint 2D PDF
        figure('visible', 'off'); box on;
        contourf(xm_joint, xp_joint, f_joint/max(f_joint(:)), num_levels, 'LineColor', 'none');
        xlabelName = strcat('$X_', num2str(m), '$');
        ylabelName = strcat('$X_', num2str(p), '$');
        xlabelHandle = xlabel(xlabelName, 'Interpreter', 'latex');
        ylabelHandle = ylabel(ylabelName, 'Interpreter', 'latex');
        colormap(slanCM('rainbow'));
        hcb = colorbar;
        ylabel(hcb, 'Normalized Joint Probability', 'FontSize', fontsize, 'FontName', fontname);
        set(xlabelHandle, 'FontSize', fontsize, 'FontName', fontname);
        set(ylabelHandle, 'FontSize', fontsize, 'FontName', fontname);
        set(0,'defaultfigurecolor',[1 1 1]);
        axesHandle = gca;
        set(axesHandle, 'FontSize', fontsize, 'FontName', fontname, 'Linewidth', line);
        set(gcf, 'Position', [100, 100, pix, pix/1.3801]);
        folderPathJoint2D = './Figures/PosteriorPDF/Joint2D/';
        % Create the folder if it doesn't exist
        if ~isfolder(folderPathJoint2D)
            mkdir(folderPathJoint2D);
        end
        filenameJoint2D_X = strcat(folderPathJoint2D, 'Joint2D_X', num2str(m), '_X', num2str(p));
        print(filenameJoint2D_X, '-dpng', '-r600');
    end
end

% ---------- Plot the marginal and joint PDFs of the numbers of carbon atoms ---------- %
customColors = lines(numComponents+1);
for m = numComponents:2*numComponents-1
    % Plot the histograms for the number of carbon atoms of the m-th component (1D marginal probability)
    figure('visible', 'off'); box on;
    num_bins_m = numel(n_ranges{m-numComponents+1});
    edges = linspace(min(n_ranges{m-numComponents+1}), max(n_ranges{m-numComponents+1}), num_bins_m + 1);  % Define edges to explicitly specify bin boundaries
    histograms = histcounts(chain_reshaped(1:(floor(t_convergence/2))*N_chains,m), edges, 'Normalization', 'probability');
    bar(histograms, 'FaceColor', customColors(m-numComponents+1, :));
    bin_centers = 1:1:num_bins_m;
    xticks(bin_centers);
    xticklabels(n_ranges{m-numComponents+1});
    xlabelName = strcat('$n_{C,', num2str(m-numComponents+1), '}$');
    ylabelName = 'Marginal Probability';
    xlabelHandle = xlabel(xlabelName, 'Interpreter', 'latex');
    ylabelHandle = ylabel(ylabelName);
    set(xlabelHandle, 'FontSize', fontsize, 'FontName', fontname);
    set(ylabelHandle, 'FontSize', fontsize, 'FontName', fontname);
    set(0,'defaultfigurecolor',[1 1 1]);
    axesHandle = gca;
    set(axesHandle, 'FontSize', fontsize, 'FontName', fontname, 'Linewidth', line);
    set(gcf, 'Position', [100, 100, pix, pix/1.3801]);
    filenameMarginal_ncm = strcat(folderPathMarginal, 'Marginal_nc', num2str(m-numComponents+1));
    print(filenameMarginal_ncm, '-dpng', '-r600');

    % 2D histogram plot representing the joint PDF of the numbers of carbon atoms of the m-th and p-th components
    for p = m+1:2*numComponents-1
        figure('visible', 'off'); box on;
        num_bins_p = numel(n_ranges{p-numComponents+1});
        m_edges = linspace(min(n_ranges{m-numComponents+1}), max(n_ranges{m-numComponents+1}), num_bins_m + 1);
        p_edges = linspace(min(n_ranges{p-numComponents+1}), max(n_ranges{p-numComponents+1}), num_bins_p + 1);
        histogram2D = histcounts2(chain_reshaped(1:(floor(t_convergence/2))*N_chains,m), chain_reshaped(1:(floor(t_convergence/2))*N_chains,p), m_edges, p_edges, 'Normalization', 'probability');
        imagesc(histogram2D');
        bin_centers_m = 1:1:num_bins_m;
        bin_centers_p = 1:1:num_bins_p;
        xticks(bin_centers_m);
        xticklabels(n_ranges{m-numComponents+1});
        yticks(bin_centers_p);
        yticklabels(n_ranges{p-numComponents+1});
        xlabelName = strcat('$n_{C,', num2str(m-numComponents+1), '}$');
        ylabelName = strcat('$n_{C,', num2str(p-numComponents+1), '}$');
        xlabelHandle = xlabel(xlabelName, 'Interpreter', 'latex');
        ylabelHandle = ylabel(ylabelName, 'Interpreter', 'latex');
        colormap(slanCM('YlGn'));
        hcb = colorbar;
        ylabel(hcb, 'Normalized Joint Probability', 'FontSize', fontsize, 'FontName', fontname);
        set(xlabelHandle, 'FontSize', fontsize, 'FontName', fontname);
        set(ylabelHandle, 'FontSize', fontsize, 'FontName', fontname);
        set(0,'defaultfigurecolor',[1 1 1]);
        axesHandle = gca;
        set(axesHandle, 'FontSize', fontsize, 'FontName', fontname, 'Linewidth', line);
        set(gcf, 'Position', [100, 100, pix, pix/1.3801]);
        folderPathJoint2D = './Figures/PosteriorPDF/Joint2D/';
        % Create the folder if it doesn't exist
        if ~isfolder(folderPathJoint2D)
            mkdir(folderPathJoint2D);
        end
        filenameJoint2D_nc = strcat(folderPathJoint2D, 'Joint2D_nc', num2str(m-numComponents+1), '_nc', num2str(p-numComponents+1));
        print(filenameJoint2D_nc, '-dpng', '-r600');
    end
end

% ---------- Plot the marginal and joint PDFs of the topochemical atom indices ---------- %
for m = 2*numComponents:3*numComponents-1
    chain_xm = chain_reshaped(1:(floor(t_convergence/2))*N_chains,m);
    % Use kernel density estimation to estimate the marginal PDF of the m-th component's topochemical atom index
    [f_x, x_values] = ksdensity(chain_xm);

    % Plot the marginal PDF
    figure('visible', 'off'); box on;
    plot(x_values, f_x, 'LineWidth', 2, 'Color',blue);
    xlabelName = strcat('$\eta_{B,', num2str(m-2*numComponents+1), '}^{*}$');
    ylabelName = 'Marginal Probability';
    xlabelHandle = xlabel(xlabelName, 'Interpreter', 'latex');
    ylabelHandle = ylabel(ylabelName);
    set(xlabelHandle, 'FontSize', fontsize, 'FontName', fontname);
    set(ylabelHandle, 'FontSize', fontsize, 'FontName', fontname);
    set(0,'defaultfigurecolor',[1 1 1]);
    axesHandle = gca;
    set(axesHandle, 'FontSize', fontsize, 'FontName', fontname, 'Linewidth', line);
    set(gcf, 'Position', [100, 100, pix, pix/1.3801]);
    % Set the x-axis limit dynamically
    xlim([max(0, min(x_values)), min(1, max(x_values))]);
    % Specify the folder path and filename
    folderPathMarginal = './Figures/PosteriorPDF/MarginalProbability/';
    % Create the folder if it doesn't exist
    if ~isfolder(folderPathMarginal)
        mkdir(folderPathMarginal);
    end
    filenameMarginalEta = strcat(folderPathMarginal, 'Marginal_eta', num2str(m-2*numComponents+1));
    print(filenameMarginalEta, '-dpng', '-r600');

    for p = m+1:3*numComponents-1
        chain_xp = chain_reshaped(1:(floor(t_convergence/2))*N_chains,p);
        % Use 2D kernel density estimation to estimate the joint PDF of the m-th and p-th components' topochemical atom indices
        [xm_joint, xp_joint] = meshgrid(linspace(min(chain_xm), max(chain_xm), 100), linspace(min(chain_xp), max(chain_xp), 100));
        positions = [xm_joint(:), xp_joint(:)];
        f_joint = ksdensity([chain_xm, chain_xp], positions);
        % Reshape the estimated PDF
        f_joint = reshape(f_joint, size(xm_joint));

        % Plot the joint 2D PDF
        figure('visible', 'off'); box on;
        contourf(xm_joint, xp_joint, f_joint/max(f_joint(:)), num_levels, 'LineColor', 'none');
        xlabelName = strcat('$\eta_{B,', num2str(m-2*numComponents+1), '}^{*}$');
        ylabelName = strcat('$\eta_{B,', num2str(p-2*numComponents+1), '}^{*}$');
        xlabelHandle = xlabel(xlabelName, 'Interpreter', 'latex');
        ylabelHandle = ylabel(ylabelName, 'Interpreter', 'latex');
        colormap(slanCM('rainbow'));
        hcb = colorbar;
        ylabel(hcb, 'Normalized Joint Probability', 'FontSize', fontsize, 'FontName', fontname);
        set(xlabelHandle, 'FontSize', fontsize, 'FontName', fontname);
        set(ylabelHandle, 'FontSize', fontsize, 'FontName', fontname);
        set(0,'defaultfigurecolor',[1 1 1]);
        axesHandle = gca;
        set(axesHandle, 'FontSize', fontsize, 'FontName', fontname, 'Linewidth', line);
        set(gcf, 'Position', [100, 100, pix, pix/1.3801]);
        folderPathJoint2D = './Figures/PosteriorPDF/Joint2D/';
        % Create the folder if it doesn't exist
        if ~isfolder(folderPathJoint2D)
            mkdir(folderPathJoint2D);
        end
        filenameJoint2D_eta = strcat(folderPathJoint2D, 'Joint2D_eta', num2str(m-2*numComponents+1), '_eta', num2str(p-2*numComponents+1));
        print(filenameJoint2D_eta, '-dpng', '-r600');
    end
end

% ---------- Plot the joint PDFs of molar fractions and numbers of carbon atoms ---------- %
for m = 1:numComponents-1
    for p = numComponents:2*numComponents-1
        figure('visible', 'off');
        map = slanCM('bwr');
        colors = map(n_ranges{p-numComponents+1}); % Get a colormap with distinct colors
        daviolinplot(chain_reshaped(1:(floor(t_convergence/2))*N_chains,m), 'groups', chain_reshaped(1:(floor(t_convergence/2))*N_chains,p), 'outliers', 0, 'xtlabels', n_ranges{p-numComponents+1}, 'color', colors, 'scatter', 0, 'jitter', 0, 'box', 2, 'boxcolors', gray, 'scattercolors', 'same', 'boxspacing', 1);
        box on;
        xlabelName = strcat('$n_{C,', num2str(p-numComponents+1), '}$');
        ylabelName = strcat('$X_', num2str(m), '$');
        xlabelHandle = xlabel(xlabelName, 'Interpreter', 'latex');
        ylabelHandle = ylabel(ylabelName, 'Interpreter', 'latex');
        set(xlabelHandle, 'FontSize', fontsize, 'FontName', fontname);
        set(ylabelHandle, 'FontSize', fontsize, 'FontName', fontname);
        set(0,'defaultfigurecolor',[1 1 1]);
        axesHandle = gca;
        set(axesHandle, 'FontSize', fontsize, 'FontName', fontname, 'Linewidth', line);
        set(gcf, 'Position', [100, 100, pix, pix/1.3801]);
        filenameJoint2D_X_nc = strcat(folderPathJoint2D, 'Joint2D_X', num2str(m), '_nc', num2str(p-numComponents+1));
        print(filenameJoint2D_X_nc, '-dpng', '-r600');
    end
end

% ---------- Plot the joint PDFs of molar fractions and topochemical atom indices ---------- %
for m = 1:numComponents-1
    chain_xm = chain_reshaped(1:(floor(t_convergence/2))*N_chains,m);
    for p = 2*numComponents:3*numComponents-1
        chain_xp = chain_reshaped(1:(floor(t_convergence/2))*N_chains,p);
        % Use 2D kernel density estimation to estimate the joint PDF of the m-th component's molar fraction and p-th component's topochemical atom index
        [xm_joint, xp_joint] = meshgrid(linspace(min(chain_xm), max(chain_xm), 100), linspace(min(chain_xp), max(chain_xp), 100));
        positions = [xm_joint(:), xp_joint(:)];
        f_joint = ksdensity([chain_xm, chain_xp], positions);
        % Reshape the estimated PDF
        f_joint = reshape(f_joint, size(xm_joint));

        % Plot the joint 2D PDF
        figure('visible', 'off'); box on;
        contourf(xm_joint, xp_joint, f_joint/max(f_joint(:)), num_levels, 'LineColor', 'none');
        xlabelName = strcat('$X_', num2str(m), '$');
        ylabelName = strcat('$\eta_{B,', num2str(p-2*numComponents+1), '}^{*}$');
        xlabelHandle = xlabel(xlabelName, 'Interpreter', 'latex');
        ylabelHandle = ylabel(ylabelName, 'Interpreter', 'latex');
        colormap(slanCM('rainbow'));
        hcb = colorbar;
        ylabel(hcb, 'Normalized Joint Probability', 'FontSize', fontsize, 'FontName', fontname);
        set(xlabelHandle, 'FontSize', fontsize, 'FontName', fontname);
        set(ylabelHandle, 'FontSize', fontsize, 'FontName', fontname);
        set(0,'defaultfigurecolor',[1 1 1]);
        axesHandle = gca;
        set(axesHandle, 'FontSize', fontsize, 'FontName', fontname, 'Linewidth', line);
        set(gcf, 'Position', [100, 100, pix, pix/1.3801]);
        filenameJoint2D_X_eta = strcat(folderPathJoint2D, 'Joint2D_X', num2str(m), '_eta', num2str(p-2*numComponents+1));
        print(filenameJoint2D_X_eta, '-dpng', '-r600');
    end
end

% ---------- Plot the joint PDFs of topochemical atom indices and numbers of carbon atoms ---------- %
for m = 2*numComponents:3*numComponents-1
    for p = numComponents:2*numComponents-1
        figure('visible', 'off');
        map = slanCM('bwr');
        colors = map(n_ranges{p-numComponents+1}); % Get a colormap with distinct colors
        daviolinplot(chain_reshaped(1:(floor(t_convergence/2))*N_chains,m), 'groups', chain_reshaped(1:(floor(t_convergence/2))*N_chains,p), 'outliers', 0, 'xtlabels', n_ranges{p-numComponents+1}, 'color', colors, 'scatter', 0, 'jitter', 0, 'box', 2, 'boxcolors', gray, 'scattercolors', 'same', 'boxspacing', 1);
        box on;
        xlabelName = strcat('$n_{C,', num2str(p-numComponents+1), '}$');
        ylabelName = strcat('$\eta_{B,', num2str(m-2*numComponents+1), '}^{*}$');
        xlabelHandle = xlabel(xlabelName, 'Interpreter', 'latex');
        ylabelHandle = ylabel(ylabelName, 'Interpreter', 'latex');
        set(xlabelHandle, 'FontSize', fontsize, 'FontName', fontname);
        set(ylabelHandle, 'FontSize', fontsize, 'FontName', fontname);
        set(0,'defaultfigurecolor',[1 1 1]);
        axesHandle = gca;
        set(axesHandle, 'FontSize', fontsize, 'FontName', fontname, 'Linewidth', line);
        set(gcf, 'Position', [100, 100, pix, pix/1.3801]);
        filenameJoint2D_eta_nc = strcat(folderPathJoint2D, 'Joint2D_eta', num2str(m-2*numComponents+1), '_nc', num2str(p-numComponents+1));
        print(filenameJoint2D_eta_nc, '-dpng', '-r600');
    end
end

% ---------- Trace plots of the molar fractions by the Markov Chain algorithm ---------- %
for m = 1:numComponents-1
    figure('visible', 'off');
    hold on;
    box on;
    colors = parula(N_chains); % Get a colormap with N_chains distinct colors
    for i = 1:N_chains
        plot(chain(:,m,i), 'LineWidth', 0.75, 'Color',colors(i, :));
    end
    xlabelName = strcat('MCMC Iterations');
    ylabelName = strcat('Trace of $X_', num2str(m), '$');
    xlabelHandle = xlabel(xlabelName);
    ylabelHandle = ylabel(ylabelName, 'Interpreter', 'latex');
    set(xlabelHandle, 'FontSize', fontsize, 'FontName', fontname);
    set(ylabelHandle, 'FontSize', fontsize, 'FontName', fontname);
    set(0,'defaultfigurecolor',[1 1 1]);
    axesHandle = gca;
    set(axesHandle, 'FontSize', fontsize, 'FontName', fontname, 'Linewidth', line);
    set(gcf, 'Position', [100, 100, pix, pix/1.3801]);
    folderPathMCMC = './Figures/MCMC_check/';
    % Create the folder if it doesn't exist
    if ~isfolder(folderPathMCMC)
        mkdir(folderPathMCMC);
    end
    filenameTraceXm = strcat(folderPathMCMC, 'Trace_X', num2str(m));
    print(filenameTraceXm, '-dpng', '-r600');
end

% ---------- Trace plots of the numbers of carbon atoms by the Markov Chain algorithm ---------- %
for m = numComponents:2*numComponents-1
    figure('visible', 'off');
    hold on;
    box on;
    colors = parula(N_chains); % Get a colormap with N_chains distinct colors
    for i = 1:N_chains
        plot(chain(:,m,i), 'LineWidth', 0.75, 'Color',colors(i, :));
    end
    xlabelName = strcat('MCMC Iterations');
    ylabelName = strcat('$n_{C,', num2str(m-numComponents+1), '}$');
    xlabelHandle = xlabel(xlabelName);
    ylabelHandle = ylabel(ylabelName, 'Interpreter', 'latex');
    set(xlabelHandle, 'FontSize', fontsize, 'FontName', fontname);
    set(ylabelHandle, 'FontSize', fontsize, 'FontName', fontname);
    set(0,'defaultfigurecolor',[1 1 1]);
    axesHandle = gca;
    set(axesHandle, 'FontSize', fontsize, 'FontName', fontname, 'Linewidth', line);
    set(gcf, 'Position', [100, 100, pix, pix/1.3801]);
    filenameTrace_ncm = strcat(folderPathMCMC, 'Trace_nc', num2str(m-numComponents+1));
    print(filenameTrace_ncm, '-dpng', '-r600');
end

% ---------- Trace plots of the topochemical atom indices by the Markov Chain algorithm ---------- %
for m = 2*numComponents:3*numComponents-1
    figure('visible', 'off');
    hold on;
    box on;
    colors = parula(N_chains); % Get a colormap with N_chains distinct colors
    for i = 1:N_chains
        plot(chain(:,m,i), 'LineWidth', 0.75, 'Color',colors(i, :));
    end
    
    xlabelName = strcat('MCMC Iterations');
    ylabelName = strcat('Trace of $\eta_{B,', num2str(m-2*numComponents+1), '}^{*}$');
    xlabelHandle = xlabel(xlabelName);
    ylabelHandle = ylabel(ylabelName, 'Interpreter', 'latex');
    set(xlabelHandle, 'FontSize', fontsize, 'FontName', fontname);
    set(ylabelHandle, 'FontSize', fontsize, 'FontName', fontname);
    set(0,'defaultfigurecolor',[1 1 1]);
    axesHandle = gca;
    set(axesHandle, 'FontSize', fontsize, 'FontName', fontname, 'Linewidth', line);
    set(gcf, 'Position', [100, 100, pix, pix/1.3801]);
    filenameTrace_ncm = strcat(folderPathMCMC, 'Trace_eta', num2str(m-2*numComponents+1));
    print(filenameTrace_ncm, '-dpng', '-r600');
end

% ---------- Plot the acceptance rate by the Markov Chain algorithm ---------- %
figure('visible', 'off');
hold on;
box on;
colors = parula(N_chains); % Get a colormap with N_chains distinct colors
for i = 1:N_chains
    plot(AR(:,i), 'LineWidth', 0.75, 'Color',colors(i, :));
end
ylim([0 100]);
xlabelName = strcat('MCMC Iterations');
ylabelName = strcat('Acceptance Rate [%]');
xlabelHandle = xlabel(xlabelName);
ylabelHandle = ylabel(ylabelName);
set(xlabelHandle, 'FontSize', fontsize, 'FontName', fontname);
set(ylabelHandle, 'FontSize', fontsize, 'FontName', fontname);
set(0,'defaultfigurecolor',[1 1 1]);
axesHandle = gca;
set(axesHandle, 'FontSize', fontsize, 'FontName', fontname, 'Linewidth', line);
set(gcf, 'Position', [100, 100, pix, pix/1.3801]);
filenameAcceptanceRate = strcat(folderPathMCMC, 'AcceptanceRate');
print(filenameAcceptanceRate, '-dpng', '-r600');

% ---------- Plot the R-hat statistic for each parameter by the Markov Chain algorithm ---------- %
figure('visible', 'off');
hold on;
box on;
abscissa = linspace(0.1*floor(maxIterations), t_convergence, t_convergence-0.1*floor(maxIterations)+1);
colors = parula(3*numComponents-1); % Get a colormap with N_chains distinct colors
for i = 1:3*numComponents-1
    plot(abscissa(1:2:end), R_hat(1:2:end,i), 'LineWidth', 0.75, 'Color',colors(i, :));
end
y_line = 1.2; % R-hat statistic threshold denoting convergence
plot(abscissa, y_line*ones(size(abscissa)), 'r', 'LineWidth', 1.5, 'LineStyle', '--');
xlim([0.1*floor(maxIterations) t_convergence]);
ylim([0 5]);
xlabelName = strcat('MCMC Iterations');
ylabelName = strcat('R-hat [-]');
xlabelHandle = xlabel(xlabelName);
ylabelHandle = ylabel(ylabelName);
set(xlabelHandle, 'FontSize', fontsize, 'FontName', fontname);
set(ylabelHandle, 'FontSize', fontsize, 'FontName', fontname);
set(0,'defaultfigurecolor',[1 1 1]);
axesHandle = gca;
set(axesHandle, 'FontSize', fontsize, 'FontName', fontname, 'Linewidth', line);
set(gcf, 'Position', [100, 100, pix, pix/1.3801]);
filenameRhat = strcat(folderPathMCMC, 'Rhat');
print(filenameRhat, '-dpng', '-r600');

% ---------- Check the performance of the MAP surrogate against experimental data ---------- %

% Indices of the surrogate components concerning the numbers of carbon
% atoms and topochemical atom indices
index_n_eta_MAP = zeros(1, numComponents-1);
% Initialize empty cell array to store the classes
classes_MAP = cell(1, numComponents-1);

for i = 1:numComponents
    classes_MAP{i} = Hydrocarbons(families{i}, n_ranges{i});
    index_n_eta_MAP(i) = findIndex_eta(classes_MAP{i}, eta_B_star_MAP(i));
end

for h = 1:numel(variable_names)

    temperature_range = linspace(minT_array(h), maxT_array(h), 25);
    property_MAP = zeros(numel(temperature_range),1);

    for i = 1:numel(temperature_range)
        if strcmp(variable_names{h}, 'rho') || strcmp(variable_names{h}, 'kappa') || strcmp(variable_names{h}, 'specificHeat') || strcmp(variable_names{h}, 'latentHeat') || strcmp(variable_names{h}, 'vaporPressure') || strcmp(variable_names{h}, 'sigma')
            property_MAP(i) = ThermophysicalProperties_LiquidMixture(variable_names{h}, temperature_range(i), x_MAP(1:numComponents-1), classes_MAP, index_n_eta_MAP, pressure_distillation);
        end
        if strcmp(variable_names{h}, 'mu') || strcmp(variable_names{h}, 'nu')
            property_MAP(i) = 1e+6*ThermophysicalProperties_LiquidMixture(variable_names{h}, temperature_range(i), x_MAP(1:numComponents-1), classes_MAP, index_n_eta_MAP, pressure_distillation);
        end
    end

    if strcmp(variable_names{h}, 'distillation')
        [volumeFraction, property_MAP] = DistillationCurve(x_MAP(1:numComponents-1), classes_MAP, index_n_eta_MAP, pressure_distillation);
    end

    numSamples = size(chain_reshaped, 1);
    mean_model = zeros(numel(temperature_range), 1);
    var_model = zeros(numel(temperature_range), 1);
    std_model = zeros(numel(temperature_range), 1);

    for i = 1:numel(temperature_range)
        model_output_samples = zeros(numSamples, 1);
        for j = 1:numSamples
            index_n_eta_samples = zeros(1,numComponents);
            for k = 1:numComponents
                index_n_eta_samples(k) = findIndex_eta(classes_MAP{k}, chain_reshaped(j, 2*numComponents+k-1));
            end
            if strcmp(variable_names{h}, 'rho')
                if j > 1 && isequal(chain_reshaped(j,:),chain_reshaped(j-1,:))
                    model_output_samples(j) = model_output_samples(j-1);
                else
                    model_output_samples(j) = ThermophysicalProperties_LiquidMixture(variable_names{h}, temperature_range(i), chain_reshaped(j, 1:numComponents-1), classes_MAP, index_n_eta_samples, pressure_distillation);
                end
            end
            if strcmp(variable_names{h}, 'mu') || strcmp(variable_names{h}, 'nu')
                if j > 1 && isequal(chain_reshaped(j,:),chain_reshaped(j-1,:))
                    model_output_samples(j) = model_output_samples(j-1);
                else
                    model_output_samples(j) = 1e+6*ThermophysicalProperties_LiquidMixture(variable_names{h}, temperature_range(i), chain_reshaped(j, 1:numComponents-1), classes_MAP, index_n_eta_samples, pressure_distillation);
                end
            end
            if strcmp(variable_names{h}, 'kappa')
                if j > 1 && isequal(chain_reshaped(j,:),chain_reshaped(j-1,:))
                    model_output_samples(j) = model_output_samples(j-1);
                else
                    model_output_samples(j) = ThermophysicalProperties_LiquidMixture(variable_names{h}, temperature_range(i), chain_reshaped(j, 1:numComponents-1), classes_MAP, index_n_eta_samples, pressure_distillation);
                end
            end
            if strcmp(variable_names{h}, 'specificHeat')
                if j > 1 && isequal(chain_reshaped(j,:),chain_reshaped(j-1,:))
                    model_output_samples(j) = model_output_samples(j-1);
                else
                    model_output_samples(j) = ThermophysicalProperties_LiquidMixture(variable_names{h}, temperature_range(i), chain_reshaped(j, 1:numComponents-1), classes_MAP, index_n_eta_samples, pressure_distillation);
                end
            end
            if strcmp(variable_names{h}, 'latentHeat')
                if j > 1 && isequal(chain_reshaped(j,:),chain_reshaped(j-1,:))
                    model_output_samples(j) = model_output_samples(j-1);
                else
                    model_output_samples(j) = ThermophysicalProperties_LiquidMixture(variable_names{h}, temperature_range(i), chain_reshaped(j, 1:numComponents-1), classes_MAP, index_n_eta_samples, pressure_distillation);
                end
            end
            if strcmp(variable_names{h}, 'vaporPressure')
                if j > 1 && isequal(chain_reshaped(j,:),chain_reshaped(j-1,:))
                    model_output_samples(j) = model_output_samples(j-1);
                else
                    model_output_samples(j) = ThermophysicalProperties_LiquidMixture(variable_names{h}, temperature_range(i), chain_reshaped(j, 1:numComponents-1), classes_MAP, index_n_eta_samples, pressure_distillation);
                end
            end
            if strcmp(variable_names{h}, 'sigma')
                if j > 1 && isequal(chain_reshaped(j,:),chain_reshaped(j-1,:))
                    model_output_samples(j) = model_output_samples(j-1);
                else
                    model_output_samples(j) = ThermophysicalProperties_LiquidMixture(variable_names{h}, temperature_range(i), chain_reshaped(j, 1:numComponents-1), classes_MAP, index_n_eta_samples, pressure_distillation);
                end
            end
        end
        mean_model(i) = mean(model_output_samples);
        var_model(i) = var(model_output_samples);
        std_model(i) = sqrt(var_model(i));
    end

    if strcmp(variable_names{h}, 'distillation')
        mean_model = zeros(numel(0:1:100), 1);
        var_model = zeros(numel(0:1:100), 1);
        std_model = zeros(numel(0:1:100), 1);
        model_output_samples = zeros(numSamples, numel(0:1:100));
        for j = 1:numSamples
            index_n_eta_samples = zeros(1,numComponents);
            for k = 1:numComponents
                index_n_eta_samples(k) = findIndex_eta(classes_MAP{k}, chain_reshaped(j, 2*numComponents+k-1));
            end
            if j > 1 && isequal(chain_reshaped(j,:),chain_reshaped(j-1,:))
                model_output_samples(j,:) = model_output_samples(j-1,:);
            else
                [volumeFraction, property] = DistillationCurve(chain_reshaped(j, 1:numComponents-1), classes_MAP, index_n_eta_samples, pressure_distillation);
                model_output_samples(j,:) = property;
            end
        end
        for i = 1:numel(0:1:100)
            mean_model(i) = mean(model_output_samples(:,i));
            var_model(i) = var(model_output_samples(:,i));
            std_model(i) = sqrt(var_model(i));
        end
    end

    figure('visible', 'off');
    xlabelName = 'Temperature [K]';
    if strcmp(variable_names{h}, 'rho')
        ylabelName = 'Density [kg/m^3]';
    end
    if strcmp(variable_names{h}, 'mu')
        ylabelName = 'Dynamic Viscosity [μPa s]';
    end
    if strcmp(variable_names{h}, 'nu')
        ylabelName = 'Kinematic Viscosity [mm^2/s]';
    end
    if strcmp(variable_names{h}, 'kappa')
        ylabelName = 'Thermal Conductivity [W/m/K]';
    end
    if strcmp(variable_names{h}, 'specificHeat')
        ylabelName = 'Specific Heat Capacity [J/kg/K]';
    end
    if strcmp(variable_names{h}, 'latentHeat')
        ylabelName = 'Enthalpy of Vaporization [J/kg]';
    end
    if strcmp(variable_names{h}, 'vaporPressure')
        ylabelName = 'Vapor Pressure [Pa]';
    end
    if strcmp(variable_names{h}, 'sigma')
        ylabelName = 'Surface Tension [N/m]';
    end
    if strcmp(variable_names{h}, 'distillation')
        ylabelName = 'Temperature [K]';
    end
    xlabelHandle = xlabel(xlabelName);
    ylabelHandle = ylabel(ylabelName);
    set(xlabelHandle, 'FontSize', fontsize, 'FontName', fontname);
    set(ylabelHandle, 'FontSize', fontsize, 'FontName', fontname);
    set(0,'defaultfigurecolor',[1 1 1]);
    axesHandle = gca;
    set(axesHandle, 'FontSize', fontsize, 'FontName', fontname, 'Linewidth', line);
    set(gcf, 'Position', [100, 100, pix, pix/1.3801]);
    hold on;
    if strcmp(variable_names{h}, 'rho') || strcmp(variable_names{h}, 'mu') || strcmp(variable_names{h}, 'nu') || strcmp(variable_names{h}, 'kappa') || strcmp(variable_names{h}, 'specificHeat') || strcmp(variable_names{h}, 'latentHeat') || strcmp(variable_names{h}, 'vaporPressure') || strcmp(variable_names{h}, 'sigma')
        xlim([0.95*temperature_range(1) 1.05*temperature_range(end)]);
        ylim([0.98*min(min(mean_model - 3*std_model),min(fullData{h}(:,2))) 1.02*max(max(mean_model + 3*std_model),max(fullData{h}(:,2)))]);
        plot(temperature_range, property_MAP, 'LineWidth', 1.1, 'Color', 'black', 'DisplayName', 'MAP Surrogate');
        fill_x = [temperature_range, fliplr(temperature_range)];
        fill_y = [mean_model' + 3*std_model', fliplr(mean_model' - 3*std_model')];
        fill(fill_x, fill_y, blue, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', '\mu \pm 3\sigma');
    elseif strcmp(variable_names{h}, 'distillation')
        xlim([volumeFraction(1) volumeFraction(end)]);
        ylim([0.99*min(min(mean_model - 3*std_model),min(fullData{h}(:,2))) 1.01*max(max(mean_model + 3*std_model),max(fullData{h}(:,2)))]);
        plot(volumeFraction, property_MAP, 'LineWidth', 1.1, 'Color','black', 'DisplayName', 'MAP Surrogate');
        fill_x = [volumeFraction, fliplr(volumeFraction)];
        fill_y = [mean_model' + 3*std_model', fliplr(mean_model' - 3*std_model')];
        fill(fill_x, fill_y, blue, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', '\mu \pm 3\sigma');
    end
    if strcmp(variable_names{h}, 'mu') || strcmp(variable_names{h}, 'nu')
        plot(fullData{h}(:,1), 1e+6*fullData{h}(:,2), 'o', 'MarkerSize', 6, 'MarkerFace', 'white', 'MarkerEdge', 'k', 'LineWidth', 1, 'DisplayName', fuel_name);
    elseif strcmp(variable_names{h}, 'distillation')
        xlabelName = 'Recovered Volume Fraction [%]';
        xlabelHandle = xlabel(xlabelName);
        set(xlabelHandle, 'FontSize', fontsize, 'FontName', fontname);
        plot(fullData{h}(:,1), fullData{h}(:,2), 'o', 'MarkerSize', 6, 'MarkerFace', 'white', 'MarkerEdge', 'k', 'LineWidth', 1, 'DisplayName', fuel_name);
    else
        plot(fullData{h}(:,1), fullData{h}(:,2), 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'white', 'MarkerEdge', 'k', 'LineWidth', 1, 'DisplayName', fuel_name);
    end
    hold off;
    box on;
    % Set the 'Interpreter' property for the legend entry you want to display differently
    h_legend = legend('show', 'Location', 'best');
    legend_entries = h_legend.String;  % Get legend entries
    entry_index = find(strcmp(legend_entries, 'MAP Surrogate'));  % Find index of 'MAP Surrogate'
    if ~isempty(entry_index)
        h_legend.EntryContainer.Children(entry_index).Label.Interpreter = 'none';  % Set interpreter for 'MAP Surrogate'
    end
    folderPathPropertyCheck = './Figures/PropertyCheck_MAP/';
    % Create the folder if it doesn't exist
    if ~isfolder(folderPathPropertyCheck)
        mkdir(folderPathPropertyCheck);
    end
    filenamePropertyCheck = strcat(folderPathPropertyCheck, variable_names{h});
    print(filenamePropertyCheck, '-dpng', '-r600');
end

end