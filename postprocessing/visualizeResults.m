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

function visualizeResults(fullData, families, classes, chain, AR, R_hat, chain_reshaped, t_convergence, t_burnin, n_ranges, numComponents, variable_names, x_MAP, nc_MAP, eta_B_star_MAP, minT_array, maxT_array, pressure_distillation, fuel_name, band_percentiles, confidence_width)

% DESCRIPTION:
% The visualizeResults function provides tools to analyze and visualize
% results from the Bayesian inference analysis.

% Auxiliary parameters:
% Nd           : number of experimental measurements for the thermophysical property under consideration  [-]
% Np           : number of thermophysical properties targeted during the formulation
% of the surrogate mixture  [-]
% N_chains     : number of DE-MC chains
% maxIterations: maximum number of iterations for each chain in the DE-MC algorithm

% Inputs:
% 1) fullData              : (1 x Np) cell array, with the i-th cell being a (Nd x 3) array
% that contains the experimental measurements about the thermophysical property under
% consideration (independent variable, dependent thermophysical property, standard deviation).
% Note that the i-th cell is a (Nd x 4) array in case the distillation curve is the i-th
% thermophysical property under consideration.
% 2) families              : (1 x numComponents) cell array, with the i-th cell containing a string
% that denotes the i-th hydrocarbon family under consideration
% 3) classes               : (1 x numComponents) cell array, with the i-th cell containing a structure array
% for the hydrocarbon family and range of number of carbon atoms of the i-th surrogate component
% 4) chain                 : (t_convergence x 3*numComponents-1 x N_chains) array containing the
% samples from the posterior PDF from each chain along the entire number of
% iterations, concerning the molar fractions, numbers of carbon atoms, and
% topochemical atom indices
% 5) AR                    : (t_convergence x N_chains) array containing the acceptance rate from each chain
% along the entire number of iterations
% 6) R_hat                 : floor((t_convergence-t_burnin)/2) x 3*numComponents-1 array containing the
% acceptance rate from each chain along the entire number of iterations
% 7) chain_reshaped        : ((t_convergence - t_burnin)*N_chains x 3*numComponents-1) array containing the
% samples from the posterior PDF from each chain discarding the burn-in period, concerning
% the molar fractions, numbers of carbon atoms, and topochemical atom indices
% 8) t_convergence         : number of iterations denoting the number of iterations the convergence of the DE-MC run
% is reached at according to the R-hat statistic ( = maxIterations in case convergence is not reached)
% 9) t_burnin              : number of iterations to be discarded as burn-in
% 10) n_ranges             : (1 x numComponents) cell array, with the i-th cell containing the range
% the number of carbon atoms of the i-th surrogate mixture component can vary within
% 11) numComponents        : number of surrogate mixture components  [-]
% 12) variable_names       : (1 x Np) cell array, with the i-th cell containing a string
% that denotes the i-th thermophysical property under consideration
% 13) x_MAP                : (1 x numComponents) array containing the molar fractions of the MAP surrogate components
% 14) nc_MAP               : (1 x numComponents) array containing the numbers of carbon atoms of the MAP surrogate components
% 15) eta_B_star_MAP       : topochemical atom indices characterizing the MAP surrogate components
% 16) minT_array           : (1 x Np) array, with the i-th element denoting the minimum
% temperature in the dataset of the i-th thermophysical property
% 17) maxT_array           : (1 x Np) array, with the i-th element denoting the maximum
% temperature in the dataset of the i-th thermophysical property
% 18) pressure_distillation: pressure the distillation curve is computed at  [Pa]
% 19) fuel_name            : plot legend string denoting real fuel experimental data
% 20) band_percentiles     : string denoting whether the band representing the 90% confidence interval is colored with
% percentiles ('True') or not ('False')
% 21) confidence_width     : width of the confidence interval

rng('default'); % to ensure reproducibility

% Plot settings
fontsize = 20;
fontname = 'Times';
line = 1.3;
pix = 550;
green = [0.45, 0.80, 0.69];
blue  = [0.0039 0.451 0.741];
orange = [0.8500, 0.3250, 0.0980];
gray = [0.7 0.7 0.7];
num_levels = 256;

% DE-MC parameters
N_chains = size(chain, 3);

% ---------- Plot the marginal and joint PDFs of the molar fractions ---------- %
for m = 1:numComponents-1

    chain_xm = chain_reshaped(:,m);
    % Use kernel density estimation to estimate the marginal PDF of the m-th component's molar fraction
    [f_x, x_values] = ksdensity(chain_xm);

    % Plot the marginal PDF
    t = tiledlayout(1,1); t.Padding = 'none'; t.TileSpacing = 'none'; nexttile
    figure('visible', 'off'); box on;
    plot(x_values, f_x, 'LineWidth', 2, 'Color',blue);
    xlabelName = strcat('$X_', num2str(m), '$');
    ylabelName = 'Marginal Probability';
    xlabelHandle = xlabel(xlabelName);
    ylabelHandle = ylabel(ylabelName);
    set(xlabelHandle, 'FontSize', fontsize, 'FontName', fontname, 'Interpreter', 'latex');
    set(ylabelHandle, 'FontSize', fontsize, 'FontName', fontname, 'Interpreter', 'latex');
    set(0,'defaultfigurecolor',[1 1 1]);
    axesHandle = gca;
    set(axesHandle, 'FontSize', fontsize, 'FontName', fontname, 'Linewidth', line);
    set(gcf, 'Position', [100, 100, pix, pix/1.3801]);
    % Get current x-axis limits
    currentLimits = xlim;
    if currentLimits(1) < 0
        % Set the new x-axis limits, keeping the upper bound the same
        xlim([0 currentLimits(2)]);
    end
    % Specify the folder path and filename
    folderPathMarginal = './Figures/PosteriorPDF/MarginalProbability/';
    % Create the folder if it doesn't exist
    if ~isfolder(folderPathMarginal)
        mkdir(folderPathMarginal);
    end
    filenameMarginalXm = strcat(folderPathMarginal, 'Marginal_X', num2str(m));
    print(filenameMarginalXm, '-dpng', '-r600');

    for p = m+1:numComponents-1
        chain_xp = chain_reshaped(:,p);
        % Use 2D kernel density estimation to estimate the joint PDF of the m-th and p-th components' molar fractions
        [xm_joint, xp_joint] = meshgrid(linspace(min(chain_xm), max(chain_xm), 100), linspace(min(chain_xp), max(chain_xp), 100));
        positions = [xm_joint(:), xp_joint(:)];
        f_joint = ksdensity([chain_xm, chain_xp], positions);
        % Reshape the estimated PDF
        f_joint = reshape(f_joint, size(xm_joint));

        % Plot the joint 2D PDF
        t = tiledlayout(1,1); t.Padding = 'none'; t.TileSpacing = 'none'; nexttile
        figure('visible', 'off'); box on;
        contourf(xm_joint, xp_joint, f_joint/max(f_joint(:)), num_levels, 'LineColor', 'none');
        xlabelName = strcat('$X_', num2str(m), '$');
        ylabelName = strcat('$X_', num2str(p), '$');
        xlabelHandle = xlabel(xlabelName);
        ylabelHandle = ylabel(ylabelName);
        colormap(flipud(slanCM('torch')));
        hcb = colorbar;
        ylabel(hcb, 'Normalized Joint Probability', 'FontSize', fontsize, 'FontName', fontname, 'Interpreter', 'latex');
        set(xlabelHandle, 'FontSize', fontsize, 'FontName', fontname, 'Interpreter', 'latex');
        set(ylabelHandle, 'FontSize', fontsize, 'FontName', fontname, 'Interpreter', 'latex');
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
    t = tiledlayout(1,1); t.Padding = 'none'; t.TileSpacing = 'none'; nexttile
    figure('visible', 'off'); box on;
    num_bins_m = numel(n_ranges{m-numComponents+1});
    edges = linspace(min(n_ranges{m-numComponents+1}), max(n_ranges{m-numComponents+1}), num_bins_m + 1);  % Define edges to explicitly specify bin boundaries
    histograms = histcounts(chain_reshaped(:,m), edges, 'Normalization', 'probability');
    bar(histograms, 'FaceColor', customColors(m-numComponents+1, :));
    bin_centers = 1:1:num_bins_m;
    xticks(bin_centers);
    xticklabels(n_ranges{m-numComponents+1});
    xlabelName = strcat('$n_{C,', num2str(m-numComponents+1), '}$');
    ylabelName = 'Marginal Probability';
    xlabelHandle = xlabel(xlabelName);
    ylabelHandle = ylabel(ylabelName);
    set(xlabelHandle, 'FontSize', fontsize, 'FontName', fontname, 'Interpreter', 'latex');
    set(ylabelHandle, 'FontSize', fontsize, 'FontName', fontname, 'Interpreter', 'latex');
    set(0,'defaultfigurecolor',[1 1 1]);
    axesHandle = gca;
    set(axesHandle, 'FontSize', fontsize, 'FontName', fontname, 'Linewidth', line);
    set(gcf, 'Position', [100, 100, pix, pix/1.3801]);
    filenameMarginal_ncm = strcat(folderPathMarginal, 'Marginal_nc', num2str(m-numComponents+1));
    print(filenameMarginal_ncm, '-dpng', '-r600');

    % 2D histogram plot representing the joint PDF of the numbers of carbon atoms of the m-th and p-th components
    for p = m+1:2*numComponents-1

        t = tiledlayout(1,1); t.Padding = 'none'; t.TileSpacing = 'none'; nexttile
        figure('visible', 'off'); box on;
        num_bins_p = numel(n_ranges{p-numComponents+1});
        m_edges = linspace(min(n_ranges{m-numComponents+1}), max(n_ranges{m-numComponents+1}), num_bins_m + 1);
        p_edges = linspace(min(n_ranges{p-numComponents+1}), max(n_ranges{p-numComponents+1}), num_bins_p + 1);
        histogram2D = histcounts2(chain_reshaped(:,m), chain_reshaped(:,p), m_edges, p_edges, 'Normalization', 'probability');
        imagesc(histogram2D');
        hold on;
        % Vertical edges
        for i = 1:num_bins_m + 1
            xline(i - 0.5, 'k', 'LineWidth', 0.2);  % Add vertical lines at the edges
        end
        % Horizontal edges
        for j = 1:num_bins_p + 1
            yline(j - 0.5, 'k', 'LineWidth', 0.2);  % Add horizontal lines at the edges
        end
        bin_centers_m = 1:1:num_bins_m;
        bin_centers_p = 1:1:num_bins_p;
        xticks(bin_centers_m);
        xticklabels(n_ranges{m-numComponents+1});
        yticks(bin_centers_p);
        yticklabels(n_ranges{p-numComponents+1});
        xlabelName = strcat('$n_{C,', num2str(m-numComponents+1), '}$');
        ylabelName = strcat('$n_{C,', num2str(p-numComponents+1), '}$');
        xlabelHandle = xlabel(xlabelName);
        ylabelHandle = ylabel(ylabelName);
        colormap(slanCM('binary'));
        hcb = colorbar;
        ylabel(hcb, 'Normalized Joint Probability', 'FontSize', fontsize, 'FontName', fontname, 'Interpreter', 'latex');
        set(xlabelHandle, 'FontSize', fontsize, 'FontName', fontname, 'Interpreter', 'latex');
        set(ylabelHandle, 'FontSize', fontsize, 'FontName', fontname, 'Interpreter', 'latex');
        set(0,'defaultfigurecolor',[1 1 1]);
        axesHandle = gca;
        set(axesHandle, 'FontSize', fontsize, 'FontName', fontname, 'Linewidth', line, 'YDir', 'normal');
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
customColors = lines(2*numComponents+1);
for m = 2*numComponents:3*numComponents-1

    % Find unique values
    [unique_values_m, ~, ~] = unique(chain_reshaped(:,m));

    % Plot the histograms for the topochemical atom indices of the m-th component (1D marginal probability)
    t = tiledlayout(1,1); t.Padding = 'none'; t.TileSpacing = 'none'; nexttile
    figure('visible', 'off'); box on;
    num_bins_m = numel(unique_values_m);
    mid_values_m = (unique_values_m(1:end-1) + unique_values_m(2:end)) / 2;
    edges = [unique_values_m(1), mid_values_m', unique_values_m(end)];  % Define edges to explicitly specify bin boundaries
    histograms = histcounts(chain_reshaped(:,m), edges, 'Normalization', 'probability');
    [~, top_indices_eta{m-2*numComponents+1}] = maxk(histograms, min(10,numel(unique_values_m)));
    top_eta{m-2*numComponents+1} = unique_values_m(top_indices_eta{m-2*numComponents+1});
    bar(histograms, 'FaceColor', customColors(m-2*numComponents+1, :));
    bin_centers_m = 1:round(num_bins_m/5):num_bins_m;
    xticks(bin_centers_m);
    xticklabels(round(unique_values_m(bin_centers_m),2));
    xlabelName = strcat('$\eta_{B,', num2str(m-2*numComponents+1), '}^{*}$');
    ylabelName = 'Marginal Probability';
    xlabelHandle = xlabel(xlabelName);
    ylabelHandle = ylabel(ylabelName);
    set(xlabelHandle, 'FontSize', fontsize, 'FontName', fontname, 'Interpreter', 'latex');
    set(ylabelHandle, 'FontSize', fontsize, 'FontName', fontname, 'Interpreter', 'latex');
    set(0,'defaultfigurecolor',[1 1 1]);
    axesHandle = gca;
    set(axesHandle, 'FontSize', fontsize, 'FontName', fontname, 'Linewidth', line);
    set(gcf, 'Position', [100, 100, pix, pix/1.3801]);
    filenameMarginalEta = strcat(folderPathMarginal, 'Marginal_eta', num2str(m-2*numComponents+1));
    print(filenameMarginalEta, '-dpng', '-r600');

    for p = m+1:3*numComponents-1

        % Find unique values
        [unique_values_p, ~, ~] = unique(chain_reshaped(:,p));

        t = tiledlayout(1,1); t.Padding = 'none'; t.TileSpacing = 'none'; nexttile
        figure('visible', 'off'); box on;
        num_bins_p = numel(unique_values_p);
        mid_values_p = (unique_values_p(1:end-1) + unique_values_p(2:end)) / 2;
        m_edges = [unique_values_m(1), mid_values_m', unique_values_m(end)];
        p_edges = [unique_values_p(1), mid_values_p', unique_values_p(end)];
        histogram2D = histcounts2(chain_reshaped(:,m), chain_reshaped(:,p), m_edges, p_edges, 'Normalization', 'probability');
        imagesc(histogram2D');
        hold on;
        if numel(unique_values_m) < 50 && numel(unique_values_p) < 50
            % Vertical edges
            for i = 1:num_bins_m + 1
                xline(i - 0.5, 'k', 'LineWidth', 0.2);  % Add vertical lines at the edges
            end
            % Horizontal edges
            for j = 1:num_bins_p + 1
                yline(j - 0.5, 'k', 'LineWidth', 0.2);  % Add horizontal lines at the edges
            end
        end
        bin_centers_m = 1:round(num_bins_m/5):num_bins_m;
        bin_centers_p = 1:round(num_bins_p/5):num_bins_p;
        xticks(bin_centers_m);
        xticklabels(round(unique_values_m(bin_centers_m),2));
        yticks(bin_centers_p);
        yticklabels(round(unique_values_p(bin_centers_p),2));
        xlabelName = strcat('$\eta_{B,', num2str(m-2*numComponents+1), '}^{*}$');
        ylabelName = strcat('$\eta_{B,', num2str(p-2*numComponents+1), '}^{*}$');
        xlabelHandle = xlabel(xlabelName);
        ylabelHandle = ylabel(ylabelName);
        colormap(slanCM('binary'));
        hcb = colorbar;
        ylabel(hcb, 'Normalized Joint Probability', 'FontSize', fontsize, 'FontName', fontname, 'Interpreter', 'latex');
        set(xlabelHandle, 'FontSize', fontsize, 'FontName', fontname, 'Interpreter', 'latex');
        set(ylabelHandle, 'FontSize', fontsize, 'FontName', fontname, 'Interpreter', 'latex');
        set(0,'defaultfigurecolor',[1 1 1]);
        axesHandle = gca;
        set(axesHandle, 'FontSize', fontsize, 'FontName', fontname, 'Linewidth', line, 'YDir', 'normal');
        set(gcf, 'Position', [100, 100, pix, pix/1.3801]);
        filenameJoint2D_eta = strcat(folderPathJoint2D, 'Joint2D_eta', num2str(m-2*numComponents+1), '_eta', num2str(p-2*numComponents+1));
        print(filenameJoint2D_eta, '-dpng', '-r600');
    end
end

% ---------- Plot the joint PDFs of molar fractions and numbers of carbon atoms ---------- %
for m = 1:numComponents-1
    for p = numComponents:2*numComponents-1
        t = tiledlayout(1,1); t.Padding = 'none'; t.TileSpacing = 'none'; nexttile
        figure('visible', 'off');
        daviolinplot(chain_reshaped(:,m), 'groups', chain_reshaped(:,p), 'outliers', 0, 'xtlabels', unique(chain_reshaped(:,p)), 'scatter', 0, 'jitter', 0, 'box', 2, 'boxcolors', gray, 'scattercolors', 'same', 'boxspacing', 1);
        box on;
        bin_centers = 1:1:num_bins_m;
        xticks(1:1:numel(unique(chain_reshaped(:,p))));
        xticklabels(unique(chain_reshaped(:,p)));
        xlabelName = strcat('$n_{C,', num2str(p-numComponents+1), '}$');
        ylabelName = strcat('$X_', num2str(m), '$');
        xlabelHandle = xlabel(xlabelName);
        ylabelHandle = ylabel(ylabelName);
        set(xlabelHandle, 'FontSize', fontsize, 'FontName', fontname, 'Interpreter', 'latex');
        set(ylabelHandle, 'FontSize', fontsize, 'FontName', fontname, 'Interpreter', 'latex');
        set(0,'defaultfigurecolor',[1 1 1]);
        axesHandle = gca;
        set(axesHandle, 'FontSize', fontsize, 'FontName', fontname, 'Linewidth', line);
        set(gcf, 'Position', [100, 100, pix, pix/1.3801]);
        % Get current y-axis limits
        currentLimits = ylim;
        if currentLimits(1) < 0
            % Set the new y-axis limits, keeping the upper bound the same
            ylim([0 currentLimits(2)]);
        end
        filenameJoint2D_X_nc = strcat(folderPathJoint2D, 'Joint2D_X', num2str(m), '_nc', num2str(p-numComponents+1));
        print(filenameJoint2D_X_nc, '-dpng', '-r600');
    end
end

% ---------- Plot the joint PDFs of molar fractions and topochemical atom indices ---------- %
for m = 1:numComponents-1
    for p = 2*numComponents:3*numComponents-1
        top_eta_indices = find(ismember(chain_reshaped(:,p), top_eta{p-2*numComponents+1}));
        chain_reshaped_m = chain_reshaped(top_eta_indices,m);
        chain_reshaped_p = chain_reshaped(top_eta_indices,p);
        t = tiledlayout(1,1); t.Padding = 'none'; t.TileSpacing = 'none'; nexttile
        figure('visible', 'off');
        daviolinplot(chain_reshaped_m, 'groups', chain_reshaped_p, 'outliers', 0, 'xtlabels', round(unique(chain_reshaped_p),2), 'scatter', 0, 'jitter', 0, 'box', 2, 'boxcolors', gray, 'scattercolors', 'same', 'boxspacing', 1);
        box on;
        xlabelName = strcat('$\eta_{B,', num2str(p-2*numComponents+1), '}^{*}$');
        ylabelName = strcat('$X_', num2str(m), '$');
        xlabelHandle = xlabel(xlabelName);
        ylabelHandle = ylabel(ylabelName);
        set(xlabelHandle, 'FontSize', fontsize, 'FontName', fontname, 'Interpreter', 'latex');
        set(ylabelHandle, 'FontSize', fontsize, 'FontName', fontname, 'Interpreter', 'latex');
        set(0,'defaultfigurecolor',[1 1 1]);
        axesHandle = gca;
        set(axesHandle, 'FontSize', fontsize, 'FontName', fontname, 'Linewidth', line);
        set(gcf, 'Position', [100, 100, pix, pix/1.3801]);
        % Get current y-axis limits
        currentLimits = ylim;
        if currentLimits(1) < 0
            % Set the new y-axis limits, keeping the upper bound the same
            ylim([0 currentLimits(2)]);
        end
        filenameJoint2D_X_eta = strcat(folderPathJoint2D, 'Joint2D_X', num2str(m), '_eta', num2str(p-2*numComponents+1));
        print(filenameJoint2D_X_eta, '-dpng', '-r600');
    end
end

% ---------- Plot the joint PDFs of topochemical atom indices and numbers of carbon atoms ---------- %
for m = 2*numComponents:3*numComponents-1

    % Find unique values
    [unique_values_m, ~, ~] = unique(chain_reshaped(:,m));

    for p = numComponents:2*numComponents-1

        % Find unique values
        [unique_values_p, ~, ~] = unique(chain_reshaped(:,p));

        t = tiledlayout(1,1); t.Padding = 'none'; t.TileSpacing = 'none'; nexttile
        figure('visible', 'off'); box on;
        num_bins_m = numel(unique_values_m);
        mid_values_m = (unique_values_m(1:end-1) + unique_values_m(2:end)) / 2;
        m_edges = [unique_values_m(1), mid_values_m', unique_values_m(end)];
        num_bins_p = numel(unique_values_p);
        p_edges = linspace(min(unique_values_p), max(unique_values_p), num_bins_p + 1);
        histogram2D = histcounts2(chain_reshaped(:,m), chain_reshaped(:,p), m_edges, p_edges, 'Normalization', 'probability');
        imagesc(histogram2D');
        hold on;
        if numel(unique_values_m) < 50
            % Vertical edges
            for i = 1:num_bins_m + 1
                xline(i - 0.5, 'k', 'LineWidth', 0.2);  % Add vertical lines at the edges
            end
            % Horizontal edges
            for j = 1:num_bins_p + 1
                yline(j - 0.5, 'k', 'LineWidth', 0.2);  % Add horizontal lines at the edges
            end
        end
        bin_centers_m = 1:round(num_bins_m/5):num_bins_m;
        xticks(bin_centers_m);
        xticklabels(round(unique_values_m(bin_centers_m),2));
        bin_centers_p = 1:1:num_bins_p;
        yticks(bin_centers_p);
        yticklabels(n_ranges{p-numComponents+1});
        xlabelName = strcat('$\eta_{B,', num2str(m-2*numComponents+1), '}^{*}$');
        ylabelName = strcat('$n_{C,', num2str(p-numComponents+1), '}$');
        xlabelHandle = xlabel(xlabelName);
        ylabelHandle = ylabel(ylabelName);
        colormap(slanCM('binary'));
        hcb = colorbar;
        ylabel(hcb, 'Normalized Joint Probability', 'FontSize', fontsize, 'FontName', fontname, 'Interpreter', 'latex');
        set(xlabelHandle, 'FontSize', fontsize, 'FontName', fontname, 'Interpreter', 'latex');
        set(ylabelHandle, 'FontSize', fontsize, 'FontName', fontname, 'Interpreter', 'latex');
        set(0,'defaultfigurecolor',[1 1 1]);
        axesHandle = gca;
        set(axesHandle, 'FontSize', fontsize, 'FontName', fontname, 'Linewidth', line, 'YDir', 'normal');
        set(gcf, 'Position', [100, 100, pix, pix/1.3801]);
        filenameJoint2D_eta_nc = strcat(folderPathJoint2D, 'Joint2D_eta', num2str(m-2*numComponents+1), '_nc', num2str(p-numComponents+1));
        print(filenameJoint2D_eta_nc, '-dpng', '-r600');
    end
end

% ---------- Trace plots of the molar fractions by the Markov Chain algorithm ---------- %
for m = 1:numComponents-1
    t = tiledlayout(1,1); t.Padding = 'none'; t.TileSpacing = 'none'; nexttile
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
    ylabelHandle = ylabel(ylabelName);
    set(xlabelHandle, 'FontSize', fontsize, 'FontName', fontname, 'Interpreter', 'latex');
    set(ylabelHandle, 'FontSize', fontsize, 'FontName', fontname, 'Interpreter', 'latex');
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
    t = tiledlayout(1,1); t.Padding = 'none'; t.TileSpacing = 'none'; nexttile
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
    ylabelHandle = ylabel(ylabelName);
    set(xlabelHandle, 'FontSize', fontsize, 'FontName', fontname, 'Interpreter', 'latex');
    set(ylabelHandle, 'FontSize', fontsize, 'FontName', fontname, 'Interpreter', 'latex');
    set(0,'defaultfigurecolor',[1 1 1]);
    axesHandle = gca;
    set(axesHandle, 'FontSize', fontsize, 'FontName', fontname, 'Linewidth', line);
    set(gcf, 'Position', [100, 100, pix, pix/1.3801]);
    filenameTrace_ncm = strcat(folderPathMCMC, 'Trace_nc', num2str(m-numComponents+1));
    print(filenameTrace_ncm, '-dpng', '-r600');
end

% ---------- Trace plots of the topochemical atom indices by the Markov Chain algorithm ---------- %
for m = 2*numComponents:3*numComponents-1
    t = tiledlayout(1,1); t.Padding = 'none'; t.TileSpacing = 'none'; nexttile
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
    ylabelHandle = ylabel(ylabelName);
    set(xlabelHandle, 'FontSize', fontsize, 'FontName', fontname, 'Interpreter', 'latex');
    set(ylabelHandle, 'FontSize', fontsize, 'FontName', fontname, 'Interpreter', 'latex');
    set(0,'defaultfigurecolor',[1 1 1]);
    axesHandle = gca;
    set(axesHandle, 'FontSize', fontsize, 'FontName', fontname, 'Linewidth', line);
    set(gcf, 'Position', [100, 100, pix, pix/1.3801]);
    filenameTrace_ncm = strcat(folderPathMCMC, 'Trace_eta', num2str(m-2*numComponents+1));
    print(filenameTrace_ncm, '-dpng', '-r600');
end

% ---------- Plot the acceptance rate by the Markov Chain algorithm ---------- %
t = tiledlayout(1,1); t.Padding = 'none'; t.TileSpacing = 'none'; nexttile
figure('visible', 'off');
hold on;
box on;
colors = parula(N_chains); % Get a colormap with N_chains distinct colors
for i = 1:N_chains
    plot(AR(:,i), 'LineWidth', 0.75, 'Color',colors(i, :));
end
xlim([0 t_convergence]);
ylim([0 100]);
xlabelName = strcat('MCMC Iterations');
ylabelName = strcat('Acceptance Rate [\%]');
xlabelHandle = xlabel(xlabelName);
ylabelHandle = ylabel(ylabelName);
set(xlabelHandle, 'FontSize', fontsize, 'FontName', fontname, 'Interpreter', 'latex');
set(ylabelHandle, 'FontSize', fontsize, 'FontName', fontname, 'Interpreter', 'latex');
set(0,'defaultfigurecolor',[1 1 1]);
axesHandle = gca;
set(axesHandle, 'FontSize', fontsize, 'FontName', fontname, 'Linewidth', line);
set(gcf, 'Position', [100, 100, pix, pix/1.3801]);
filenameAcceptanceRate = strcat(folderPathMCMC, 'AcceptanceRate');
print(filenameAcceptanceRate, '-dpng', '-r600');

% ---------- Plot the R-hat statistic for each parameter by the Markov Chain algorithm ---------- %
t = tiledlayout(1,1); t.Padding = 'none'; t.TileSpacing = 'none'; nexttile
figure('visible', 'off');
hold on;
box on;
abscissa = linspace(2*t_burnin, t_convergence, t_convergence-2*t_burnin+1);
colors = parula(3*numComponents-1); % Get a colormap with N_chains distinct colors
for i = 1:3*numComponents-1
    plot(abscissa(1:2:end), R_hat(:,i), 'LineWidth', 0.75, 'Color',colors(i, :));
end
y_line = 1.2; % R-hat statistic threshold denoting convergence
plot(abscissa, y_line*ones(size(abscissa)), 'r', 'LineWidth', 1.5, 'LineStyle', '--');
xlim([2*t_burnin t_convergence]);
ylim([0 5]);
xlabelName = strcat('MCMC Iterations');
ylabelName = strcat('R-hat [-]');
xlabelHandle = xlabel(xlabelName);
ylabelHandle = ylabel(ylabelName);
set(xlabelHandle, 'FontSize', fontsize, 'FontName', fontname, 'Interpreter', 'latex');
set(ylabelHandle, 'FontSize', fontsize, 'FontName', fontname, 'Interpreter', 'latex');
set(0,'defaultfigurecolor',[1 1 1]);
axesHandle = gca;
set(axesHandle, 'FontSize', fontsize, 'FontName', fontname, 'Linewidth', line);
set(gcf, 'Position', [100, 100, pix, pix/1.3801]);
filenameRhat = strcat(folderPathMCMC, 'Rhat');
print(filenameRhat, '-dpng', '-r600');

% ---------- Check the performance of the MAP surrogate against experimental data ---------- %

% Indices of the surrogate components concerning the numbers of carbon atoms and topochemical atom indices
index_n_eta_MAP = zeros(1, numComponents-1);

for i = 1:numComponents
    index_n_eta_MAP(i) = findIndex_eta(classes{i}, nc_MAP(i), eta_B_star_MAP(i));
end

counter = 1;
lumped_variables = {};

for h = 1:numel(variable_names)

    if ~strcmp(variable_names{h}, 'molWeight') && ~strcmp(variable_names{h}, 'HC') && ~strcmp(variable_names{h}, 'DCN') && ~strcmp(variable_names{h}, 'flash') && ~strcmp(variable_names{h}, 'freezing') && ~strcmp(variable_names{h}, 'LHV')

        temperature_range = linspace(minT_array(h), maxT_array(h), 25);
        property_MAP = zeros(numel(temperature_range),1);

        if ~strcmp(variable_names{h}, 'distillation')

            for i = 1:numel(temperature_range)
                if strcmp(variable_names{h}, 'mu') || strcmp(variable_names{h}, 'nu')
                    property_MAP(i) = 1e+06*ThermophysicalProperties_LiquidMixture(variable_names{h}, temperature_range(i), x_MAP(1:numComponents-1), classes, index_n_eta_MAP, pressure_distillation);
                elseif strcmp(variable_names{h}, 'specificHeat')
                    property_MAP(i) = 1e-03*ThermophysicalProperties_LiquidMixture(variable_names{h}, temperature_range(i), x_MAP(1:numComponents-1), classes, index_n_eta_MAP, pressure_distillation);
                else
                    property_MAP(i) = ThermophysicalProperties_LiquidMixture(variable_names{h}, temperature_range(i), x_MAP(1:numComponents-1), classes, index_n_eta_MAP, pressure_distillation);
                end
            end

        else
            [volumeFraction_MAP, property_MAP] = DistillationCurve(x_MAP(1:numComponents-1), classes, index_n_eta_MAP, pressure_distillation);
        end

        percentiles_list = (0 + (1 - confidence_width)/2 ):0.01:(1 - (1 - confidence_width)/2 ); % confidence interval

        mean_model = zeros(numel(temperature_range), 1);
        percentiles_model = zeros(numel(temperature_range), size(percentiles_list,2));

        model_output_samples = zeros(size(chain_reshaped,1), numel(temperature_range));

        % Identify unique chain samples and their indices
        [unique_rows, ~, unique_idx] = unique(round(chain_reshaped, 6), 'rows');

        if ~strcmp(variable_names{h}, 'distillation')

            % Preallocate for storing results of unique chain samples
            unique_results = zeros(size(unique_rows, 1), numel(0:1:100));

            for i = 1:numel(temperature_range)
                parfor j = 1:size(unique_rows,1)
                    index_n_eta_samples = zeros(1,numComponents);
                    for k = 1:numComponents
                        index_n_eta_samples(k) = findIndex_eta(classes{k}, unique_rows(j, numComponents+k-1), unique_rows(j, 2*numComponents+k-1));
                    end
                    if strcmp(variable_names{h}, 'mu') || strcmp(variable_names{h}, 'nu')
                        unique_results(j,i) = 1e+06*ThermophysicalProperties_LiquidMixture(variable_names{h}, temperature_range(i), unique_rows(j, 1:numComponents-1), classes, index_n_eta_samples, pressure_distillation);
                    elseif strcmp(variable_names{h}, 'specificHeat')
                        unique_results(j,i) = 1e-03*ThermophysicalProperties_LiquidMixture(variable_names{h}, temperature_range(i), unique_rows(j, 1:numComponents-1), classes, index_n_eta_samples, pressure_distillation);
                    else
                        unique_results(j,i) = ThermophysicalProperties_LiquidMixture(variable_names{h}, temperature_range(i), unique_rows(j, 1:numComponents-1), classes, index_n_eta_samples, pressure_distillation);
                    end
                end
                % Map results to all samples
                for j = 1:size(model_output_samples,1)
                    model_output_samples(j, i) = unique_results(unique_idx(j), i);
                end
                mean_model(i) = mean(model_output_samples(:,i));
                for l = 1:size(percentiles_list,2)
                    percentiles_model(i,l) = quantile(model_output_samples(:,i), percentiles_list(l));
                end
            end

        else
            volFrac_interp = linspace(1, 99, 99);
            volFrac_interp = [0.1, volFrac_interp, 99.9];
            mean_model = zeros(numel(volFrac_interp), 1);
            percentiles_model = zeros(numel(volFrac_interp), size(percentiles_list,2));
            model_output_samples = zeros(size(chain_reshaped,1), numel(volFrac_interp));
            % Preallocate for storing results of unique chain samples
            unique_results = zeros(size(unique_rows, 1), numel(volFrac_interp));
            parfor j = 1:size(unique_rows,1)
                index_n_eta_samples = zeros(1,numComponents);
                for k = 1:numComponents
                    index_n_eta_samples(k) = findIndex_eta(classes{k}, unique_rows(j, numComponents+k-1), unique_rows(j, 2*numComponents+k-1));
                end
                [volumeFraction, property] = DistillationCurve(unique_rows(j, 1:numComponents-1), classes, index_n_eta_samples, pressure_distillation);
                T_interpolated = interp1(volumeFraction, property, volFrac_interp, 'linear');
                unique_results(j,:) = T_interpolated;
            end
            % Map results to all samples
            for j = 1:size(model_output_samples,1)
                model_output_samples(j, :) = unique_results(unique_idx(j), :);
            end
            for i = 1:numel(volFrac_interp)
                mean_model(i) = mean(model_output_samples(:,i));
                for l = 1:size(percentiles_list,2)
                    percentiles_model(i,l) = quantile(model_output_samples(:,i), percentiles_list(l));
                end
            end
        end

    elseif strcmp(variable_names{h}, 'molWeight')

        Wi_MAP = zeros(1, numComponents);
        for i = 1:numComponents
            Wi_MAP(i) = classes{i}(index_n_eta_MAP(i)).molWeight;
        end

        property_MAP = molWeight(x_MAP, Wi_MAP);

        percentiles_list = (0 + (1 - confidence_width)/2 ):0.01:(1 - (1 - confidence_width)/2 );
        percentiles_model = zeros(1, size(percentiles_list,2));
        model_output_samples = zeros(size(chain_reshaped,1), 1);

        for j = 1:size(chain_reshaped,1)
            index_n_eta_samples = zeros(1,numComponents);
            for k = 1:numComponents
                index_n_eta_samples(k) = findIndex_eta(classes{k}, chain_reshaped(j, numComponents+k-1), chain_reshaped(j, 2*numComponents+k-1));
            end
            if j > 1 && isequal(chain_reshaped(j,:),chain_reshaped(j-1,:))
                model_output_samples(j) = model_output_samples(j-1);
            else
                Wl = zeros(1, numComponents);
                for l = 1:numComponents
                    Wl(l) = classes{l}(index_n_eta_samples(l)).molWeight;
                end
                model_output_samples(j) = molWeight([chain_reshaped(j, 1:numComponents-1) 1-sum(chain_reshaped(j, 1:numComponents-1))], Wl);
            end
        end
        mean_model = mean(model_output_samples);
        for l = 1:size(percentiles_list,2)
            percentiles_model(l) = quantile(model_output_samples, percentiles_list(l));
        end

    elseif strcmp(variable_names{h}, 'HC')

        % Compute the molecular formula of the MAP surrogate mixture
        nCarbon_MAP = round(sum(x_MAP.*nc_MAP),2); % number of carbon atoms denoting the molecular formula
        nHydrogen_MAP = 0.0; % number of hydrogen atoms denoting the molecular formula
        for l = 1:numComponents
            if strcmp(families(l), 'nparaffins') || strcmp(families(l), 'isoparaffins')
                nHydrogen_MAP = nHydrogen_MAP + x_MAP(l)*(2*nc_MAP(l)+2);
            elseif strcmp(families(l), 'cycloparaffins')
                nHydrogen_MAP = nHydrogen_MAP + x_MAP(l)*2*nc_MAP(l);
            elseif strcmp(families(l), 'dicycloparaffins')
                nHydrogen_MAP = nHydrogen_MAP + x_MAP(l)*(2*nc_MAP(l)-2);
            elseif strcmp(families(l), 'alkylbenzenes')
                nHydrogen_MAP = nHydrogen_MAP + x_MAP(l)*(2*nc_MAP(l)-6);
            elseif strcmp(families(l), 'alkylnaphtalenes')
                nHydrogen_MAP = nHydrogen_MAP + x_MAP(l)*(2*nc_MAP(l)-12);
            elseif strcmp(families(l), 'cycloaromatics')
                nHydrogen_MAP = nHydrogen_MAP + x_MAP(l)*(2*nc_MAP(l)-8);
            end
        end

        nHydrogen_MAP = round(nHydrogen_MAP, 2);

        % Hydrogen-to-carbon ratio (H/C) of the MAP surrogate mixture
        property_MAP = round(nHydrogen_MAP/nCarbon_MAP, 3);

        percentiles_list = (0 + (1 - confidence_width)/2 ):0.01:(1 - (1 - confidence_width)/2 );
        percentiles_model = zeros(1, size(percentiles_list,2));
        model_output_samples = zeros(size(chain_reshaped,1), 1);

        for j = 1:size(chain_reshaped,1)

            if j > 1 && isequal(chain_reshaped(j,:),chain_reshaped(j-1,:))
                model_output_samples(j) = model_output_samples(j-1);
            else
                molFrac_j = [chain_reshaped(j,1:numComponents-1) 1-sum(chain_reshaped(j,1:numComponents-1))];
                nC_j = chain_reshaped(j,numComponents:2*numComponents-1);
                nCarbon = round(sum(molFrac_j.*chain_reshaped(j,numComponents:2*numComponents-1)),2);
                nHydrogen = 0.0;
                for l = 1:numComponents
                    if strcmp(families(l), 'nparaffins') || strcmp(families(l), 'isoparaffins')
                        nHydrogen = nHydrogen + molFrac_j(l)*(2*nC_j(l)+2);
                    elseif strcmp(families(l), 'cycloparaffins')
                        nHydrogen = nHydrogen + molFrac_j(l)*2*nC_j(l);
                    elseif strcmp(families(l), 'dicycloparaffins')
                        nHydrogen = nHydrogen + molFrac_j(l)*(2*nC_j(l)-2);
                    elseif strcmp(families(l), 'alkylbenzenes')
                        nHydrogen = nHydrogen + molFrac_j(l)*(2*nC_j(l)-6);
                    elseif strcmp(families(l), 'alkylnaphtalenes')
                        nHydrogen = nHydrogen + molFrac_j(l)*(2*nC_j(l)-12);
                    elseif strcmp(families(l), 'cycloaromatics')
                        nHydrogen = nHydrogen + molFrac_j(l)*(2*nC_j(l)-8);
                    end
                end
                nHydrogen = round(nHydrogen, 2);
                model_output_samples(j) = round(nHydrogen/nCarbon, 3);
            end

        end

        mean_model = mean(model_output_samples);
        for l = 1:size(percentiles_list,2)
            percentiles_model(l) = quantile(model_output_samples, percentiles_list(l));
        end

    elseif strcmp(variable_names{h}, 'DCN')

        % +++ Calculate the derived cetane number (DCN) of the MAP surrogate mixture
        rho_i = zeros(1,numComponents);
        numerator = zeros(1,numComponents);
        DCN_i = zeros(1,numComponents);

        Wi_MAP = zeros(1, numComponents);
        for i = 1:numComponents
            Wi_MAP(i) = classes{i}(index_n_eta_MAP(i)).molWeight;
        end

        y_MAP = mol2mass(x_MAP, Wi_MAP);

        for i = 1:numComponents

            rho_i(i) = ThermophysicalProperties_SingleLiquid('rho', 300, classes{i}(index_n_eta_MAP(i)).molWeight, classes{i}(index_n_eta_MAP(i)).Tc, classes{i}(index_n_eta_MAP(i)).coeffRho(1), classes{i}(index_n_eta_MAP(i)).coeffRho(2), classes{i}(index_n_eta_MAP(i)).coeffRho(3), classes{i}(index_n_eta_MAP(i)).coeffMu(1), classes{i}(index_n_eta_MAP(i)).coeffMu(2), classes{i}(index_n_eta_MAP(i)).coeffMu(3), classes{i}(index_n_eta_MAP(i)).coeffMu(4), classes{i}(index_n_eta_MAP(i)).coeffK(1), classes{i}(index_n_eta_MAP(i)).coeffK(2), classes{i}(index_n_eta_MAP(i)).coeffK(3), classes{i}(index_n_eta_MAP(i)).coeffCl(1), classes{i}(index_n_eta_MAP(i)).coeffCl(2), classes{i}(index_n_eta_MAP(i)).coeffCl(3), classes{i}(index_n_eta_MAP(i)).coeffCl(4), classes{i}(index_n_eta_MAP(i)).coeffHv(1), classes{i}(index_n_eta_MAP(i)).coeffHv(2), classes{i}(index_n_eta_MAP(i)).coeffPsat(1), classes{i}(index_n_eta_MAP(i)).coeffPsat(2), classes{i}(index_n_eta_MAP(i)).coeffPsat(3), classes{i}(index_n_eta_MAP(i)).coeffPsat(4), classes{i}(index_n_eta_MAP(i)).coeffPsat(5), classes{i}(index_n_eta_MAP(i)).coeffSigma(1), classes{i}(index_n_eta_MAP(i)).coeffSigma(2), 101325);
            numerator(i) = x_MAP(i)*classes{i}(index_n_eta_MAP(i)).molWeight;
            DCN_i(i) = classes{i}(index_n_eta_MAP(i)).DCN;

        end

        % Compute the MAP surrogate mixture density
        rho_mixture = sum(numerator, 2)./sum(numerator./ rho_i, 2);

        % Compute the MAP surrogate mixture volume fractions
        Vi = rho_mixture*y_MAP'./rho_i';

        property_MAP = sum(Vi.*DCN_i');

        percentiles_list = (0 + (1 - confidence_width)/2 ):0.01:(1 - (1 - confidence_width)/2 );
        percentiles_model = zeros(1, size(percentiles_list,2));
        model_output_samples = zeros(size(chain_reshaped,1), 1);

        for j = 1:size(chain_reshaped,1)
            index_n_eta_samples = zeros(1,numComponents);
            for k = 1:numComponents
                index_n_eta_samples(k) = findIndex_eta(classes{k}, chain_reshaped(j, numComponents+k-1), chain_reshaped(j, 2*numComponents+k-1));
            end
            if j > 1 && isequal(chain_reshaped(j,:),chain_reshaped(j-1,:))
                model_output_samples(j) = model_output_samples(j-1);
            else
                rho_l = zeros(1,numComponents);
                numerator = zeros(1,numComponents);
                DCN_l = zeros(1,numComponents);

                x_l = [chain_reshaped(j,1:numComponents-1) 1-sum(chain_reshaped(j,1:numComponents-1))];
                Wl = zeros(1, numComponents);

                for l = 1:numComponents

                    rho_l(l) = ThermophysicalProperties_SingleLiquid('rho', 300, classes{l}(index_n_eta_samples(l)).molWeight, classes{l}(index_n_eta_samples(l)).Tc, classes{l}(index_n_eta_samples(l)).coeffRho(1), classes{l}(index_n_eta_samples(l)).coeffRho(2), classes{l}(index_n_eta_samples(l)).coeffRho(3), classes{l}(index_n_eta_samples(l)).coeffMu(1), classes{l}(index_n_eta_samples(l)).coeffMu(2), classes{l}(index_n_eta_samples(l)).coeffMu(3), classes{l}(index_n_eta_samples(l)).coeffMu(4), classes{l}(index_n_eta_samples(l)).coeffK(1), classes{l}(index_n_eta_samples(l)).coeffK(2), classes{l}(index_n_eta_samples(l)).coeffK(3), classes{l}(index_n_eta_samples(l)).coeffCl(1), classes{l}(index_n_eta_samples(l)).coeffCl(2), classes{l}(index_n_eta_samples(l)).coeffCl(3), classes{l}(index_n_eta_samples(l)).coeffCl(4), classes{l}(index_n_eta_samples(l)).coeffHv(1), classes{l}(index_n_eta_samples(l)).coeffHv(2), classes{l}(index_n_eta_samples(l)).coeffPsat(1), classes{l}(index_n_eta_samples(l)).coeffPsat(2), classes{l}(index_n_eta_samples(l)).coeffPsat(3), classes{l}(index_n_eta_samples(l)).coeffPsat(4), classes{l}(index_n_eta_samples(l)).coeffPsat(5), classes{l}(index_n_eta_samples(l)).coeffSigma(1), classes{l}(index_n_eta_samples(l)).coeffSigma(2), 101325);
                    Wl(l) = classes{l}(index_n_eta_samples(l)).molWeight;
                    numerator(l) = x_l(l)*classes{l}(index_n_eta_samples(l)).molWeight;
                    DCN_l(l) = classes{l}(index_n_eta_samples(l)).DCN;

                end

                y_l = mol2mass(x_l, Wl);

                rho_mixture = sum(numerator, 2)./sum(numerator./ rho_l, 2);
                Vl = rho_mixture*y_l'./rho_l';
                model_output_samples(j) = sum(Vl.*DCN_l');
            end
        end

        mean_model = mean(model_output_samples);
        for l = 1:size(percentiles_list,2)
            percentiles_model(l) = quantile(model_output_samples, percentiles_list(l));
        end

    elseif strcmp(variable_names{h}, 'flash')

        % +++ Calculate the flash point in [K] of the MAP surrogate mixture
        flash_i = zeros(1,numComponents);
        BI_flash_i = zeros(1,numComponents);

        Wi_MAP = zeros(1, numComponents);
        for i = 1:numComponents
            Wi_MAP(i) = classes{i}(index_n_eta_MAP(i)).molWeight;
        end

        y_MAP = mol2mass(x_MAP, Wi_MAP);

        for i = 1:numComponents

            flash_i(i) = classes{i}(index_n_eta_MAP(i)).Tf; % flash point of the i-th MAP surrogate in [K]
            flash_F = (flash_i(i)-273.15)*9./5+32; % flash point of the i-th MAP surrogate in [°F]
            BI_flash_i(i) = 51708*exp((log(flash_F)-2.6287)^2/(-0.91725));

        end

        BI_blend = sum(y_MAP.*BI_flash_i);

        property_MAP_F = exp(((-0.91725)*log(BI_blend/51708))^0.5+2.6287); % MAP surrogate mixture flash point in [°F]
        property_MAP = (property_MAP_F-32)*5./9+273.15; % MAP surrogate mixture flash point in [K]

        percentiles_list = (0 + (1 - confidence_width)/2 ):0.01:(1 - (1 - confidence_width)/2 );
        percentiles_model = zeros(1, size(percentiles_list,2));
        model_output_samples = zeros(size(chain_reshaped,1), 1);

        for j = 1:size(chain_reshaped,1)
            index_n_eta_samples = zeros(1,numComponents);
            for k = 1:numComponents
                index_n_eta_samples(k) = findIndex_eta(classes{k}, chain_reshaped(j, numComponents+k-1), chain_reshaped(j, 2*numComponents+k-1));
            end
            if j > 1 && isequal(chain_reshaped(j,:),chain_reshaped(j-1,:))
                model_output_samples(j) = model_output_samples(j-1);
            else
                flash_l = zeros(1,numComponents);
                BI_flash_l = zeros(1,numComponents);

                x_l = [chain_reshaped(j,1:numComponents-1) 1-sum(chain_reshaped(j,1:numComponents-1))];
                Wl = zeros(1, numComponents);

                for l = 1:numComponents

                    Wl(l) = classes{l}(index_n_eta_samples(l)).molWeight;
                    flash_l(l) = classes{l}(index_n_eta_samples(l)).Tf;
                    flash_F = (flash_l(l)-273.15)*9./5+32;
                    BI_flash_l(l) = 51708*exp((log(flash_F)-2.6287)^2/(-0.91725));

                end

                y_l = mol2mass(x_l, Wl);

                BI_blend = sum(y_l.*BI_flash_l);

                model_output_samples_F = exp(((-0.91725)*log(BI_blend/51708))^0.5+2.6287); % mixture flash point [°F]
                model_output_samples(j) = (model_output_samples_F-32)*5./9+273.15; % mixture flash point [K]
            end
        end

        mean_model = mean(model_output_samples);
        for l = 1:size(percentiles_list,2)
            percentiles_model(l) = quantile(model_output_samples, percentiles_list(l));
        end

    elseif strcmp(variable_names{h}, 'freezing')

        % +++ Calculate the freezing point in [K] of the MAP surrogate mixture
        freezing_i = zeros(1,numComponents);
        BI_freezing_i = zeros(1,numComponents);

        Wi_MAP = zeros(1, numComponents);
        rho_i = zeros(1,numComponents);
        numerator = zeros(1,numComponents);
        for i = 1:numComponents
            rho_i(i) = ThermophysicalProperties_SingleLiquid('rho', 300, classes{i}(index_n_eta_MAP(i)).molWeight, classes{i}(index_n_eta_MAP(i)).Tc, classes{i}(index_n_eta_MAP(i)).coeffRho(1), classes{i}(index_n_eta_MAP(i)).coeffRho(2), classes{i}(index_n_eta_MAP(i)).coeffRho(3), classes{i}(index_n_eta_MAP(i)).coeffMu(1), classes{i}(index_n_eta_MAP(i)).coeffMu(2), classes{i}(index_n_eta_MAP(i)).coeffMu(3), classes{i}(index_n_eta_MAP(i)).coeffMu(4), classes{i}(index_n_eta_MAP(i)).coeffK(1), classes{i}(index_n_eta_MAP(i)).coeffK(2), classes{i}(index_n_eta_MAP(i)).coeffK(3), classes{i}(index_n_eta_MAP(i)).coeffCl(1), classes{i}(index_n_eta_MAP(i)).coeffCl(2), classes{i}(index_n_eta_MAP(i)).coeffCl(3), classes{i}(index_n_eta_MAP(i)).coeffCl(4), classes{i}(index_n_eta_MAP(i)).coeffHv(1), classes{i}(index_n_eta_MAP(i)).coeffHv(2), classes{i}(index_n_eta_MAP(i)).coeffPsat(1), classes{i}(index_n_eta_MAP(i)).coeffPsat(2), classes{i}(index_n_eta_MAP(i)).coeffPsat(3), classes{i}(index_n_eta_MAP(i)).coeffPsat(4), classes{i}(index_n_eta_MAP(i)).coeffPsat(5), classes{i}(index_n_eta_MAP(i)).coeffSigma(1), classes{i}(index_n_eta_MAP(i)).coeffSigma(2), 101325);
            Wi_MAP(i) = classes{i}(index_n_eta_MAP(i)).molWeight;
            numerator(i) = x_MAP(i)*classes{i}(index_n_eta_MAP(i)).molWeight;
        end

        y_MAP = mol2mass(x_MAP, Wi_MAP);

        for i = 1:numComponents

            freezing_i(i) = classes{i}(index_n_eta_MAP(i)).Tfz; % freezing point of the i-th MAP surrogate in [K]
            BI_freezing_i(i) = freezing_i(i)^(1/0.05);

        end

        rho_MAP = sum(numerator, 2)./sum(numerator./ rho_i, 2);
        Vi_MAP = rho_MAP*y_MAP./rho_i; % volume fractions

        BI_blend = sum(Vi_MAP.*BI_freezing_i);

        property_MAP = BI_blend^0.05; % MAP surrogate mixture freezing point in [K]

        percentiles_list = (0 + (1 - confidence_width)/2 ):0.01:(1 - (1 - confidence_width)/2 );
        percentiles_model = zeros(1, size(percentiles_list,2));
        model_output_samples = zeros(size(chain_reshaped,1), 1);

        for j = 1:size(chain_reshaped,1)
            index_n_eta_samples = zeros(1,numComponents);
            for k = 1:numComponents
                index_n_eta_samples(k) = findIndex_eta(classes{k}, chain_reshaped(j, numComponents+k-1), chain_reshaped(j, 2*numComponents+k-1));
            end
            if j > 1 && isequal(chain_reshaped(j,:),chain_reshaped(j-1,:))
                model_output_samples(j) = model_output_samples(j-1);
            else
                freezing_l = zeros(1,numComponents);
                BI_freezing_l = zeros(1,numComponents);

                x_l = [chain_reshaped(j,1:numComponents-1) 1-sum(chain_reshaped(j,1:numComponents-1))];
                Wl = zeros(1, numComponents);
                rho_l = zeros(1,numComponents);
                numerator = zeros(1,numComponents);

                for l = 1:numComponents

                    Wl(l) = classes{l}(index_n_eta_samples(l)).molWeight;
                    rho_l(l) = ThermophysicalProperties_SingleLiquid('rho', 300, classes{l}(index_n_eta_samples(l)).molWeight, classes{l}(index_n_eta_samples(l)).Tc, classes{l}(index_n_eta_samples(l)).coeffRho(1), classes{l}(index_n_eta_samples(l)).coeffRho(2), classes{l}(index_n_eta_samples(l)).coeffRho(3), classes{l}(index_n_eta_samples(l)).coeffMu(1), classes{l}(index_n_eta_samples(l)).coeffMu(2), classes{l}(index_n_eta_samples(l)).coeffMu(3), classes{l}(index_n_eta_samples(l)).coeffMu(4), classes{l}(index_n_eta_samples(l)).coeffK(1), classes{l}(index_n_eta_samples(l)).coeffK(2), classes{l}(index_n_eta_samples(l)).coeffK(3), classes{l}(index_n_eta_samples(l)).coeffCl(1), classes{l}(index_n_eta_samples(l)).coeffCl(2), classes{l}(index_n_eta_samples(l)).coeffCl(3), classes{l}(index_n_eta_samples(l)).coeffCl(4), classes{l}(index_n_eta_samples(l)).coeffHv(1), classes{l}(index_n_eta_samples(l)).coeffHv(2), classes{l}(index_n_eta_samples(l)).coeffPsat(1), classes{l}(index_n_eta_samples(l)).coeffPsat(2), classes{l}(index_n_eta_samples(l)).coeffPsat(3), classes{l}(index_n_eta_samples(l)).coeffPsat(4), classes{l}(index_n_eta_samples(l)).coeffPsat(5), classes{l}(index_n_eta_samples(l)).coeffSigma(1), classes{l}(index_n_eta_samples(l)).coeffSigma(2), 101325);
                    freezing_l(l) = classes{l}(index_n_eta_samples(l)).Tfz;
                    BI_freezing_l(l) = freezing_l(l)^(1/0.05);
                    numerator(l) = x_l(l)*classes{l}(index_n_eta_samples(l)).molWeight;

                end

                y_l = mol2mass(x_l, Wl);
                rho_mixture = sum(numerator, 2)./sum(numerator./ rho_l, 2);
                Vl = rho_mixture*y_l./rho_l; % volume fractions

                BI_blend = sum(Vl.*BI_freezing_l);

                model_output_samples(j) = BI_blend^0.05; % mixture freezing point [K]
            end
        end

        mean_model = mean(model_output_samples);
        for l = 1:size(percentiles_list,2)
            percentiles_model(l) = quantile(model_output_samples, percentiles_list(l));
        end

    elseif strcmp(variable_names{h}, 'LHV')

        % +++ Calculate the LHV in [MJ/kg] of the MAP surrogate mixture
        LHV_i = zeros(1,numComponents);

        Wi_MAP = zeros(1, numComponents);
        for i = 1:numComponents
            Wi_MAP(i) = classes{i}(index_n_eta_MAP(i)).molWeight;
        end

        y_MAP = mol2mass(x_MAP, Wi_MAP);

        for i = 1:numComponents

            LHV_i(i) = classes{i}(index_n_eta_MAP(i)).Hc; % LHV of the i-th MAP surrogate in [MJ/kg]

        end

        property_MAP = sum(y_MAP.*LHV_i); % MAP LHV in [MJ/kg]

        percentiles_list = (0 + (1 - confidence_width)/2 ):0.01:(1 - (1 - confidence_width)/2 );
        percentiles_model = zeros(1, size(percentiles_list,2));
        model_output_samples = zeros(size(chain_reshaped,1), 1);

        for j = 1:size(chain_reshaped,1)
            index_n_eta_samples = zeros(1,numComponents);
            for k = 1:numComponents
                index_n_eta_samples(k) = findIndex_eta(classes{k}, chain_reshaped(j, numComponents+k-1), chain_reshaped(j, 2*numComponents+k-1));
            end
            if j > 1 && isequal(chain_reshaped(j,:),chain_reshaped(j-1,:))
                model_output_samples(j) = model_output_samples(j-1);
            else
                LHV_l = zeros(1,numComponents);

                x_l = [chain_reshaped(j,1:numComponents-1) 1-sum(chain_reshaped(j,1:numComponents-1))];
                Wl = zeros(1, numComponents);

                for l = 1:numComponents

                    Wl(l) = classes{l}(index_n_eta_samples(l)).molWeight;
                    LHV_l(l) = classes{l}(index_n_eta_samples(l)).Hc;

                end

                y_l = mol2mass(x_l, Wl);

                model_output_samples(j) = sum(y_l.*LHV_l); % mixture LHV [MJ/kg]
            end
        end

        mean_model = mean(model_output_samples);
        for l = 1:size(percentiles_list,2)
            percentiles_model(l) = quantile(model_output_samples, percentiles_list(l));
        end

    end

    if ~strcmp(variable_names{h}, 'molWeight') && ~strcmp(variable_names{h}, 'HC') && ~strcmp(variable_names{h}, 'DCN') && ~strcmp(variable_names{h}, 'flash') && ~strcmp(variable_names{h}, 'freezing') && ~strcmp(variable_names{h}, 'LHV')

        t = tiledlayout(1,1); t.Padding = 'none'; t.TileSpacing = 'none'; nexttile
        figure('visible', 'off');
        xlabelName = 'Temperature [K]';
        if strcmp(variable_names{h}, 'rho')
            ylabelName = 'Density [kg/m$^3$]';
        end
        if strcmp(variable_names{h}, 'mu')
            ylabelName = 'Dynamic Viscosity [$\mu$Pa s]';
        end
        if strcmp(variable_names{h}, 'nu')
            ylabelName = 'Kinematic Viscosity [mm$^2$/s]';
        end
        if strcmp(variable_names{h}, 'kappa')
            ylabelName = 'Thermal Conductivity [W/m/K]';
        end
        if strcmp(variable_names{h}, 'specificHeat')
            ylabelName = 'Specific Heat Capacity [kJ/kg/K]';
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
        set(xlabelHandle, 'FontSize', fontsize, 'FontName', fontname, 'Interpreter', 'latex');
        set(ylabelHandle, 'FontSize', fontsize, 'FontName', fontname, 'Interpreter', 'latex');
        set(0,'defaultfigurecolor',[1 1 1]);
        axesHandle = gca;
        set(axesHandle, 'FontSize', fontsize, 'FontName', fontname, 'Linewidth', line);
        set(gcf, 'Position', [100, 100, pix, pix/1.3801]);
        if strcmp(band_percentiles, 'True')
            cmap = slanCM('bwr',size(percentiles_model,2));
        end
        hold on;
        if strcmp(variable_names{h}, 'rho') || strcmp(variable_names{h}, 'mu') || strcmp(variable_names{h}, 'nu') || strcmp(variable_names{h}, 'kappa') || strcmp(variable_names{h}, 'specificHeat') || strcmp(variable_names{h}, 'latentHeat') || strcmp(variable_names{h}, 'vaporPressure') || strcmp(variable_names{h}, 'sigma')
            if strcmp(variable_names{h}, 'specificHeat')
                ylim([0.98*min(min(percentiles_model(:,1)),1e-03*min(fullData{h}(:,2))) 1.02*max(max(percentiles_model(:,end)),1e-03*max(fullData{h}(:,2)))]);
            else
                ylim([0.98*min(min(percentiles_model(:,1)),min(fullData{h}(:,2))) 1.02*max(max(percentiles_model(:,end)),max(fullData{h}(:,2)))]);
            end
            xlim([0.95*temperature_range(1) 1.05*temperature_range(end)]);
            legendAdded = false;  % Flag to check if the legend is already added
            for i = 1:length(temperature_range)-1
                for j = 1:numel(percentiles_model(i,:))-1
                    if strcmp(band_percentiles, 'True')
                        fill([temperature_range(i) temperature_range(i) temperature_range(i+1) temperature_range(i+1)], ...
                            [percentiles_model(i,j) percentiles_model(i,j+1) percentiles_model(i+1,j+1) percentiles_model(i+1,j)], ...
                            cmap(j, :), 'EdgeColor', 'none', 'HandleVisibility', 'off');
                    else
                        if ~legendAdded
                            fill([temperature_range(i) temperature_range(i) temperature_range(i+1) temperature_range(i+1)], ...
                                [percentiles_model(i,j) percentiles_model(i,j+1) percentiles_model(i+1,j+1) percentiles_model(i+1,j)], ...
                                blue, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', '90% Confidence Interval');
                            legendAdded = true;  % Set flag to true after adding the legend
                        else
                            fill([temperature_range(i) temperature_range(i) temperature_range(i+1) temperature_range(i+1)], ...
                                [percentiles_model(i,j) percentiles_model(i,j+1) percentiles_model(i+1,j+1) percentiles_model(i+1,j)], ...
                                blue, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                        end
                    end
                end
            end
            plot(temperature_range, property_MAP, 'LineWidth', 2, 'Color', 'black', 'DisplayName', 'MAP Surrogate');
            if strcmp(band_percentiles, 'True')
                colormap(cmap);
                hcb = colorbar;
                ylabel(hcb, 'Percentile [\%]', 'FontSize', fontsize, 'FontName', fontname, 'Interpreter', 'latex');
                clim([0 1]);
                hcb.Ticks = [0 0.2895 0.5 0.7105 1];
                hcb.TickLabels = {'2.5', '30', '50', '70', '97.5'};
            end
        elseif strcmp(variable_names{h}, 'distillation')
            xlim([volumeFraction_MAP(1) volumeFraction_MAP(end)]);
            ylim([0.99*min(min(percentiles_model(:,1)),min(fullData{h}(:,2))) 1.01*max(max(percentiles_model(:,end)),max(fullData{h}(:,2)))]);
            volumeFraction = linspace(1, 99, 99);
            volumeFraction = [0.1, volumeFraction, 99.9];
            legendAdded = false;  % Flag to check if the legend is already added
            for i = 1:length(volumeFraction)-1
                for j = 1:numel(percentiles_model(i,:))-1
                    if strcmp(band_percentiles, 'True')
                        fill([volumeFraction(i) volumeFraction(i) volumeFraction(i+1) volumeFraction(i+1)], [percentiles_model(i,j) percentiles_model(i,j+1) percentiles_model(i+1,j+1) percentiles_model(i+1,j)], cmap(j, :), 'EdgeColor', 'none', 'HandleVisibility', 'off');
                    else
                        if ~legendAdded
                            fill([volumeFraction(i) volumeFraction(i) volumeFraction(i+1) volumeFraction(i+1)], [percentiles_model(i,j) percentiles_model(i,j+1) percentiles_model(i+1,j+1) percentiles_model(i+1,j)], blue, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', '90% Confidence Interval');
                            legendAdded = true;  % Set flag to true after adding the legend
                        else
                            fill([volumeFraction(i) volumeFraction(i) volumeFraction(i+1) volumeFraction(i+1)], [percentiles_model(i,j) percentiles_model(i,j+1) percentiles_model(i+1,j+1) percentiles_model(i+1,j)], blue, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                        end
                    end
                end
            end
            plot(volumeFraction_MAP, property_MAP, 'LineWidth', 2, 'Color','black', 'DisplayName', 'MAP Surrogate');
            if strcmp(band_percentiles, 'True')
                colormap(cmap);
                hcb = colorbar;
                ylabel(hcb, 'Percentile [\%]', 'FontSize', fontsize, 'FontName', fontname, 'Interpreter', 'latex');
                clim([0 1]);
                hcb.Ticks = [0 0.277 0.5 0.722 1];
                hcb.TickLabels = {'5', '30', '50', '70', '95'};
            end
        end
        if strcmp(variable_names{h}, 'mu') || strcmp(variable_names{h}, 'nu')
            plot(fullData{h}(:,1), 1e+06*fullData{h}(:,2), 'd', 'MarkerSize', 8, 'MarkerFace', gray, 'MarkerEdge', 'k', 'LineWidth', 1, 'DisplayName', fuel_name);
        elseif strcmp(variable_names{h}, 'specificHeat')
            plot(fullData{h}(:,1), 1e-03*fullData{h}(:,2), 'd', 'MarkerSize', 8, 'MarkerFaceColor', gray, 'MarkerEdge', 'k', 'LineWidth', 1, 'DisplayName', fuel_name);
        elseif strcmp(variable_names{h}, 'distillation')
            xlabelName = 'Recovered Volume Fraction [\%]';
            xlabelHandle = xlabel(xlabelName);
            set(xlabelHandle, 'FontSize', fontsize, 'FontName', fontname, 'Interpreter', 'latex');
            plot(fullData{h}(:,1), fullData{h}(:,2), 'd', 'MarkerSize', 8, 'MarkerFace', gray, 'MarkerEdge', 'k', 'LineWidth', 1, 'DisplayName', fuel_name);
        else
            plot(fullData{h}(:,1), fullData{h}(:,2), 'd', 'MarkerSize', 8, 'MarkerFaceColor', gray, 'MarkerEdge', 'k', 'LineWidth', 1, 'DisplayName', fuel_name);
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


    else

        % Use kernel density estimation to estimate the PDF of the lumped variable under consideration
        [f_x, x_values] = ksdensity(model_output_samples);

        % Plot the marginal PDF
        t = tiledlayout(1,1); t.Padding = 'none'; t.TileSpacing = 'none'; nexttile
        figure('visible', 'off'); box on;
        plot(x_values, f_x, 'LineWidth', 2, 'Color',blue);
        if strcmp(variable_names{h}, 'molWeight')
            xlabelName = 'MW [g/mol]';
        elseif strcmp(variable_names{h}, 'HC')
            xlabelName = 'HC [-]';
        elseif strcmp(variable_names{h}, 'DCN')
            xlabelName = 'DCN [-]';
        elseif strcmp(variable_names{h}, 'flash')
            xlabelName = 'Flash Point [K]';
        elseif strcmp(variable_names{h}, 'freezing')
            xlabelName = 'Freezing Point [K]';
        elseif strcmp(variable_names{h}, 'LHV')
            xlabelName = 'Lower Heating Value [MJ/kg]';
        end
        ylabelName = 'Probability Density';
        xlabelHandle = xlabel(xlabelName);
        ylabelHandle = ylabel(ylabelName);
        set(xlabelHandle, 'FontSize', fontsize, 'FontName', fontname, 'Interpreter', 'latex');
        set(ylabelHandle, 'FontSize', fontsize, 'FontName', fontname, 'Interpreter', 'latex');
        set(0,'defaultfigurecolor',[1 1 1]);
        axesHandle = gca;
        set(axesHandle, 'FontSize', fontsize, 'FontName', fontname, 'Linewidth', line);
        set(gcf, 'Position', [100, 100, pix, pix/1.3801]);
        if strcmp(variable_names{h}, 'DCN')
            xlim([0 100]);
        else
            % Set the x-axis limit dynamically
            xlim([min(x_values), max(x_values)]);
        end
        folderPathPropertyCheck = './Figures/PropertyCheck_MAP/';
        % Create the folder if it doesn't exist
        if ~isfolder(folderPathPropertyCheck)
            mkdir(folderPathPropertyCheck);
        end
        filenamePropertyCheck = strcat(folderPathPropertyCheck, variable_names{h});
        print(filenamePropertyCheck, '-dpng', '-r600');

    end

    % ---------- Sensitivity analysis through grouped Sobol' first-order indices (3 groups: molar fractions, numbers of carbon atoms, branching indices) ---------- %

    if ~strcmp(variable_names{h}, 'molWeight') && ~strcmp(variable_names{h}, 'HC') && ~strcmp(variable_names{h}, 'DCN') && ~strcmp(variable_names{h}, 'flash') && ~strcmp(variable_names{h}, 'freezing') && ~strcmp(variable_names{h}, 'LHV')

        if ~strcmp(variable_names{h}, 'distillation')
            sobol_idx = sobol(classes, numComponents, temperature_range, chain_reshaped, model_output_samples, variable_names{h}, pressure_distillation);
        else
            sobol_idx = sobol_distillation(classes, numComponents, chain_reshaped, model_output_samples, pressure_distillation);
        end

        % Subset the temperature_range or recovered volume fraction and sobol_idx arrays
        if ~strcmp(variable_names{h}, 'distillation')
            selected_indices = [1, 7, 13, 19, 25];
            selected_temperature_range = temperature_range(selected_indices);
        else
            selected_indices = [2, 26, 51, 76, 100];
        end
        selected_sobol_idx = sobol_idx(:, selected_indices);

        t = tiledlayout(1,1); t.Padding = 'none'; t.TileSpacing = 'none'; nexttile
        figure('visible', 'off');
        barWidth = 1;
        if ~strcmp(variable_names{h}, 'distillation')
            b = bar(round(selected_temperature_range,0), [selected_sobol_idx(1,:); selected_sobol_idx(2,:); selected_sobol_idx(3,:)]', 'grouped', 'BarWidth', barWidth);
        else
            b = bar([1 25 50 75 99], [selected_sobol_idx(1,:); selected_sobol_idx(2,:); selected_sobol_idx(3,:)]', 'grouped', 'BarWidth', barWidth);
        end
        % Set colors for bars
        bar_colors = [green; blue; orange];
        % Set face and edge color and line width
        for k = 1:length(b)
            b(k).FaceColor = bar_colors(k,:);
            b(k).EdgeColor = 'black'; % Black edge color
            b(k).LineWidth = 1; % Line width of 1
        end
        if ~strcmp(variable_names{h}, 'distillation')
            xlabelName = 'Temperature [K]';
        else
            xlabelName = 'Recovered Volume Fraction [\%]';
        end
        xlabelHandle = xlabel(xlabelName);
        set(xlabelHandle, 'FontSize', fontsize, 'FontName', fontname, 'Interpreter', 'latex');
        if strcmp(variable_names{h}, 'rho')
            ylabelName = '$S_{i,\rho}$ [-]';
        end
        if strcmp(variable_names{h}, 'mu')
            ylabelName = '$S_{i,\mu}$ [-]';
        end
        if strcmp(variable_names{h}, 'nu')
            ylabelName = '$S_{i,\nu}$ [-]';
        end
        if strcmp(variable_names{h}, 'kappa')
            ylabelName = '$S_{i,k}$ [-]';
        end
        if strcmp(variable_names{h}, 'specificHeat')
            ylabelName = '$S_{i,c_{p}}$ [-]';
        end
        if strcmp(variable_names{h}, 'latentHeat')
            ylabelName = '$S_{i,L_{v}}$ [-]';
        end
        if strcmp(variable_names{h}, 'vaporPressure')
            ylabelName = '$S_{i,p_{v}}$ [-]';
        end
        if strcmp(variable_names{h}, 'sigma')
            ylabelName = '$S_{i,\sigma}$ [-]';
        end
        if strcmp(variable_names{h}, 'distillation')
            ylabelName = '$S_{i,dist}$ [-]';
        end
        ylabelHandle = ylabel(ylabelName);
        set(ylabelHandle, 'FontSize', fontsize, 'FontName', fontname, 'Interpreter', 'latex');
        ylim([0 1.05]);
        axesHandle = gca;
        set(axesHandle, 'FontSize', 20, 'FontName', fontname, 'Linewidth', line);
        legend({'$X_{i}$', '$n_{C,i}$', '$\eta_{B,i}^{*}$'}, 'FontSize', 20, 'FontName', fontname, 'Interpreter', 'latex', 'Location', 'northoutside', 'NumColumns', 3);
        set(gcf, 'Position', [100, 100, pix, pix/1.3801]);
        box on;
        folderPathSobol = './Figures/GSA/';
        % Create the folder if it doesn't exist
        if ~isfolder(folderPathSobol)
            mkdir(folderPathSobol);
        end
        filenameSobol = strcat(folderPathSobol, 'Si_', variable_names{h});
        print(filenameSobol, '-dpng', '-r600');

    else

        sobol_idx_lumped(:, counter) = sobol_lumped(families, classes, numComponents, chain_reshaped, model_output_samples, variable_names{h});
        lumped_variables{counter} = variable_names{h};

        counter = counter + 1;

    end

end

if numel(lumped_variables) > 0

    t = tiledlayout(1,1); t.Padding = 'none'; t.TileSpacing = 'none'; nexttile
    figure('visible', 'off');
    barWidth = 1;
    b = bar([sobol_idx_lumped(1,:); sobol_idx_lumped(2,:); sobol_idx_lumped(3,:)]', 'grouped', 'BarWidth', barWidth);
    % Set colors for bars
    bar_colors = [green; blue; orange];
    % Set face and edge color and line width
    for k = 1:length(b)
        b(k).FaceColor = bar_colors(k,:);
        b(k).EdgeColor = 'black'; % Black edge color
        b(k).LineWidth = 1; % Line width of 1
    end
    xticks(1:numel(lumped_variables));
    xticklabels(lumped_variables);
    ylabelName = '$S_{i,lumped}$ [-]';
    ylabelHandle = ylabel(ylabelName);
    set(ylabelHandle, 'FontSize', fontsize, 'FontName', fontname, 'Interpreter', 'latex');
    ylim([0 1.05]);
    axesHandle = gca;
    set(axesHandle, 'FontSize', 20, 'FontName', fontname, 'Linewidth', line);
    legend({'$X_{i}$', '$n_{C,i}$', '$\eta_{B,i}^{*}$'}, 'FontSize', 20, 'FontName', fontname, 'Interpreter', 'latex', 'Location', 'northoutside', 'NumColumns', 3);
    set(gcf, 'Position', [100, 100, pix, pix/1.3801]);
    box on;
    folderPathSobol = './Figures/GSA/';
    % Create the folder if it doesn't exist
    if ~isfolder(folderPathSobol)
        mkdir(folderPathSobol);
    end
    filenameSobol = strcat(folderPathSobol, 'Si_lumped');
    print(filenameSobol, '-dpng', '-r600');

end

end