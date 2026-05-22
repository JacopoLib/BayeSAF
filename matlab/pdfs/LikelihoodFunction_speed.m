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
% LikelihoodFunction_speed evaluates the vectorized Gaussian log-likelihood
% for multiple MCMC samples simultaneously. Supported observable types:
% molWeight, HC (H/C ratio), DCN, LHV, flash point, freezing point,
% distillation curve, deltaT_dist, and generic thermophysical properties
% (routed to ThermophysicalProperties_LiquidMixture_speed).
%
% Inputs:
% 1) families        : (1 x Nc) cell array of family strings
% 2) classes         : (1 x Nc) cell array of struct arrays
% 3) data            : experimental data matrix
% 4) stdData         : standard deviations of experimental data
% 5) X_all           : (numSamples x Nc-1) molar fractions  [-]
% 6) n_all           : (numSamples x Nc) carbon-atom numbers  [-]
% 7) eta_B_star_all  : (numSamples x Nc) normalised topochemical indices  [-]
% 8) variable        : string identifier of the observable
% 9) pressure        : pressure  [Pa]
%
% Outputs:
% 1) L: (numSamples x 1) log-likelihood values  [-]
% ------------------------------------------------------------------------

function L = LikelihoodFunction_speed(families, classes, data, stdData, X_all, n_all, eta_B_star_all, variable, pressure)

[numSamples, Nc_minus1] = size(X_all);
Nc = Nc_minus1 + 1;

molFracMat = zeros(numSamples, Nc);
molFracMat(:, 1:Nc-1) = X_all;
molFracMat(:, Nc)     = 1 - sum(X_all, 2);

index_n_eta = zeros(numSamples, Nc);
for j = 1:Nc
    for k = 1:numSamples
        index_n_eta(k,j) = findIndex_eta(classes{j}, n_all(k,j), eta_B_star_all(k,j));
    end
end

persistent cache

switch variable
    case 'molWeight'
        WiMat = zeros(numSamples, Nc);
        for j = 1:Nc
            for k = 1:numSamples
                WiMat(k,j) = classes{j}(index_n_eta(k,j)).molWeight;
            end
        end
        phiModel = sum(molFracMat .* WiMat, 2);

        stdData  = reshape(stdData, 1, []);
        data_mat = repmat(data(:,1)', numSamples, 1);
        std_mat  = repmat(stdData,    numSamples, 1);
        arrayLikelihood = log(2*pi*std_mat.^2) + ((phiModel - data_mat).^2) ./ (std_mat.^2);
        L = -0.5 * sum(arrayLikelihood, 2);

    case 'HC'
        nCarbon   = round(sum(molFracMat .* n_all, 2), 2);
        nHydrogen = zeros(numSamples, 1);
        for j = 1:Nc
            switch families{j}
                case {'nparaffins','isoparaffins'}
                    nHydrogen = nHydrogen + molFracMat(:,j) .* (2*n_all(:,j)+2);
                case 'cycloparaffins'
                    nHydrogen = nHydrogen + molFracMat(:,j) .* (2*n_all(:,j));
                case 'dicycloparaffins'
                    nHydrogen = nHydrogen + molFracMat(:,j) .* (2*n_all(:,j)-2);
                case 'alkylbenzenes'
                    nHydrogen = nHydrogen + molFracMat(:,j) .* (2*n_all(:,j)-6);
                case 'alkylnaphtalenes'
                    nHydrogen = nHydrogen + molFracMat(:,j) .* (2*n_all(:,j)-12);
                case 'cycloaromatics'
                    nHydrogen = nHydrogen + molFracMat(:,j) .* (2*n_all(:,j)-8);
            end
        end
        phiModel = round(nHydrogen ./ nCarbon, 3);

        stdData  = reshape(stdData, 1, []);
        data_mat = repmat(data(:,1)', numSamples, 1);
        std_mat  = repmat(stdData,    numSamples, 1);
        arrayLikelihood = log(2*pi*std_mat.^2) + ((phiModel - data_mat).^2) ./ (std_mat.^2);
        L = -0.5 * sum(arrayLikelihood, 2);

    case 'DCN'
        WiMat  = zeros(numSamples, Nc);
        DCNMat = zeros(numSamples, Nc);
        rhoMat = zeros(numSamples, Nc);
        for j = 1:Nc
            for k = 1:numSamples
                c = classes{j}(index_n_eta(k,j));
                WiMat(k,j)  = c.molWeight;
                DCNMat(k,j) = c.DCN;
                rhoMat(k,j) = ThermophysicalProperties_SingleLiquid('rho', 300, c.molWeight, c.Tc, ...
                    c.coeffRho(1), c.coeffRho(2), c.coeffRho(3), ...
                    c.coeffMu(1),  c.coeffMu(2),  c.coeffMu(3),  c.coeffMu(4), ...
                    c.coeffK(1),   c.coeffK(2),   c.coeffK(3), ...
                    c.coeffCl(1),  c.coeffCl(2),  c.coeffCl(3),  c.coeffCl(4), ...
                    c.coeffHv(1),  c.coeffHv(2), ...
                    c.coeffPsat(1), c.coeffPsat(2), c.coeffPsat(3), c.coeffPsat(4), c.coeffPsat(5), ...
                    c.coeffSigma(1), c.coeffSigma(2), pressure);
            end
        end
        rho_mixture = sum(molFracMat .* WiMat, 2) ./ sum((molFracMat .* WiMat) ./ rhoMat, 2);
        Yi = zeros(numSamples, Nc);
        for jj = 1:numSamples
            Yi(jj,:) = mol2mass(molFracMat(jj,:), WiMat(jj,:));
        end
        Vi       = rho_mixture .* Yi ./ rhoMat;
        phiModel = sum(Vi .* DCNMat, 2);

        stdData  = reshape(stdData, 1, []);
        data_mat = repmat(data(:,1)', numSamples, 1);
        std_mat  = repmat(stdData,    numSamples, 1);
        arrayLikelihood = log(2*pi*std_mat.^2) + ((phiModel - data_mat).^2) ./ (std_mat.^2);
        L = -0.5 * sum(arrayLikelihood, 2);

    case 'LHV'
        WiMat  = zeros(numSamples, Nc);
        LHVMat = zeros(numSamples, Nc);
        for j = 1:Nc
            for k = 1:numSamples
                c = classes{j}(index_n_eta(k,j));
                WiMat(k,j)  = c.molWeight;
                LHVMat(k,j) = c.Hc;
            end
        end
        Yi = zeros(numSamples, Nc);
        for jj = 1:numSamples
            Yi(jj,:) = mol2mass(molFracMat(jj,:), WiMat(jj,:));
        end
        phiModel = sum(Yi .* LHVMat, 2);

        stdData  = reshape(stdData, 1, []);
        data_mat = repmat(data(:,1)', numSamples, 1);
        std_mat  = repmat(stdData,    numSamples, 1);
        arrayLikelihood = log(2*pi*std_mat.^2) + ((phiModel - data_mat).^2) ./ (std_mat.^2);
        L = -0.5 * sum(arrayLikelihood, 2);

    case 'flash'
        WiMat    = zeros(numSamples, Nc);
        BI_flash = zeros(numSamples, Nc);
        for j = 1:Nc
            for k = 1:numSamples
                c = classes{j}(index_n_eta(k,j));
                WiMat(k,j)    = c.molWeight;
                flashF        = (c.Tf - 273.15) * 9/5 + 32;
                BI_flash(k,j) = 51708 * exp((log(flashF) - 2.6287)^2 / (-0.91725));
            end
        end
        Yi = zeros(numSamples, Nc);
        for jj = 1:numSamples
            Yi(jj,:) = mol2mass(molFracMat(jj,:), WiMat(jj,:));
        end
        BI_blend    = sum(Yi .* BI_flash, 2);
        phiModel_F  = exp(((-0.91725) * log(BI_blend/51708)).^0.5 + 2.6287);
        phiModel    = (phiModel_F - 32) * 5/9 + 273.15;

        stdData  = reshape(stdData, 1, []);
        data_mat = repmat(data(:,1)', numSamples, 1);
        std_mat  = repmat(stdData,    numSamples, 1);
        arrayLikelihood = log(2*pi*std_mat.^2) + ((phiModel - data_mat).^2) ./ (std_mat.^2);
        L = -0.5 * sum(arrayLikelihood, 2);

    case 'freezing'
        WiMat       = zeros(numSamples, Nc);
        rhoMat      = zeros(numSamples, Nc);
        BI_freezing = zeros(numSamples, Nc);
        for j = 1:Nc
            for k = 1:numSamples
                c = classes{j}(index_n_eta(k,j));
                WiMat(k,j)       = c.molWeight;
                rhoMat(k,j)      = ThermophysicalProperties_SingleLiquid('rho', 300, c.molWeight, c.Tc, ...
                    c.coeffRho(1), c.coeffRho(2), c.coeffRho(3), ...
                    c.coeffMu(1),  c.coeffMu(2),  c.coeffMu(3),  c.coeffMu(4), ...
                    c.coeffK(1),   c.coeffK(2),   c.coeffK(3), ...
                    c.coeffCl(1),  c.coeffCl(2),  c.coeffCl(3),  c.coeffCl(4), ...
                    c.coeffHv(1),  c.coeffHv(2), ...
                    c.coeffPsat(1), c.coeffPsat(2), c.coeffPsat(3), c.coeffPsat(4), c.coeffPsat(5), ...
                    c.coeffSigma(1), c.coeffSigma(2), pressure);
                BI_freezing(k,j) = c.Tfz^(1/0.05);
            end
        end
        Yi = zeros(numSamples, Nc);
        for jj = 1:numSamples
            Yi(jj,:) = mol2mass(molFracMat(jj,:), WiMat(jj,:));
        end
        rho_mixture = sum(molFracMat .* WiMat, 2) ./ sum((molFracMat .* WiMat) ./ rhoMat, 2);
        Vi          = rho_mixture .* Yi ./ rhoMat;
        BI_blend    = max(0, sum(Vi .* BI_freezing, 2));
        phiModel    = BI_blend.^0.05;

        stdData  = reshape(stdData, 1, []);
        data_mat = repmat(data(:,1)', numSamples, 1);
        std_mat  = repmat(stdData,    numSamples, 1);
        arrayLikelihood = log(2*pi*std_mat.^2) + ((phiModel - data_mat).^2) ./ (std_mat.^2);
        L = -0.5 * sum(arrayLikelihood, 2);

    case 'distillation'
        [volFrac_all, T_dist_all] = DistillationCurve_speed(X_all, classes, index_n_eta, pressure);

        cache.volFrac_all = volFrac_all;
        cache.T_dist_all  = T_dist_all;
        cache.pressure    = pressure;
        cache.index_n_eta = index_n_eta;
        cache.X_all       = X_all;

        T_target  = data(:, 1);
        N_exp     = numel(T_target);
        N_samples = size(volFrac_all, 1);

        T_interpolated = zeros(N_samples, N_exp);
        parfor s = 1:N_samples
            vf_s = real(volFrac_all(s,:));
            T_s  = real(T_dist_all(s,:));
            valid = isfinite(vf_s) & isfinite(T_s);
            if sum(valid) >= 2
                T_interpolated(s,:) = interp1(vf_s(valid), T_s(valid), T_target, 'linear', 'extrap');
            else
                T_interpolated(s,:) = NaN;
            end
        end

        data_mat = repmat(data(:,2)',  N_samples, 1);
        std_mat  = repmat(reshape(stdData, 1, []), N_samples, 1);
        arrayLikelihood = log(2*pi*std_mat.^2) + ((T_interpolated - data_mat).^2) ./ (std_mat.^2);
        L = -0.5 * sum(arrayLikelihood, 2);
        L(any(~isfinite(T_interpolated), 2)) = -Inf;

    case 'deltaT_dist'
        if isempty(cache) || ~isequal(cache.pressure, pressure)
            [volFrac_all, T_dist_all] = DistillationCurve_speed(X_all, classes, index_n_eta, pressure);
        else
            volFrac_all = cache.volFrac_all;
            T_dist_all  = cache.T_dist_all;
        end

        N_samples = size(volFrac_all, 1);
        T10 = zeros(N_samples, 1);
        T50 = zeros(N_samples, 1);
        T90 = zeros(N_samples, 1);
        parfor s = 1:N_samples
            vf_s = real(volFrac_all(s,:));
            T_s  = real(T_dist_all(s,:));
            valid = isfinite(vf_s) & isfinite(T_s);
            if sum(valid) >= 2
                T10(s) = interp1(vf_s(valid), T_s(valid), 10, 'linear', 'extrap');
                T50(s) = interp1(vf_s(valid), T_s(valid), 50, 'linear', 'extrap');
                T90(s) = interp1(vf_s(valid), T_s(valid), 90, 'linear', 'extrap');
            else
                T10(s) = NaN;
                T50(s) = NaN;
                T90(s) = NaN;
            end
        end

        T5010 = T50 - T10;
        T9010 = T90 - T10;

        phiModel = zeros(N_samples, numel(data(:,1)));
        for i = 1:numel(data(:,1))
            if data(i,1) == 5010
                phiModel(:,i) = T5010;
            elseif data(i,1) == 9010
                phiModel(:,i) = T9010;
            end
        end

        stdData  = reshape(stdData, 1, []);
        data_mat = repmat(data(:,2)', N_samples, 1);
        std_mat  = repmat(stdData,    N_samples, 1);
        arrayLikelihood = log(2*pi*std_mat.^2) + ((phiModel - data_mat).^2) ./ (std_mat.^2);
        L = -0.5 * sum(arrayLikelihood, 2);

    otherwise
        phiModel = ThermophysicalProperties_LiquidMixture_speed(variable, data(:,1), X_all, classes, index_n_eta, pressure);

        stdData  = reshape(stdData, 1, []);
        data_mat = repmat(data(:,2)', size(phiModel,1), 1);
        std_mat  = repmat(stdData,    size(phiModel,1), 1);
        arrayLikelihood = log(2*pi*std_mat.^2) + ((phiModel - data_mat).^2) ./ (std_mat.^2);
        L = -0.5 * sum(arrayLikelihood, 2);
end

end
