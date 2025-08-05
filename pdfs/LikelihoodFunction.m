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

function L = LikelihoodFunction(families, classes, data, std, X, n, eta_B_star, variable, pressure)

% DESCRIPTION:
% The LikelihoodFunction function computes a Gaussian approximation for marginalized log-likelihoods.

% --- References ---

% Hu, J. and Burns, A., 1970.
% Index predicts cloud, pour and flash points of distillates fuel blends.
% Oil Gas J. 68(45):66.

% Riazi, R., 2005.
% Characterization and Properties of Petroleum Fractions.
% ASTM International, Philadelphia.

% Auxiliary parameters:
% Nc: number of surrogate mixture components  [-]
% Nd: number of experimental measurements for the thermophysical property under consideration  [-]

% Inputs:
% 1) families  : (1 x Nc) cell array, with the i-th cell containing a string
% that denotes the i-th hydrocarbon family under consideration
% 2) classes   : (1 x Nc) cell array, with the i-th cell containing a structure array for the
% hydrocarbon family and range of number of carbon atoms of the i-th surrogate component
% 3) data      : (Nd x 3) array containing the experimental measurements about the
% thermophysical property under consideration (independent variable,
% dependent thermophysical property, standard deviation). Note that "data"
% is a (Nd x 4) array in case the distillation curve is addressed in the
% thermophysical characterization of the surrogate mixture, with the 4th
% column denoting the operating pressure the distillation curve is computed at.
% In contrast, "data" is a (1 x 2) array in case the molecular weight or hydrogen-carbon
% ratio is addressed in the chemical characterization of the surrogate mixture, with
% the 1th and 2nd columns denoting the molecular weight or hydrogen-carbon ratio of
% the real fuel and related uncertainty, respectively.
% 4) std       : (Nd x 1) array containing the standard deviation associated with the
% experimental measurements about the thermophysical property under consideration
% 5) X         : (1 x Nc-1) array containing the molar fractions of the surrogate
% components in the liquid phase before the distillation process begins  [-]
% 6) n         : (1 x Nc) array containing the numbers of carbon atoms of the
% surrogate components  [-]
% 7) eta_B_star: (1 x Nc) array containing the topochemical atom indices
% of the surrogate components  [-]
% 8) variable  : string denoting the thermophysical property under
% consideration
% 9) pressure  : operating pressure  [Pa]

% Outputs:
% 1) L: Gaussian log-likelihood

index_n_eta = zeros(1, numel(n)); % indices of the components

for i = 1:numel(n)
    index_n_eta(i) = findIndex_eta(classes{i}, n(i), eta_B_star(i));
end

if ~strcmp(variable,'distillation') && ~strcmp(variable,'molWeight') && ~strcmp(variable,'HC') && ~strcmp(variable,'DCN') && ~strcmp(variable,'flash') && ~strcmp(variable,'freezing') && ~strcmp(variable,'LHV')

    phiModel = ThermophysicalProperties_LiquidMixture(variable, data(:,1), X, classes, index_n_eta, pressure);
    arrayLikelihood = log(2.*pi.*std.^2) + (phiModel - data(:, 2)).^2./(std.^2);
    L = - 1./(2*numel(arrayLikelihood)) * sum(arrayLikelihood);

elseif strcmp(variable,'distillation')

    [volumeFraction, T_distillation]  = DistillationCurve(X, classes, index_n_eta, pressure);
    % Check if all the temperature values along the distillation curve are real
    if all(isreal(T_distillation)) && all(isreal(volumeFraction))
        T_interpolated = interp1(volumeFraction, T_distillation, data(:, 1), 'linear');
        arrayLikelihood = log(2.*pi.*std.^2) + (T_interpolated - data(:, 2)).^2./(std.^2);
    else
        arrayLikelihood = Inf;
    end
    L = - 1./(2*numel(arrayLikelihood)) * sum(arrayLikelihood);

elseif strcmp(variable,'molWeight')

    molFrac = zeros(numel(classes),1);
    Wi = zeros(numel(classes),1);
    for i = 1:numel(classes)
        if i < numel(classes)
            molFrac(i) = X(i);
        else
            molFrac(i) = 1-sum(X);
        end
        Wi(i) = classes{i}(index_n_eta(i)).molWeight;
    end
    phiModel = molWeight(molFrac, Wi);
    arrayLikelihood = log(2.*pi.*std.^2) + (phiModel - data(:, 1)).^2./(std.^2);
    L = - 1./(2*numel(arrayLikelihood)) * sum(arrayLikelihood);

elseif strcmp(variable,'HC')

    molFrac = zeros(numel(classes),1);
    for i = 1:numel(classes)
        if i < numel(classes)
            molFrac(i) = X(i);
        else
            molFrac(i) = 1-sum(X);
        end
    end
    % Compute the molecular formula of the surrogate mixture under consideration
    nCarbon = round(sum(molFrac'.*n),2); % number of carbon atoms denoting the molecular formula
    nHydrogen = 0.0; % number of hydrogen atoms denoting the molecular formula
    for j = 1:numel(n)
        if strcmp(families(j), 'nparaffins') || strcmp(families(j), 'isoparaffins')
            nHydrogen = nHydrogen + molFrac(j)*(2*n(j)+2);
        elseif strcmp(families(j), 'cycloparaffins')
            nHydrogen = nHydrogen + molFrac(j)*2*n(j);
        elseif strcmp(families(j), 'dicycloparaffins')
            nHydrogen = nHydrogen + molFrac(j)*(2*n(j)-2);
        elseif strcmp(families(j), 'alkylbenzenes')
            nHydrogen = nHydrogen + molFrac(j)*(2*n(j)-6);
        elseif strcmp(families(j), 'alkylnaphtalenes')
            nHydrogen = nHydrogen + molFrac(j)*(2*n(j)-12);
        elseif strcmp(families(j), 'cycloaromatics')
            nHydrogen = nHydrogen + molFrac(j)*(2*n(j)-8);
        end
    end
    nHydrogen = round(nHydrogen, 2);
    phiModel = round(nHydrogen/nCarbon, 3);
    arrayLikelihood = log(2.*pi.*std.^2) + (phiModel - data(:, 1)).^2./(std.^2);
    L = - 1./(2*numel(arrayLikelihood)) * sum(arrayLikelihood);

elseif strcmp(variable,'DCN')

    molFrac = zeros(numel(classes),1);
    Wi = zeros(numel(classes),1);
    DCN_i = zeros(numel(classes),1);
    for i = 1:numel(classes)
        if i < numel(classes)
            molFrac(i) = X(i);
        else
            molFrac(i) = 1-sum(X);
        end
        Wi(i) = classes{i}(index_n_eta(i)).molWeight;
        DCN_i(i) = classes{i}(index_n_eta(i)).DCN;
    end

    rho_i = zeros(1,numel(classes));
    numerator = zeros(1,numel(classes));

    for i = 1:numel(classes)
        rho_i(i) = ThermophysicalProperties_SingleLiquid('rho', 300, classes{i}(index_n_eta(i)).molWeight, classes{i}(index_n_eta(i)).Tc, classes{i}(index_n_eta(i)).coeffRho(1), classes{i}(index_n_eta(i)).coeffRho(2), classes{i}(index_n_eta(i)).coeffRho(3), classes{i}(index_n_eta(i)).coeffMu(1), classes{i}(index_n_eta(i)).coeffMu(2), classes{i}(index_n_eta(i)).coeffMu(3), classes{i}(index_n_eta(i)).coeffMu(4), classes{i}(index_n_eta(i)).coeffK(1), classes{i}(index_n_eta(i)).coeffK(2), classes{i}(index_n_eta(i)).coeffK(3), classes{i}(index_n_eta(i)).coeffCl(1), classes{i}(index_n_eta(i)).coeffCl(2), classes{i}(index_n_eta(i)).coeffCl(3), classes{i}(index_n_eta(i)).coeffCl(4), classes{i}(index_n_eta(i)).coeffHv(1), classes{i}(index_n_eta(i)).coeffHv(2), classes{i}(index_n_eta(i)).coeffPsat(1), classes{i}(index_n_eta(i)).coeffPsat(2), classes{i}(index_n_eta(i)).coeffPsat(3), classes{i}(index_n_eta(i)).coeffPsat(4), classes{i}(index_n_eta(i)).coeffPsat(5), classes{i}(index_n_eta(i)).coeffSigma(1), classes{i}(index_n_eta(i)).coeffSigma(2), pressure);
        numerator(i) = molFrac(i)*classes{i}(index_n_eta(i)).molWeight;
    end

    % Compute the mixture density
    rho_mixture = sum(numerator, 2)./sum(numerator./ rho_i, 2);

    % Convert mass fractions to volume fractions
    Yi = mol2mass(molFrac, Wi);
    Vi = rho_mixture*Yi./rho_i';

    phiModel = sum(Vi.*DCN_i);
    arrayLikelihood = log(2.*pi.*std.^2) + (phiModel - data(:, 1)).^2./(std.^2);
    L = - 1./(2*numel(arrayLikelihood)) * sum(arrayLikelihood);

elseif strcmp(variable,'flash')

    % --- Blending rule from Gary and Handwerk (2001) to estimate the mixture flash point --- %
    molFrac = zeros(numel(classes),1);
    Wi = zeros(numel(classes),1);
    flash_i = zeros(numel(classes),1);
    BI_flash_i = zeros(numel(classes),1);
    for i = 1:numel(classes)
        if i < numel(classes)
            molFrac(i) = X(i);
        else
            molFrac(i) = 1-sum(X);
        end
        Wi(i) = classes{i}(index_n_eta(i)).molWeight;
        flash_i(i) = classes{i}(index_n_eta(i)).Tf;
        flash_F = (flash_i(i)-273.15)*9./5+32;
        BI_flash_i(i) = 51708*exp((log(flash_F)-2.6287)^2/(-0.91725));
    end

    % Convert molar fractions to mass fractions
    Yi = mol2mass(molFrac, Wi);

    BI_blend = sum(Yi.*BI_flash_i);

    phiModel_F = exp(((-0.91725)*log(BI_blend/51708))^0.5+2.6287); % mixture flash point in [°F]
    phiModel = (phiModel_F-32)*5./9+273.15; % mixture flash point in [K]

    arrayLikelihood = log(2.*pi.*std.^2) + (phiModel - data(:, 1)).^2./(std.^2);
    L = - 1./(2*numel(arrayLikelihood)) * sum(arrayLikelihood);

elseif strcmp(variable,'freezing')

    % --- Blending rule from Hu and Burns (1970), Riazi (2005) to estimate the mixture freezing point --- %
    molFrac = zeros(numel(classes),1);
    Wi = zeros(numel(classes),1);
    freezing_i = zeros(numel(classes),1);
    BI_freezing_i = zeros(numel(classes),1);
    rho_i = zeros(1,numel(classes));
    numerator = zeros(1,numel(classes));
    for i = 1:numel(classes)
        if i < numel(classes)
            molFrac(i) = X(i);
        else
            molFrac(i) = 1-sum(X);
        end
        Wi(i) = classes{i}(index_n_eta(i)).molWeight;
        rho_i(i) = ThermophysicalProperties_SingleLiquid('rho', 300, classes{i}(index_n_eta(i)).molWeight, classes{i}(index_n_eta(i)).Tc, classes{i}(index_n_eta(i)).coeffRho(1), classes{i}(index_n_eta(i)).coeffRho(2), classes{i}(index_n_eta(i)).coeffRho(3), classes{i}(index_n_eta(i)).coeffMu(1), classes{i}(index_n_eta(i)).coeffMu(2), classes{i}(index_n_eta(i)).coeffMu(3), classes{i}(index_n_eta(i)).coeffMu(4), classes{i}(index_n_eta(i)).coeffK(1), classes{i}(index_n_eta(i)).coeffK(2), classes{i}(index_n_eta(i)).coeffK(3), classes{i}(index_n_eta(i)).coeffCl(1), classes{i}(index_n_eta(i)).coeffCl(2), classes{i}(index_n_eta(i)).coeffCl(3), classes{i}(index_n_eta(i)).coeffCl(4), classes{i}(index_n_eta(i)).coeffHv(1), classes{i}(index_n_eta(i)).coeffHv(2), classes{i}(index_n_eta(i)).coeffPsat(1), classes{i}(index_n_eta(i)).coeffPsat(2), classes{i}(index_n_eta(i)).coeffPsat(3), classes{i}(index_n_eta(i)).coeffPsat(4), classes{i}(index_n_eta(i)).coeffPsat(5), classes{i}(index_n_eta(i)).coeffSigma(1), classes{i}(index_n_eta(i)).coeffSigma(2), 101325);
        numerator(i) = molFrac(i)*classes{i}(index_n_eta(i)).molWeight;
        freezing_i(i) = classes{i}(index_n_eta(i)).Tfz;
        BI_freezing_i(i) = freezing_i(i)^(1/0.05);
    end

    % Convert molar fractions to mass fractions
    Yi = mol2mass(molFrac, Wi);
    rho_mixture = sum(numerator, 2)./sum(numerator./ rho_i, 2);
    Vi = rho_mixture*Yi./rho_i'; % volume fractions

    BI_blend = max(0, sum(Vi.*BI_freezing_i));

    phiModel = BI_blend^0.05; % mixture freezing point [K]

    arrayLikelihood = log(2.*pi.*std.^2) + (phiModel - data(:, 1)).^2./(std.^2);
    L = - 1./(2*numel(arrayLikelihood)) * sum(arrayLikelihood);

elseif strcmp(variable,'LHV')

    molFrac = zeros(numel(classes),1);
    Wi = zeros(numel(classes),1);
    LHV_i= zeros(1,numel(classes));

    for i = 1:numel(classes)

        if i < numel(classes)
            molFrac(i) = X(i);
        else
            molFrac(i) = 1-sum(X);
        end
        Wi(i) = classes{i}(index_n_eta(i)).molWeight;

        LHV_i(i) = classes{i}(index_n_eta(i)).Hc; % LHV of the i-th surrogate component [MJ/kg]

    end

    % Convert molar fractions to mass fractions
    Yi = mol2mass(molFrac, Wi);

    phiModel = sum(Yi'.*LHV_i);

    arrayLikelihood = log(2.*pi.*std.^2) + (phiModel - data(:, 1)).^2./(std.^2);
    L = - 1./(2*numel(arrayLikelihood)) * sum(arrayLikelihood);

end

end
