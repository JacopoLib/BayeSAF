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

function sobol_idx = sobol_lumped(families, classes, numComponents, chain_reshaped, model_output_samples, variable_name)

% DESCRIPTION:
% The sobol_lumped function calculates grouped first-order Sobol' indices based on the posterior PDF about three sets of
% input parameters and their impact on lumped mixture properties (molecular weight, hydrogen-to-carbon ratio, and derived
% cetane number): i) molar fractions, ii) numbers of carbon atoms, and iii) topochemical atom indices.
% Note that first-order Sobol' indices are calculated according to the algorithm proposed in Sobol' and Myshetskaya (2008).

% --- References ---

% Sobol', I.M. and Myshetskaya, E.E., 2008.
% Monte Carlo estimators for small sensitivity indices.
% Monte Carlo Methods Appl 13 (5-6), 455–65.

% Hu, J. and Burns, A., 1970.
% Index predicts cloud, pour and flash points of distillates fuel blends.
% Oil Gas J. 68(45):66.

% Riazi, R., 2005.
% Characterization and Properties of Petroleum Fractions.
% ASTM International, Philadelphia.

% Auxiliary parameters:
% N_chains     : number of DE-MC chains
% t_convergence: number of iterations denoting the number of iterations the convergence of the DE-MC run
% is reached at according to the R-hat statistic ( = maxIterations in case convergence is not reached)
% t_burnin     : number of iterations to be discarded as burn-in
% maxIterations: maximum number of iterations for each chain in the DE-MC algorithm

% Inputs:
% 1) families            : (1 x numComponents) cell array, with the i-th cell containing a string
% that denotes the i-th hydrocarbon family under consideration
% 2) classes             : (1 x numComponents) cell array, with the i-th cell containing a structure array
% for the hydrocarbon family and range of number of carbon atoms of the i-th surrogate component
% 3) numComponents       : number of surrogate mixture components  [-]
% 4) chain_reshaped      : ((t_convergence - t_burnin)*N_chains x 3*numComponents-1) array containing the
% samples from the posterior PDF from each chain discarding the burn-in period, concerning
% the molar fractions, numbers of carbon atoms, and topochemical atom indices
% 5) model_output_samples: (t_convergence - t_burnin)*N_chains x 1) array describing the evolution of the
% lumped property under consideration for each sample in chain_reshaped
% 6) variable_name       : string denoting the lumped property under consideration

% Outputs:
% 1) sobol_idx: (3 x 1) array containing the grouped first-order Sobol' indices

rng('default'); % to ensure reproducibility

% Remove duplicate samples
[samples,samples_indices]=unique(round(chain_reshaped,6),'rows');
% Remove model evaluations corresponding to duplicate samples
new_output = model_output_samples(samples_indices);

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
y_A = new_output(indicesA);
y_B = new_output(indicesB);

% Variance and mean of the overall model evaluations
Vy = var([y_A; y_B]);
c = mean([y_A; y_B]);
% Initialize vector containing first-order Sobol' indices
sobol_idx = zeros(3, 1);

% +++ First-order Sobol' indices for molar fractions group +++ %
A_B = B;
y_AB = zeros(size(A_B,1), 1);
A_B(:, 1:numComponents-1) = A(:, 1:numComponents-1);
molFrac_AB = zeros(size(A_B,1), numComponents);
nC_AB = zeros(size(A_B,1), numComponents);
for p = 1:size(A_B,1)
    for q = 1:numComponents
        index_n_eta_AB(q) = findIndex_eta(classes{q}, A_B(p, numComponents+q-1), A_B(p, 2*numComponents+q-1));
    end

    molFrac_AB(p,:) = [A_B(p,1:numComponents-1) 1-sum(A_B(p,1:numComponents-1))];

    if strcmp(variable_name, 'molWeight')
        Wl = zeros(1, numComponents);
        for l = 1:numComponents
            Wl(l) = classes{l}(index_n_eta_AB(l)).molWeight;
        end
        y_AB(p) = molWeight(molFrac_AB(p,:), Wl);

    elseif strcmp(variable_name, 'HC')
        nC_AB(p,:) = A_B(p,numComponents:2*numComponents-1);
        nCarbon = round(sum(molFrac_AB(p,:).*nC_AB(p,:)),2);
        nHydrogen = 0.0;
        for l = 1:numComponents
            if strcmp(families(l), 'nparaffins') || strcmp(families(l), 'isoparaffins')
                nHydrogen = nHydrogen + molFrac_AB(p,l)*(2*nC_AB(p,l)+2);
            elseif strcmp(families(l), 'cycloparaffins')
                nHydrogen = nHydrogen + molFrac_AB(p,l)*2*nC_AB(p,l);
            elseif strcmp(families(l), 'dicycloparaffins')
                nHydrogen = nHydrogen + molFrac_AB(p,l)*(2*nC_AB(p,l)-2);
            elseif strcmp(families(l), 'alkylbenzenes')
                nHydrogen = nHydrogen + molFrac_AB(p,l)*(2*nC_AB(p,l)-6);
            elseif strcmp(families(l), 'alkylnaphtalenes')
                nHydrogen = nHydrogen + molFrac_AB(p,l)*(2*nC_AB(p,l)-12);
            elseif strcmp(families(l), 'cycloaromatics')
                nHydrogen = nHydrogen + molFrac_AB(p,l)*(2*nC_AB(p,l)-8);
            end
        end
        nHydrogen = round(nHydrogen, 2);
        y_AB(p) = round(nHydrogen/nCarbon, 3);

    elseif strcmp(variable_name, 'DCN')
        rho_l = zeros(1,numComponents);
        numerator = zeros(1,numComponents);
        DCN_l = zeros(1,numComponents);

        Wl = zeros(1, numComponents);

        for l = 1:numComponents

            rho_l(l) = ThermophysicalProperties_SingleLiquid('rho', 300, classes{l}(index_n_eta_AB(l)).molWeight, classes{l}(index_n_eta_AB(l)).Tc, classes{l}(index_n_eta_AB(l)).coeffRho(1), classes{l}(index_n_eta_AB(l)).coeffRho(2), classes{l}(index_n_eta_AB(l)).coeffRho(3), classes{l}(index_n_eta_AB(l)).coeffMu(1), classes{l}(index_n_eta_AB(l)).coeffMu(2), classes{l}(index_n_eta_AB(l)).coeffMu(3), classes{l}(index_n_eta_AB(l)).coeffMu(4), classes{l}(index_n_eta_AB(l)).coeffK(1), classes{l}(index_n_eta_AB(l)).coeffK(2), classes{l}(index_n_eta_AB(l)).coeffK(3), classes{l}(index_n_eta_AB(l)).coeffCl(1), classes{l}(index_n_eta_AB(l)).coeffCl(2), classes{l}(index_n_eta_AB(l)).coeffCl(3), classes{l}(index_n_eta_AB(l)).coeffCl(4), classes{l}(index_n_eta_AB(l)).coeffHv(1), classes{l}(index_n_eta_AB(l)).coeffHv(2), classes{l}(index_n_eta_AB(l)).coeffPsat(1), classes{l}(index_n_eta_AB(l)).coeffPsat(2), classes{l}(index_n_eta_AB(l)).coeffPsat(3), classes{l}(index_n_eta_AB(l)).coeffPsat(4), classes{l}(index_n_eta_AB(l)).coeffPsat(5), classes{l}(index_n_eta_AB(l)).coeffSigma(1), classes{l}(index_n_eta_AB(l)).coeffSigma(2), 101325);
            Wl(l) = classes{l}(index_n_eta_AB(l)).molWeight;
            numerator(l) = molFrac_AB(p,l)*classes{l}(index_n_eta_AB(l)).molWeight;
            DCN_l(l) = classes{l}(index_n_eta_AB(l)).DCN;

        end

        y_l = mol2mass(molFrac_AB(p,:), Wl);

        rho_mixture = sum(numerator, 2)./sum(numerator./ rho_l, 2);
        Vl = rho_mixture*y_l'./rho_l';
        y_AB(p) = sum(Vl.*DCN_l');

    elseif strcmp(variable_name, 'flash')
        flash_l = zeros(1,numComponents);
        BI_flash_l = zeros(1,numComponents); % blending index for flash point from Gary and Handwerk (2001)
        Wl = zeros(1, numComponents);

        for l = 1:numComponents

            Wl(l) = classes{l}(index_n_eta_AB(l)).molWeight;
            flash_l(l) = classes{l}(index_n_eta_AB(l)).Tf;
            flash_F = (flash_l(l)-273.15)*9./5+32;
            BI_flash_l(l) = 51708*exp((log(flash_F)-2.6287)^2/(-0.91725));

        end

        y_l = mol2mass(molFrac_AB(p,:), Wl);

        BI_blend = sum(y_l.*BI_flash_l);

        y_AB_F = exp(((-0.91725)*log(BI_blend/51708))^0.5+2.6287); % mixture flash point [°F]
        y_AB(p) = (y_AB_F-32)*5./9+273.15; % mixture flash point [K]

    elseif strcmp(variable_name, 'freezing')
        freezing_l = zeros(1,numComponents);
        BI_freezing_l = zeros(1,numComponents); % blending index for freezing point from Gary and Handwerk (2001)
        Wl = zeros(1, numComponents);
        rho_l = zeros(1,numComponents);
        numerator = zeros(1,numComponents);

        for l = 1:numComponents

            rho_l(l) = ThermophysicalProperties_SingleLiquid('rho', 300, classes{l}(index_n_eta_AB(l)).molWeight, classes{l}(index_n_eta_AB(l)).Tc, classes{l}(index_n_eta_AB(l)).coeffRho(1), classes{l}(index_n_eta_AB(l)).coeffRho(2), classes{l}(index_n_eta_AB(l)).coeffRho(3), classes{l}(index_n_eta_AB(l)).coeffMu(1), classes{l}(index_n_eta_AB(l)).coeffMu(2), classes{l}(index_n_eta_AB(l)).coeffMu(3), classes{l}(index_n_eta_AB(l)).coeffMu(4), classes{l}(index_n_eta_AB(l)).coeffK(1), classes{l}(index_n_eta_AB(l)).coeffK(2), classes{l}(index_n_eta_AB(l)).coeffK(3), classes{l}(index_n_eta_AB(l)).coeffCl(1), classes{l}(index_n_eta_AB(l)).coeffCl(2), classes{l}(index_n_eta_AB(l)).coeffCl(3), classes{l}(index_n_eta_AB(l)).coeffCl(4), classes{l}(index_n_eta_AB(l)).coeffHv(1), classes{l}(index_n_eta_AB(l)).coeffHv(2), classes{l}(index_n_eta_AB(l)).coeffPsat(1), classes{l}(index_n_eta_AB(l)).coeffPsat(2), classes{l}(index_n_eta_AB(l)).coeffPsat(3), classes{l}(index_n_eta_AB(l)).coeffPsat(4), classes{l}(index_n_eta_AB(l)).coeffPsat(5), classes{l}(index_n_eta_AB(l)).coeffSigma(1), classes{l}(index_n_eta_AB(l)).coeffSigma(2), 101325);
            Wl(l) = classes{l}(index_n_eta_AB(l)).molWeight;
            numerator(l) = molFrac_AB(p,l)*classes{l}(index_n_eta_AB(l)).molWeight;
            freezing_l(l) = classes{l}(index_n_eta_AB(l)).Tfz; % freezing point of the l-th surrogate component [K]
            BI_freezing_l(l) = freezing_l(l)^(1/0.05); % blending index as from Hu and Burns (1970), Riazi (2005)

        end

        y_l = mol2mass(molFrac_AB(p,:), Wl); % mass fractions
        rho_mixture = sum(numerator, 2)./sum(numerator./ rho_l, 2);
        Vl = rho_mixture*y_l./rho_l; % volume fractions

        BI_blend = sum(Vl.*BI_freezing_l);

        y_AB(p) = BI_blend^0.05; % mixture freezing point [K]

    elseif strcmp(variable_name, 'LHV')
        LHV_l = zeros(1,numComponents);
        Wl = zeros(1, numComponents);

        for l = 1:numComponents

            Wl(l) = classes{l}(index_n_eta_AB(l)).molWeight;
            LHV_l(l) = classes{l}(index_n_eta_AB(l)).Hc; % LHV of the l-th surrogate component [MJ/kg]

        end

        y_l = mol2mass(molFrac_AB(p,:), Wl); % mass fractions

        y_AB(p) = sum(y_l.*LHV_l);

    end
end

V_1 = abs(mean((y_A - c).*(y_AB - y_B)));
sobol_idx(1) = V_1 / Vy;

% +++ First-order Sobol' indices for number of carbon atoms group +++ %
A_B = B;
y_AB = zeros(size(A_B,1), 1);
A_B(:, numComponents:2*numComponents-1) = A(:, numComponents:2*numComponents-1);
for p = 1:size(A_B,1)
    for q = 1:numComponents
        index_n_eta_AB(q) = findIndex_eta(classes{q}, A_B(p, numComponents+q-1), A_B(p, 2*numComponents+q-1));
    end

    molFrac_AB(p,:) = [A_B(p,1:numComponents-1) 1-sum(A_B(p,1:numComponents-1))];

    if strcmp(variable_name, 'molWeight')
        Wl = zeros(1, numComponents);
        for l = 1:numComponents
            Wl(l) = classes{l}(index_n_eta_AB(l)).molWeight;
        end
        y_AB(p) = molWeight(molFrac_AB(p,:), Wl);

    elseif strcmp(variable_name, 'HC')
        nC_AB(p,:) = A_B(p,numComponents:2*numComponents-1);
        nCarbon = round(sum(molFrac_AB(p,:).*nC_AB(p,:)),2);
        nHydrogen = 0.0;
        for l = 1:numComponents
            if strcmp(families(l), 'nparaffins') || strcmp(families(l), 'isoparaffins')
                nHydrogen = nHydrogen + molFrac_AB(p,l)*(2*nC_AB(p,l)+2);
            elseif strcmp(families(l), 'cycloparaffins')
                nHydrogen = nHydrogen + molFrac_AB(p,l)*2*nC_AB(p,l);
            elseif strcmp(families(l), 'dicycloparaffins')
                nHydrogen = nHydrogen + molFrac_AB(p,l)*(2*nC_AB(p,l)-2);
            elseif strcmp(families(l), 'alkylbenzenes')
                nHydrogen = nHydrogen + molFrac_AB(p,l)*(2*nC_AB(p,l)-6);
            elseif strcmp(families(l), 'alkylnaphtalenes')
                nHydrogen = nHydrogen + molFrac_AB(p,l)*(2*nC_AB(p,l)-12);
            elseif strcmp(families(l), 'cycloaromatics')
                nHydrogen = nHydrogen + molFrac_AB(p,l)*(2*nC_AB(p,l)-8);
            end
        end
        nHydrogen = round(nHydrogen, 2);
        y_AB(p) = round(nHydrogen/nCarbon, 3);

    elseif strcmp(variable_name, 'DCN')
        rho_l = zeros(1,numComponents);
        numerator = zeros(1,numComponents);
        DCN_l = zeros(1,numComponents);

        Wl = zeros(1, numComponents);

        for l = 1:numComponents

            rho_l(l) = ThermophysicalProperties_SingleLiquid('rho', 300, classes{l}(index_n_eta_AB(l)).molWeight, classes{l}(index_n_eta_AB(l)).Tc, classes{l}(index_n_eta_AB(l)).coeffRho(1), classes{l}(index_n_eta_AB(l)).coeffRho(2), classes{l}(index_n_eta_AB(l)).coeffRho(3), classes{l}(index_n_eta_AB(l)).coeffMu(1), classes{l}(index_n_eta_AB(l)).coeffMu(2), classes{l}(index_n_eta_AB(l)).coeffMu(3), classes{l}(index_n_eta_AB(l)).coeffMu(4), classes{l}(index_n_eta_AB(l)).coeffK(1), classes{l}(index_n_eta_AB(l)).coeffK(2), classes{l}(index_n_eta_AB(l)).coeffK(3), classes{l}(index_n_eta_AB(l)).coeffCl(1), classes{l}(index_n_eta_AB(l)).coeffCl(2), classes{l}(index_n_eta_AB(l)).coeffCl(3), classes{l}(index_n_eta_AB(l)).coeffCl(4), classes{l}(index_n_eta_AB(l)).coeffHv(1), classes{l}(index_n_eta_AB(l)).coeffHv(2), classes{l}(index_n_eta_AB(l)).coeffPsat(1), classes{l}(index_n_eta_AB(l)).coeffPsat(2), classes{l}(index_n_eta_AB(l)).coeffPsat(3), classes{l}(index_n_eta_AB(l)).coeffPsat(4), classes{l}(index_n_eta_AB(l)).coeffPsat(5), classes{l}(index_n_eta_AB(l)).coeffSigma(1), classes{l}(index_n_eta_AB(l)).coeffSigma(2), 101325);
            Wl(l) = classes{l}(index_n_eta_AB(l)).molWeight;
            numerator(l) = molFrac_AB(p,l)*classes{l}(index_n_eta_AB(l)).molWeight;
            DCN_l(l) = classes{l}(index_n_eta_AB(l)).DCN;

        end

        y_l = mol2mass(molFrac_AB(p,:), Wl);

        rho_mixture = sum(numerator, 2)./sum(numerator./ rho_l, 2);
        Vl = rho_mixture*y_l'./rho_l';
        y_AB(p) = sum(Vl.*DCN_l');

    elseif strcmp(variable_name, 'flash')
        flash_l = zeros(1,numComponents);
        BI_flash_l = zeros(1,numComponents); % blending index for flash point from Gary and Handwerk (2001)
        Wl = zeros(1, numComponents);

        for l = 1:numComponents

            Wl(l) = classes{l}(index_n_eta_AB(l)).molWeight;
            flash_l(l) = classes{l}(index_n_eta_AB(l)).Tf;
            flash_F = (flash_l(l)-273.15)*9./5+32;
            BI_flash_l(l) = 51708*exp((log(flash_F)-2.6287)^2/(-0.91725));

        end

        y_l = mol2mass(molFrac_AB(p,:), Wl);

        BI_blend = sum(y_l.*BI_flash_l);

        y_AB_F = exp(((-0.91725)*log(BI_blend/51708))^0.5+2.6287); % mixture flash point [°F]
        y_AB(p) = (y_AB_F-32)*5./9+273.15; % mixture flash point [K]

    elseif strcmp(variable_name, 'freezing')
        freezing_l = zeros(1,numComponents);
        BI_freezing_l = zeros(1,numComponents); % blending index for freezing point from Gary and Handwerk (2001)
        Wl = zeros(1, numComponents);
        rho_l = zeros(1,numComponents);
        numerator = zeros(1,numComponents);

        for l = 1:numComponents

            rho_l(l) = ThermophysicalProperties_SingleLiquid('rho', 300, classes{l}(index_n_eta_AB(l)).molWeight, classes{l}(index_n_eta_AB(l)).Tc, classes{l}(index_n_eta_AB(l)).coeffRho(1), classes{l}(index_n_eta_AB(l)).coeffRho(2), classes{l}(index_n_eta_AB(l)).coeffRho(3), classes{l}(index_n_eta_AB(l)).coeffMu(1), classes{l}(index_n_eta_AB(l)).coeffMu(2), classes{l}(index_n_eta_AB(l)).coeffMu(3), classes{l}(index_n_eta_AB(l)).coeffMu(4), classes{l}(index_n_eta_AB(l)).coeffK(1), classes{l}(index_n_eta_AB(l)).coeffK(2), classes{l}(index_n_eta_AB(l)).coeffK(3), classes{l}(index_n_eta_AB(l)).coeffCl(1), classes{l}(index_n_eta_AB(l)).coeffCl(2), classes{l}(index_n_eta_AB(l)).coeffCl(3), classes{l}(index_n_eta_AB(l)).coeffCl(4), classes{l}(index_n_eta_AB(l)).coeffHv(1), classes{l}(index_n_eta_AB(l)).coeffHv(2), classes{l}(index_n_eta_AB(l)).coeffPsat(1), classes{l}(index_n_eta_AB(l)).coeffPsat(2), classes{l}(index_n_eta_AB(l)).coeffPsat(3), classes{l}(index_n_eta_AB(l)).coeffPsat(4), classes{l}(index_n_eta_AB(l)).coeffPsat(5), classes{l}(index_n_eta_AB(l)).coeffSigma(1), classes{l}(index_n_eta_AB(l)).coeffSigma(2), 101325);
            Wl(l) = classes{l}(index_n_eta_AB(l)).molWeight;
            numerator(l) = molFrac_AB(p,l)*classes{l}(index_n_eta_AB(l)).molWeight;
            freezing_l(l) = classes{l}(index_n_eta_AB(l)).Tfz; % freezing point of the l-th surrogate component [K]
            BI_freezing_l(l) = freezing_l(l)^(1/0.05); % blending index as from Hu and Burns (1970), Riazi (2005)

        end

        y_l = mol2mass(molFrac_AB(p,:), Wl); % mass fractions
        rho_mixture = sum(numerator, 2)./sum(numerator./ rho_l, 2);
        Vl = rho_mixture*y_l./rho_l; % volume fractions

        BI_blend = sum(Vl.*BI_freezing_l);

        y_AB(p) = BI_blend^0.05; % mixture freezing point [K]

    elseif strcmp(variable_name, 'LHV')
        LHV_l = zeros(1,numComponents);
        Wl = zeros(1,numComponents);

        for l = 1:numComponents

            Wl(l) = classes{l}(index_n_eta_AB(l)).molWeight;
            LHV_l(l) = classes{l}(index_n_eta_AB(l)).Hc; % LHV of the l-th surrogate component [MJ/kg]

        end

        y_l = mol2mass(molFrac_AB(p,:), Wl); % mass fractions

        y_AB(p) = sum(y_l.*LHV_l);

    end
end

V_2 = abs(mean((y_A - c).*(y_AB - y_B)));
sobol_idx(2) = V_2 / Vy;

% +++ First-order Sobol' indices for branching indices group +++ %
A_B = B;
y_AB = zeros(size(A_B,1), 1);
A_B(:, 2*numComponents:3*numComponents-1) = A(:, 2*numComponents:3*numComponents-1);
for p = 1:size(A_B,1)
    for q = 1:numComponents
        index_n_eta_AB(q) = findIndex_eta(classes{q}, A_B(p, numComponents+q-1), A_B(p, 2*numComponents+q-1));
    end

    molFrac_AB(p,:) = [A_B(p,1:numComponents-1) 1-sum(A_B(p,1:numComponents-1))];

    if strcmp(variable_name, 'molWeight')
        Wl = zeros(1, numComponents);
        for l = 1:numComponents
            Wl(l) = classes{l}(index_n_eta_AB(l)).molWeight;
        end
        y_AB(p) = molWeight(molFrac_AB(p,:), Wl);

    elseif strcmp(variable_name, 'HC')
        nC_AB(p,:) = A_B(p,numComponents:2*numComponents-1);
        nCarbon = round(sum(molFrac_AB(p,:).*nC_AB(p,:)),2);
        nHydrogen = 0.0;
        for l = 1:numComponents
            if strcmp(families(l), 'nparaffins') || strcmp(families(l), 'isoparaffins')
                nHydrogen = nHydrogen + molFrac_AB(p,l)*(2*nC_AB(p,l)+2);
            elseif strcmp(families(l), 'cycloparaffins')
                nHydrogen = nHydrogen + molFrac_AB(p,l)*2*nC_AB(p,l);
            elseif strcmp(families(l), 'dicycloparaffins')
                nHydrogen = nHydrogen + molFrac_AB(p,l)*(2*nC_AB(p,l)-2);
            elseif strcmp(families(l), 'alkylbenzenes')
                nHydrogen = nHydrogen + molFrac_AB(p,l)*(2*nC_AB(p,l)-6);
            elseif strcmp(families(l), 'alkylnaphtalenes')
                nHydrogen = nHydrogen + molFrac_AB(p,l)*(2*nC_AB(p,l)-12);
            elseif strcmp(families(l), 'cycloaromatics')
                nHydrogen = nHydrogen + molFrac_AB(p,l)*(2*nC_AB(p,l)-8);
            end
        end
        nHydrogen = round(nHydrogen, 2);
        y_AB(p) = round(nHydrogen/nCarbon, 3);

    elseif strcmp(variable_name, 'DCN')
        rho_l = zeros(1,numComponents);
        numerator = zeros(1,numComponents);
        DCN_l = zeros(1,numComponents);

        Wl = zeros(1, numComponents);

        for l = 1:numComponents

            rho_l(l) = ThermophysicalProperties_SingleLiquid('rho', 300, classes{l}(index_n_eta_AB(l)).molWeight, classes{l}(index_n_eta_AB(l)).Tc, classes{l}(index_n_eta_AB(l)).coeffRho(1), classes{l}(index_n_eta_AB(l)).coeffRho(2), classes{l}(index_n_eta_AB(l)).coeffRho(3), classes{l}(index_n_eta_AB(l)).coeffMu(1), classes{l}(index_n_eta_AB(l)).coeffMu(2), classes{l}(index_n_eta_AB(l)).coeffMu(3), classes{l}(index_n_eta_AB(l)).coeffMu(4), classes{l}(index_n_eta_AB(l)).coeffK(1), classes{l}(index_n_eta_AB(l)).coeffK(2), classes{l}(index_n_eta_AB(l)).coeffK(3), classes{l}(index_n_eta_AB(l)).coeffCl(1), classes{l}(index_n_eta_AB(l)).coeffCl(2), classes{l}(index_n_eta_AB(l)).coeffCl(3), classes{l}(index_n_eta_AB(l)).coeffCl(4), classes{l}(index_n_eta_AB(l)).coeffHv(1), classes{l}(index_n_eta_AB(l)).coeffHv(2), classes{l}(index_n_eta_AB(l)).coeffPsat(1), classes{l}(index_n_eta_AB(l)).coeffPsat(2), classes{l}(index_n_eta_AB(l)).coeffPsat(3), classes{l}(index_n_eta_AB(l)).coeffPsat(4), classes{l}(index_n_eta_AB(l)).coeffPsat(5), classes{l}(index_n_eta_AB(l)).coeffSigma(1), classes{l}(index_n_eta_AB(l)).coeffSigma(2), 101325);
            Wl(l) = classes{l}(index_n_eta_AB(l)).molWeight;
            numerator(l) = molFrac_AB(p,l)*classes{l}(index_n_eta_AB(l)).molWeight;
            DCN_l(l) = classes{l}(index_n_eta_AB(l)).DCN;

        end

        y_l = mol2mass(molFrac_AB(p,:), Wl);

        rho_mixture = sum(numerator, 2)./sum(numerator./ rho_l, 2);
        Vl = rho_mixture*y_l'./rho_l';
        y_AB(p) = sum(Vl.*DCN_l');

    elseif strcmp(variable_name, 'flash')
        flash_l = zeros(1,numComponents);
        BI_flash_l = zeros(1,numComponents); % blending index for flash point from Gary and Handwerk (2001)
        Wl = zeros(1, numComponents);

        for l = 1:numComponents

            Wl(l) = classes{l}(index_n_eta_AB(l)).molWeight;
            flash_l(l) = classes{l}(index_n_eta_AB(l)).Tf;
            flash_F = (flash_l(l)-273.15)*9./5+32;
            BI_flash_l(l) = 51708*exp((log(flash_F)-2.6287)^2/(-0.91725));

        end

        y_l = mol2mass(molFrac_AB(p,:), Wl);

        BI_blend = sum(y_l.*BI_flash_l);

        y_AB_F = exp(((-0.91725)*log(BI_blend/51708))^0.5+2.6287); % mixture flash point [°F]
        y_AB(p) = (y_AB_F-32)*5./9+273.15; % mixture flash point [K]

    elseif strcmp(variable_name, 'freezing')
        freezing_l = zeros(1,numComponents);
        BI_freezing_l = zeros(1,numComponents); % blending index for freezing point from Gary and Handwerk (2001)
        Wl = zeros(1, numComponents);
        rho_l = zeros(1,numComponents);
        numerator = zeros(1,numComponents);

        for l = 1:numComponents

            rho_l(l) = ThermophysicalProperties_SingleLiquid('rho', 300, classes{l}(index_n_eta_AB(l)).molWeight, classes{l}(index_n_eta_AB(l)).Tc, classes{l}(index_n_eta_AB(l)).coeffRho(1), classes{l}(index_n_eta_AB(l)).coeffRho(2), classes{l}(index_n_eta_AB(l)).coeffRho(3), classes{l}(index_n_eta_AB(l)).coeffMu(1), classes{l}(index_n_eta_AB(l)).coeffMu(2), classes{l}(index_n_eta_AB(l)).coeffMu(3), classes{l}(index_n_eta_AB(l)).coeffMu(4), classes{l}(index_n_eta_AB(l)).coeffK(1), classes{l}(index_n_eta_AB(l)).coeffK(2), classes{l}(index_n_eta_AB(l)).coeffK(3), classes{l}(index_n_eta_AB(l)).coeffCl(1), classes{l}(index_n_eta_AB(l)).coeffCl(2), classes{l}(index_n_eta_AB(l)).coeffCl(3), classes{l}(index_n_eta_AB(l)).coeffCl(4), classes{l}(index_n_eta_AB(l)).coeffHv(1), classes{l}(index_n_eta_AB(l)).coeffHv(2), classes{l}(index_n_eta_AB(l)).coeffPsat(1), classes{l}(index_n_eta_AB(l)).coeffPsat(2), classes{l}(index_n_eta_AB(l)).coeffPsat(3), classes{l}(index_n_eta_AB(l)).coeffPsat(4), classes{l}(index_n_eta_AB(l)).coeffPsat(5), classes{l}(index_n_eta_AB(l)).coeffSigma(1), classes{l}(index_n_eta_AB(l)).coeffSigma(2), 101325);
            Wl(l) = classes{l}(index_n_eta_AB(l)).molWeight;
            numerator(l) = molFrac_AB(p,l)*classes{l}(index_n_eta_AB(l)).molWeight;
            freezing_l(l) = classes{l}(index_n_eta_AB(l)).Tfz; % freezing point of the l-th surrogate component [K]
            BI_freezing_l(l) = freezing_l(l)^(1/0.05); % blending index as from Hu and Burns (1970), Riazi (2005)

        end

        y_l = mol2mass(molFrac_AB(p,:), Wl); % mass fractions
        rho_mixture = sum(numerator, 2)./sum(numerator./ rho_l, 2);
        Vl = rho_mixture*y_l./rho_l; % volume fractions

        BI_blend = sum(Vl.*BI_freezing_l);

        y_AB(p) = BI_blend^0.05; % mixture freezing point [K]

    elseif strcmp(variable_name, 'LHV')
        LHV_l = zeros(1,numComponents);
        Wl = zeros(1,numComponents);

        for l = 1:numComponents

            Wl(l) = classes{l}(index_n_eta_AB(l)).molWeight;
            LHV_l(l) = classes{l}(index_n_eta_AB(l)).Hc; % LHV of the l-th surrogate component [MJ/kg]

        end

        y_l = mol2mass(molFrac_AB(p,:), Wl); % mass fractions

        y_AB(p) = sum(y_l.*LHV_l);

    end

end

V_3 = abs(mean((y_A - c).*(y_AB - y_B)));
sobol_idx(3) = V_3 / Vy;

end