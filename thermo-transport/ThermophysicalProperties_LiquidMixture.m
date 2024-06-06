function phi = ThermophysicalProperties_LiquidMixture(property, T, X, classes, index_n_eta, pressure)

% DESCRIPTION:
% The ThermophysicalProperties_LiquidMixture function calculates the value of a thermophysical
% property for a liquid mixture at given temperature, pressure, and composition.

% Auxiliary parameters:
% Nc: number of surrogate mixture components  [-]

% Inputs:
% 1) property   : string denoting the thermophysical property under consideration
% 2) T          : temperature  [K]
% 3) X          : (1 x Nc-1) array containing the molar fractions of the surrogate
% components in the liquid phase before the distillation process begins  [-]
% 4) classes    : (1 x Nc) cell array, with the i-th cell containing a structure
% array related to the i-th hydrocarbon family, range of number of carbon
% atoms, and range of topochemical atom index
% 5) index_n_eta: (1 x Nc) array, with the i-th element representing the
% index of the species in the i-th cell of "classes" displaying the
% closest topochemical atom index to the index currently investigated for
% the i-th surrogate mixture component
% 6) pressure   : pressure  [Pa]

% Outputs:
% 1) phi: value of the thermophysical property under consideration

% ------------------------------------------------------------------------%
% Contributors/Copyright
% 2024 Jacopo Liberatori, jacopo.liberatori@uniroma1.it
% 2024 Davide Cavalieri,  davide.cavalieri@uniroma1.it
% Department of Mechanical and Aerospace Engineering (DIMA)
% Sapienza University of Rome
% ------------------------------------------------------------------------%

if strcmp(property,'rho')

    numerator = zeros(numel(classes),1);
    rhoi = zeros(numel(classes),1);
    for i = 1:numel(classes)
        if i < numel(classes)
            numerator(i) = X(i)*classes{i}(index_n_eta(i)).molWeight;
        else
            numerator(i) = (1-sum(X))*classes{i}(index_n_eta(i)).molWeight;
        end
        rhoi(i) = ThermophysicalProperties_SingleLiquid(property, T, classes{i}(index_n_eta(i)).molWeight, classes{i}(index_n_eta(i)).Tc, classes{i}(index_n_eta(i)).coeffRho(1), classes{i}(index_n_eta(i)).coeffRho(2), classes{i}(index_n_eta(i)).coeffRho(3), classes{i}(index_n_eta(i)).coeffMu(1), classes{i}(index_n_eta(i)).coeffMu(2), classes{i}(index_n_eta(i)).coeffMu(3), classes{i}(index_n_eta(i)).coeffMu(4), classes{i}(index_n_eta(i)).coeffK(1), classes{i}(index_n_eta(i)).coeffK(2), classes{i}(index_n_eta(i)).coeffK(3), classes{i}(index_n_eta(i)).coeffCl(1), classes{i}(index_n_eta(i)).coeffCl(2), classes{i}(index_n_eta(i)).coeffCl(3), classes{i}(index_n_eta(i)).coeffCl(4), classes{i}(index_n_eta(i)).coeffHv(1), classes{i}(index_n_eta(i)).coeffHv(2), classes{i}(index_n_eta(i)).coeffPsat(1), classes{i}(index_n_eta(i)).coeffPsat(2), classes{i}(index_n_eta(i)).coeffPsat(3), classes{i}(index_n_eta(i)).coeffPsat(4), classes{i}(index_n_eta(i)).coeffPsat(5), classes{i}(index_n_eta(i)).coeffSigma(1), classes{i}(index_n_eta(i)).coeffSigma(2), pressure);
    end

    phi = sum(numerator)/sum(numerator./rhoi);

elseif strcmp(property,'mu')

    mui = zeros(numel(classes),1);
    for i = 1:numel(classes)
        if i < numel(classes)
            mui(i) = X(i)*log(ThermophysicalProperties_SingleLiquid(property, T, classes{i}(index_n_eta(i)).molWeight, classes{i}(index_n_eta(i)).Tc, classes{i}(index_n_eta(i)).coeffRho(1), classes{i}(index_n_eta(i)).coeffRho(2), classes{i}(index_n_eta(i)).coeffRho(3), classes{i}(index_n_eta(i)).coeffMu(1), classes{i}(index_n_eta(i)).coeffMu(2), classes{i}(index_n_eta(i)).coeffMu(3), classes{i}(index_n_eta(i)).coeffMu(4), classes{i}(index_n_eta(i)).coeffK(1), classes{i}(index_n_eta(i)).coeffK(2), classes{i}(index_n_eta(i)).coeffK(3), classes{i}(index_n_eta(i)).coeffCl(1), classes{i}(index_n_eta(i)).coeffCl(2), classes{i}(index_n_eta(i)).coeffCl(3), classes{i}(index_n_eta(i)).coeffCl(4), classes{i}(index_n_eta(i)).coeffHv(1), classes{i}(index_n_eta(i)).coeffHv(2), classes{i}(index_n_eta(i)).coeffPsat(1), classes{i}(index_n_eta(i)).coeffPsat(2), classes{i}(index_n_eta(i)).coeffPsat(3), classes{i}(index_n_eta(i)).coeffPsat(4), classes{i}(index_n_eta(i)).coeffPsat(5), classes{i}(index_n_eta(i)).coeffSigma(1), classes{i}(index_n_eta(i)).coeffSigma(2), pressure));
        else
            mui(i) = (1-sum(X))*log(ThermophysicalProperties_SingleLiquid(property, T, classes{i}(index_n_eta(i)).molWeight, classes{i}(index_n_eta(i)).Tc, classes{i}(index_n_eta(i)).coeffRho(1), classes{i}(index_n_eta(i)).coeffRho(2), classes{i}(index_n_eta(i)).coeffRho(3), classes{i}(index_n_eta(i)).coeffMu(1), classes{i}(index_n_eta(i)).coeffMu(2), classes{i}(index_n_eta(i)).coeffMu(3), classes{i}(index_n_eta(i)).coeffMu(4), classes{i}(index_n_eta(i)).coeffK(1), classes{i}(index_n_eta(i)).coeffK(2), classes{i}(index_n_eta(i)).coeffK(3), classes{i}(index_n_eta(i)).coeffCl(1), classes{i}(index_n_eta(i)).coeffCl(2), classes{i}(index_n_eta(i)).coeffCl(3), classes{i}(index_n_eta(i)).coeffCl(4), classes{i}(index_n_eta(i)).coeffHv(1), classes{i}(index_n_eta(i)).coeffHv(2), classes{i}(index_n_eta(i)).coeffPsat(1), classes{i}(index_n_eta(i)).coeffPsat(2), classes{i}(index_n_eta(i)).coeffPsat(3), classes{i}(index_n_eta(i)).coeffPsat(4), classes{i}(index_n_eta(i)).coeffPsat(5), classes{i}(index_n_eta(i)).coeffSigma(1), classes{i}(index_n_eta(i)).coeffSigma(2), pressure));
        end
    end

    phi = exp(sum(mui));

elseif strcmp(property,'nu')

    numerator = zeros(numel(classes),1);
    rhoi = zeros(numel(classes),1);
    for i = 1:numel(classes)
        if i < numel(classes)
            numerator(i) = X(i)*classes{i}(index_n_eta(i)).molWeight;
        else
            numerator(i) = (1-sum(X))*classes{i}(index_n_eta(i)).molWeight;
        end
        rhoi(i) = ThermophysicalProperties_SingleLiquid('rho', T, classes{i}(index_n_eta(i)).molWeight, classes{i}(index_n_eta(i)).Tc, classes{i}(index_n_eta(i)).coeffRho(1), classes{i}(index_n_eta(i)).coeffRho(2), classes{i}(index_n_eta(i)).coeffRho(3), classes{i}(index_n_eta(i)).coeffMu(1), classes{i}(index_n_eta(i)).coeffMu(2), classes{i}(index_n_eta(i)).coeffMu(3), classes{i}(index_n_eta(i)).coeffMu(4), classes{i}(index_n_eta(i)).coeffK(1), classes{i}(index_n_eta(i)).coeffK(2), classes{i}(index_n_eta(i)).coeffK(3), classes{i}(index_n_eta(i)).coeffCl(1), classes{i}(index_n_eta(i)).coeffCl(2), classes{i}(index_n_eta(i)).coeffCl(3), classes{i}(index_n_eta(i)).coeffCl(4), classes{i}(index_n_eta(i)).coeffHv(1), classes{i}(index_n_eta(i)).coeffHv(2), classes{i}(index_n_eta(i)).coeffPsat(1), classes{i}(index_n_eta(i)).coeffPsat(2), classes{i}(index_n_eta(i)).coeffPsat(3), classes{i}(index_n_eta(i)).coeffPsat(4), classes{i}(index_n_eta(i)).coeffPsat(5), classes{i}(index_n_eta(i)).coeffSigma(1), classes{i}(index_n_eta(i)).coeffSigma(2), pressure);
    end

    rhoMixture = sum(numerator)/sum(numerator./rhoi);

    mui = zeros(numel(classes),1);
    for i = 1:numel(classes)
        if i < numel(classes)
            mui(i) = X(i)*log(ThermophysicalProperties_SingleLiquid('mu', T, classes{i}(index_n_eta(i)).molWeight, classes{i}(index_n_eta(i)).Tc, classes{i}(index_n_eta(i)).coeffRho(1), classes{i}(index_n_eta(i)).coeffRho(2), classes{i}(index_n_eta(i)).coeffRho(3), classes{i}(index_n_eta(i)).coeffMu(1), classes{i}(index_n_eta(i)).coeffMu(2), classes{i}(index_n_eta(i)).coeffMu(3), classes{i}(index_n_eta(i)).coeffMu(4), classes{i}(index_n_eta(i)).coeffK(1), classes{i}(index_n_eta(i)).coeffK(2), classes{i}(index_n_eta(i)).coeffK(3), classes{i}(index_n_eta(i)).coeffCl(1), classes{i}(index_n_eta(i)).coeffCl(2), classes{i}(index_n_eta(i)).coeffCl(3), classes{i}(index_n_eta(i)).coeffCl(4), classes{i}(index_n_eta(i)).coeffHv(1), classes{i}(index_n_eta(i)).coeffHv(2), classes{i}(index_n_eta(i)).coeffPsat(1), classes{i}(index_n_eta(i)).coeffPsat(2), classes{i}(index_n_eta(i)).coeffPsat(3), classes{i}(index_n_eta(i)).coeffPsat(4), classes{i}(index_n_eta(i)).coeffPsat(5), classes{i}(index_n_eta(i)).coeffSigma(1), classes{i}(index_n_eta(i)).coeffSigma(2), pressure));
        else
            mui(i) = (1-sum(X))*log(ThermophysicalProperties_SingleLiquid('mu', T, classes{i}(index_n_eta(i)).molWeight, classes{i}(index_n_eta(i)).Tc, classes{i}(index_n_eta(i)).coeffRho(1), classes{i}(index_n_eta(i)).coeffRho(2), classes{i}(index_n_eta(i)).coeffRho(3), classes{i}(index_n_eta(i)).coeffMu(1), classes{i}(index_n_eta(i)).coeffMu(2), classes{i}(index_n_eta(i)).coeffMu(3), classes{i}(index_n_eta(i)).coeffMu(4), classes{i}(index_n_eta(i)).coeffK(1), classes{i}(index_n_eta(i)).coeffK(2), classes{i}(index_n_eta(i)).coeffK(3), classes{i}(index_n_eta(i)).coeffCl(1), classes{i}(index_n_eta(i)).coeffCl(2), classes{i}(index_n_eta(i)).coeffCl(3), classes{i}(index_n_eta(i)).coeffCl(4), classes{i}(index_n_eta(i)).coeffHv(1), classes{i}(index_n_eta(i)).coeffHv(2), classes{i}(index_n_eta(i)).coeffPsat(1), classes{i}(index_n_eta(i)).coeffPsat(2), classes{i}(index_n_eta(i)).coeffPsat(3), classes{i}(index_n_eta(i)).coeffPsat(4), classes{i}(index_n_eta(i)).coeffPsat(5), classes{i}(index_n_eta(i)).coeffSigma(1), classes{i}(index_n_eta(i)).coeffSigma(2), pressure));
        end
    end

    muMixture = exp(sum(mui));

    phi = muMixture/rhoMixture;

elseif strcmp(property,'kappa')

    % Initialize array of superficial volume fractions
    psi_array = zeros(numel(classes),1);
    for i = 1:numel(classes)
        if i < numel(classes)
            psi_array(i) = X(i)*classes{i}(index_n_eta(i)).molWeight/ThermophysicalProperties_SingleLiquid('rho', T, classes{i}(index_n_eta(i)).molWeight, classes{i}(index_n_eta(i)).Tc, classes{i}(index_n_eta(i)).coeffRho(1), classes{i}(index_n_eta(i)).coeffRho(2), classes{i}(index_n_eta(i)).coeffRho(3), classes{i}(index_n_eta(i)).coeffMu(1), classes{i}(index_n_eta(i)).coeffMu(2), classes{i}(index_n_eta(i)).coeffMu(3), classes{i}(index_n_eta(i)).coeffMu(4), classes{i}(index_n_eta(i)).coeffK(1), classes{i}(index_n_eta(i)).coeffK(2), classes{i}(index_n_eta(i)).coeffK(3), classes{i}(index_n_eta(i)).coeffCl(1), classes{i}(index_n_eta(i)).coeffCl(2), classes{i}(index_n_eta(i)).coeffCl(3), classes{i}(index_n_eta(i)).coeffCl(4), classes{i}(index_n_eta(i)).coeffHv(1), classes{i}(index_n_eta(i)).coeffHv(2), classes{i}(index_n_eta(i)).coeffPsat(1), classes{i}(index_n_eta(i)).coeffPsat(2), classes{i}(index_n_eta(i)).coeffPsat(3), classes{i}(index_n_eta(i)).coeffPsat(4), classes{i}(index_n_eta(i)).coeffPsat(5), classes{i}(index_n_eta(i)).coeffSigma(1), classes{i}(index_n_eta(i)).coeffSigma(2), pressure);
        else
            psi_array(i) = (1-sum(X))*classes{i}(index_n_eta(i)).molWeight/ThermophysicalProperties_SingleLiquid('rho', T, classes{i}(index_n_eta(i)).molWeight, classes{i}(index_n_eta(i)).Tc, classes{i}(index_n_eta(i)).coeffRho(1), classes{i}(index_n_eta(i)).coeffRho(2), classes{i}(index_n_eta(i)).coeffRho(3), classes{i}(index_n_eta(i)).coeffMu(1), classes{i}(index_n_eta(i)).coeffMu(2), classes{i}(index_n_eta(i)).coeffMu(3), classes{i}(index_n_eta(i)).coeffMu(4), classes{i}(index_n_eta(i)).coeffK(1), classes{i}(index_n_eta(i)).coeffK(2), classes{i}(index_n_eta(i)).coeffK(3), classes{i}(index_n_eta(i)).coeffCl(1), classes{i}(index_n_eta(i)).coeffCl(2), classes{i}(index_n_eta(i)).coeffCl(3), classes{i}(index_n_eta(i)).coeffCl(4), classes{i}(index_n_eta(i)).coeffHv(1), classes{i}(index_n_eta(i)).coeffHv(2), classes{i}(index_n_eta(i)).coeffPsat(1), classes{i}(index_n_eta(i)).coeffPsat(2), classes{i}(index_n_eta(i)).coeffPsat(3), classes{i}(index_n_eta(i)).coeffPsat(4), classes{i}(index_n_eta(i)).coeffPsat(5), classes{i}(index_n_eta(i)).coeffSigma(1), classes{i}(index_n_eta(i)).coeffSigma(2), pressure);
        end
    end

    pSum = 0.0;
    for i=1:numel(psi_array)
        pSum = pSum + psi_array(i);
    end

    psi_array = psi_array/pSum;

    phi = 0.0;

    for j=1:numel(psi_array)
        kappa_j = ThermophysicalProperties_SingleLiquid(property, T, classes{j}(index_n_eta(j)).molWeight, classes{j}(index_n_eta(j)).Tc, classes{j}(index_n_eta(j)).coeffRho(1), classes{j}(index_n_eta(j)).coeffRho(2), classes{j}(index_n_eta(j)).coeffRho(3), classes{j}(index_n_eta(j)).coeffMu(1), classes{j}(index_n_eta(j)).coeffMu(2), classes{j}(index_n_eta(j)).coeffMu(3), classes{j}(index_n_eta(j)).coeffMu(4), classes{j}(index_n_eta(j)).coeffK(1), classes{j}(index_n_eta(j)).coeffK(2), classes{j}(index_n_eta(j)).coeffK(3), classes{j}(index_n_eta(j)).coeffCl(1), classes{j}(index_n_eta(j)).coeffCl(2), classes{j}(index_n_eta(j)).coeffCl(3), classes{j}(index_n_eta(j)).coeffCl(4), classes{j}(index_n_eta(j)).coeffHv(1), classes{j}(index_n_eta(j)).coeffHv(2), classes{j}(index_n_eta(j)).coeffPsat(1), classes{j}(index_n_eta(j)).coeffPsat(2), classes{j}(index_n_eta(j)).coeffPsat(3), classes{j}(index_n_eta(j)).coeffPsat(4), classes{j}(index_n_eta(j)).coeffPsat(5), classes{j}(index_n_eta(j)).coeffSigma(1), classes{j}(index_n_eta(j)).coeffSigma(2), pressure);
        for k=1:numel(psi_array)
            kappa_k = ThermophysicalProperties_SingleLiquid(property, T, classes{k}(index_n_eta(k)).molWeight, classes{k}(index_n_eta(k)).Tc, classes{k}(index_n_eta(k)).coeffRho(1), classes{k}(index_n_eta(k)).coeffRho(2), classes{k}(index_n_eta(k)).coeffRho(3), classes{k}(index_n_eta(k)).coeffMu(1), classes{k}(index_n_eta(k)).coeffMu(2), classes{k}(index_n_eta(k)).coeffMu(3), classes{k}(index_n_eta(k)).coeffMu(4), classes{k}(index_n_eta(k)).coeffK(1), classes{k}(index_n_eta(k)).coeffK(2), classes{k}(index_n_eta(k)).coeffK(3), classes{k}(index_n_eta(k)).coeffCl(1), classes{k}(index_n_eta(k)).coeffCl(2), classes{k}(index_n_eta(k)).coeffCl(3), classes{k}(index_n_eta(k)).coeffCl(4), classes{k}(index_n_eta(k)).coeffHv(1), classes{k}(index_n_eta(k)).coeffHv(2), classes{k}(index_n_eta(k)).coeffPsat(1), classes{k}(index_n_eta(k)).coeffPsat(2), classes{k}(index_n_eta(k)).coeffPsat(3), classes{k}(index_n_eta(k)).coeffPsat(4), classes{k}(index_n_eta(k)).coeffPsat(5), classes{k}(index_n_eta(k)).coeffSigma(1), classes{k}(index_n_eta(k)).coeffSigma(2), pressure);
            Kij = 2.0/(1.0/kappa_j + 1.0/kappa_k);
            phi = phi + psi_array(j)*psi_array(k)*Kij;
        end
    end

elseif strcmp(property,'specificHeat')

    numerator = zeros(numel(classes),1);
    denominator = zeros(numel(classes),1);
    for i = 1:numel(classes)
        if i < numel(classes)
            numerator(i) = X(i)*classes{i}(index_n_eta(i)).molWeight*ThermophysicalProperties_SingleLiquid(property, T, classes{i}(index_n_eta(i)).molWeight, classes{i}(index_n_eta(i)).Tc, classes{i}(index_n_eta(i)).coeffRho(1), classes{i}(index_n_eta(i)).coeffRho(2), classes{i}(index_n_eta(i)).coeffRho(3), classes{i}(index_n_eta(i)).coeffMu(1), classes{i}(index_n_eta(i)).coeffMu(2), classes{i}(index_n_eta(i)).coeffMu(3), classes{i}(index_n_eta(i)).coeffMu(4), classes{i}(index_n_eta(i)).coeffK(1), classes{i}(index_n_eta(i)).coeffK(2), classes{i}(index_n_eta(i)).coeffK(3), classes{i}(index_n_eta(i)).coeffCl(1), classes{i}(index_n_eta(i)).coeffCl(2), classes{i}(index_n_eta(i)).coeffCl(3), classes{i}(index_n_eta(i)).coeffCl(4), classes{i}(index_n_eta(i)).coeffHv(1), classes{i}(index_n_eta(i)).coeffHv(2), classes{i}(index_n_eta(i)).coeffPsat(1), classes{i}(index_n_eta(i)).coeffPsat(2), classes{i}(index_n_eta(i)).coeffPsat(3), classes{i}(index_n_eta(i)).coeffPsat(4), classes{i}(index_n_eta(i)).coeffPsat(5), classes{i}(index_n_eta(i)).coeffSigma(1), classes{i}(index_n_eta(i)).coeffSigma(2), pressure);
            denominator(i) = X(i)*classes{i}(index_n_eta(i)).molWeight;
        else
            numerator(i) = (1-sum(X))*classes{i}(index_n_eta(i)).molWeight*ThermophysicalProperties_SingleLiquid(property, T, classes{i}(index_n_eta(i)).molWeight, classes{i}(index_n_eta(i)).Tc, classes{i}(index_n_eta(i)).coeffRho(1), classes{i}(index_n_eta(i)).coeffRho(2), classes{i}(index_n_eta(i)).coeffRho(3), classes{i}(index_n_eta(i)).coeffMu(1), classes{i}(index_n_eta(i)).coeffMu(2), classes{i}(index_n_eta(i)).coeffMu(3), classes{i}(index_n_eta(i)).coeffMu(4), classes{i}(index_n_eta(i)).coeffK(1), classes{i}(index_n_eta(i)).coeffK(2), classes{i}(index_n_eta(i)).coeffK(3), classes{i}(index_n_eta(i)).coeffCl(1), classes{i}(index_n_eta(i)).coeffCl(2), classes{i}(index_n_eta(i)).coeffCl(3), classes{i}(index_n_eta(i)).coeffCl(4), classes{i}(index_n_eta(i)).coeffHv(1), classes{i}(index_n_eta(i)).coeffHv(2), classes{i}(index_n_eta(i)).coeffPsat(1), classes{i}(index_n_eta(i)).coeffPsat(2), classes{i}(index_n_eta(i)).coeffPsat(3), classes{i}(index_n_eta(i)).coeffPsat(4), classes{i}(index_n_eta(i)).coeffPsat(5), classes{i}(index_n_eta(i)).coeffSigma(1), classes{i}(index_n_eta(i)).coeffSigma(2), pressure);
            denominator(i) = (1-sum(X))*classes{i}(index_n_eta(i)).molWeight;
        end
    end

    phi = sum(numerator)/sum(denominator);

elseif strcmp(property,'latentHeat')

    hli = zeros(numel(classes),1);
    Yi = zeros(numel(classes),1);
    for i = 1:numel(classes)
        if i < numel(classes)
            Yi(i) = X(i)*classes{i}(index_n_eta(i)).molWeight;
            hli(i) = Yi(i)*ThermophysicalProperties_SingleLiquid(property, T, classes{i}(index_n_eta(i)).molWeight, classes{i}(index_n_eta(i)).Tc, classes{i}(index_n_eta(i)).coeffRho(1), classes{i}(index_n_eta(i)).coeffRho(2), classes{i}(index_n_eta(i)).coeffRho(3), classes{i}(index_n_eta(i)).coeffMu(1), classes{i}(index_n_eta(i)).coeffMu(2), classes{i}(index_n_eta(i)).coeffMu(3), classes{i}(index_n_eta(i)).coeffMu(4), classes{i}(index_n_eta(i)).coeffK(1), classes{i}(index_n_eta(i)).coeffK(2), classes{i}(index_n_eta(i)).coeffK(3), classes{i}(index_n_eta(i)).coeffCl(1), classes{i}(index_n_eta(i)).coeffCl(2), classes{i}(index_n_eta(i)).coeffCl(3), classes{i}(index_n_eta(i)).coeffCl(4), classes{i}(index_n_eta(i)).coeffHv(1), classes{i}(index_n_eta(i)).coeffHv(2), classes{i}(index_n_eta(i)).coeffPsat(1), classes{i}(index_n_eta(i)).coeffPsat(2), classes{i}(index_n_eta(i)).coeffPsat(3), classes{i}(index_n_eta(i)).coeffPsat(4), classes{i}(index_n_eta(i)).coeffPsat(5), classes{i}(index_n_eta(i)).coeffSigma(1), classes{i}(index_n_eta(i)).coeffSigma(2), pressure);
        else
            Yi(i) = (1-sum(X))*classes{i}(index_n_eta(i)).molWeight;
            hli(i) = Yi(i)*ThermophysicalProperties_SingleLiquid(property, T, classes{i}(index_n_eta(i)).molWeight, classes{i}(index_n_eta(i)).Tc, classes{i}(index_n_eta(i)).coeffRho(1), classes{i}(index_n_eta(i)).coeffRho(2), classes{i}(index_n_eta(i)).coeffRho(3), classes{i}(index_n_eta(i)).coeffMu(1), classes{i}(index_n_eta(i)).coeffMu(2), classes{i}(index_n_eta(i)).coeffMu(3), classes{i}(index_n_eta(i)).coeffMu(4), classes{i}(index_n_eta(i)).coeffK(1), classes{i}(index_n_eta(i)).coeffK(2), classes{i}(index_n_eta(i)).coeffK(3), classes{i}(index_n_eta(i)).coeffCl(1), classes{i}(index_n_eta(i)).coeffCl(2), classes{i}(index_n_eta(i)).coeffCl(3), classes{i}(index_n_eta(i)).coeffCl(4), classes{i}(index_n_eta(i)).coeffHv(1), classes{i}(index_n_eta(i)).coeffHv(2), classes{i}(index_n_eta(i)).coeffPsat(1), classes{i}(index_n_eta(i)).coeffPsat(2), classes{i}(index_n_eta(i)).coeffPsat(3), classes{i}(index_n_eta(i)).coeffPsat(4), classes{i}(index_n_eta(i)).coeffPsat(5), classes{i}(index_n_eta(i)).coeffSigma(1), classes{i}(index_n_eta(i)).coeffSigma(2), pressure);
        end
    end

    phi = sum(hli)/sum(Yi);

elseif strcmp(property,'vaporPressure')

    pSati = zeros(numel(classes),1);
    Yi = zeros(numel(classes),1);
    for i = 1:numel(classes)
        if i < numel(classes)
            Yi(i) = X(i)*classes{i}(index_n_eta(i)).molWeight;
            pSati(i) = Yi(i)*ThermophysicalProperties_SingleLiquid(property, T, classes{i}(index_n_eta(i)).molWeight, classes{i}(index_n_eta(i)).Tc, classes{i}(index_n_eta(i)).coeffRho(1), classes{i}(index_n_eta(i)).coeffRho(2), classes{i}(index_n_eta(i)).coeffRho(3), classes{i}(index_n_eta(i)).coeffMu(1), classes{i}(index_n_eta(i)).coeffMu(2), classes{i}(index_n_eta(i)).coeffMu(3), classes{i}(index_n_eta(i)).coeffMu(4), classes{i}(index_n_eta(i)).coeffK(1), classes{i}(index_n_eta(i)).coeffK(2), classes{i}(index_n_eta(i)).coeffK(3), classes{i}(index_n_eta(i)).coeffCl(1), classes{i}(index_n_eta(i)).coeffCl(2), classes{i}(index_n_eta(i)).coeffCl(3), classes{i}(index_n_eta(i)).coeffCl(4), classes{i}(index_n_eta(i)).coeffHv(1), classes{i}(index_n_eta(i)).coeffHv(2), classes{i}(index_n_eta(i)).coeffPsat(1), classes{i}(index_n_eta(i)).coeffPsat(2), classes{i}(index_n_eta(i)).coeffPsat(3), classes{i}(index_n_eta(i)).coeffPsat(4), classes{i}(index_n_eta(i)).coeffPsat(5), classes{i}(index_n_eta(i)).coeffSigma(1), classes{i}(index_n_eta(i)).coeffSigma(2), pressure);
        else
            Yi(i) = (1-sum(X))*classes{i}(index_n_eta(i)).molWeight;
            pSati(i) = Yi(i)*ThermophysicalProperties_SingleLiquid(property, T, classes{i}(index_n_eta(i)).molWeight, classes{i}(index_n_eta(i)).Tc, classes{i}(index_n_eta(i)).coeffRho(1), classes{i}(index_n_eta(i)).coeffRho(2), classes{i}(index_n_eta(i)).coeffRho(3), classes{i}(index_n_eta(i)).coeffMu(1), classes{i}(index_n_eta(i)).coeffMu(2), classes{i}(index_n_eta(i)).coeffMu(3), classes{i}(index_n_eta(i)).coeffMu(4), classes{i}(index_n_eta(i)).coeffK(1), classes{i}(index_n_eta(i)).coeffK(2), classes{i}(index_n_eta(i)).coeffK(3), classes{i}(index_n_eta(i)).coeffCl(1), classes{i}(index_n_eta(i)).coeffCl(2), classes{i}(index_n_eta(i)).coeffCl(3), classes{i}(index_n_eta(i)).coeffCl(4), classes{i}(index_n_eta(i)).coeffHv(1), classes{i}(index_n_eta(i)).coeffHv(2), classes{i}(index_n_eta(i)).coeffPsat(1), classes{i}(index_n_eta(i)).coeffPsat(2), classes{i}(index_n_eta(i)).coeffPsat(3), classes{i}(index_n_eta(i)).coeffPsat(4), classes{i}(index_n_eta(i)).coeffPsat(5), classes{i}(index_n_eta(i)).coeffSigma(1), classes{i}(index_n_eta(i)).coeffSigma(2), pressure);
        end
    end

    phi = sum(pSati)/sum(Yi);

elseif strcmp(property,'sigma')

    Xsi = zeros(numel(classes),1); % surface molar fractions estimated by Raoult's law
    sigmai = zeros(numel(classes),1);
    XsSum = 0.0;
    for i = 1:numel(classes)
        sigmai(i) = ThermophysicalProperties_SingleLiquid(property, T, classes{i}(index_n_eta(i)).molWeight, classes{i}(index_n_eta(i)).Tc, classes{i}(index_n_eta(i)).coeffRho(1), classes{i}(index_n_eta(i)).coeffRho(2), classes{i}(index_n_eta(i)).coeffRho(3), classes{i}(index_n_eta(i)).coeffMu(1), classes{i}(index_n_eta(i)).coeffMu(2), classes{i}(index_n_eta(i)).coeffMu(3), classes{i}(index_n_eta(i)).coeffMu(4), classes{i}(index_n_eta(i)).coeffK(1), classes{i}(index_n_eta(i)).coeffK(2), classes{i}(index_n_eta(i)).coeffK(3), classes{i}(index_n_eta(i)).coeffCl(1), classes{i}(index_n_eta(i)).coeffCl(2), classes{i}(index_n_eta(i)).coeffCl(3), classes{i}(index_n_eta(i)).coeffCl(4), classes{i}(index_n_eta(i)).coeffHv(1), classes{i}(index_n_eta(i)).coeffHv(2), classes{i}(index_n_eta(i)).coeffPsat(1), classes{i}(index_n_eta(i)).coeffPsat(2), classes{i}(index_n_eta(i)).coeffPsat(3), classes{i}(index_n_eta(i)).coeffPsat(4), classes{i}(index_n_eta(i)).coeffPsat(5), classes{i}(index_n_eta(i)).coeffSigma(1), classes{i}(index_n_eta(i)).coeffSigma(2), pressure);
        pSati = ThermophysicalProperties_SingleLiquid('vaporPressure', T, classes{i}(index_n_eta(i)).molWeight, classes{i}(index_n_eta(i)).Tc, classes{i}(index_n_eta(i)).coeffRho(1), classes{i}(index_n_eta(i)).coeffRho(2), classes{i}(index_n_eta(i)).coeffRho(3), classes{i}(index_n_eta(i)).coeffMu(1), classes{i}(index_n_eta(i)).coeffMu(2), classes{i}(index_n_eta(i)).coeffMu(3), classes{i}(index_n_eta(i)).coeffMu(4), classes{i}(index_n_eta(i)).coeffK(1), classes{i}(index_n_eta(i)).coeffK(2), classes{i}(index_n_eta(i)).coeffK(3), classes{i}(index_n_eta(i)).coeffCl(1), classes{i}(index_n_eta(i)).coeffCl(2), classes{i}(index_n_eta(i)).coeffCl(3), classes{i}(index_n_eta(i)).coeffCl(4), classes{i}(index_n_eta(i)).coeffHv(1), classes{i}(index_n_eta(i)).coeffHv(2), classes{i}(index_n_eta(i)).coeffPsat(1), classes{i}(index_n_eta(i)).coeffPsat(2), classes{i}(index_n_eta(i)).coeffPsat(3), classes{i}(index_n_eta(i)).coeffPsat(4), classes{i}(index_n_eta(i)).coeffPsat(5), classes{i}(index_n_eta(i)).coeffSigma(1), classes{i}(index_n_eta(i)).coeffSigma(2), pressure);
        if i < numel(classes)
            Xsi(i) = X(i)*pSati/pressure;
        else
            Xsi(i) = (1-sum(X))*pSati/pressure;
        end
        XsSum = XsSum + Xsi(i);
    end

    Xsi = Xsi./XsSum;

    phi = sum(Xsi.*sigmai);

end

end