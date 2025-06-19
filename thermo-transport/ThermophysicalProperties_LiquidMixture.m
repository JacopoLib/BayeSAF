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

molFrac = zeros(numel(classes),1);

for i = 1:numel(classes)
    if i < numel(classes)
        molFrac(i) = X(i);
    else
        molFrac(i) = 1-sum(X);
    end
end

if strcmp(property,'rho')

    rho = zeros(size(T,1),numel(classes));
    numerator = zeros(size(T,1),numel(classes),1);
    
    for i = 1:numel(classes)
        rho(:,i) = ThermophysicalProperties_SingleLiquid(property, T, classes{i}(index_n_eta(i)).molWeight, classes{i}(index_n_eta(i)).Tc, classes{i}(index_n_eta(i)).coeffRho(1), classes{i}(index_n_eta(i)).coeffRho(2), classes{i}(index_n_eta(i)).coeffRho(3), classes{i}(index_n_eta(i)).coeffMu(1), classes{i}(index_n_eta(i)).coeffMu(2), classes{i}(index_n_eta(i)).coeffMu(3), classes{i}(index_n_eta(i)).coeffMu(4), classes{i}(index_n_eta(i)).coeffK(1), classes{i}(index_n_eta(i)).coeffK(2), classes{i}(index_n_eta(i)).coeffK(3), classes{i}(index_n_eta(i)).coeffCl(1), classes{i}(index_n_eta(i)).coeffCl(2), classes{i}(index_n_eta(i)).coeffCl(3), classes{i}(index_n_eta(i)).coeffCl(4), classes{i}(index_n_eta(i)).coeffHv(1), classes{i}(index_n_eta(i)).coeffHv(2), classes{i}(index_n_eta(i)).coeffPsat(1), classes{i}(index_n_eta(i)).coeffPsat(2), classes{i}(index_n_eta(i)).coeffPsat(3), classes{i}(index_n_eta(i)).coeffPsat(4), classes{i}(index_n_eta(i)).coeffPsat(5), classes{i}(index_n_eta(i)).coeffSigma(1), classes{i}(index_n_eta(i)).coeffSigma(2), pressure);
        numerator(:,i) = molFrac(i)*classes{i}(index_n_eta(i)).molWeight;
    end
    
    phi = sum(numerator, 2)./sum(numerator./ rho, 2);

elseif strcmp(property,'mu')

    mu = zeros(size(T,1),numel(classes));
    
    for i = 1:numel(classes)
        mu(:,i) = molFrac(i)*log(ThermophysicalProperties_SingleLiquid(property, T, classes{i}(index_n_eta(i)).molWeight, classes{i}(index_n_eta(i)).Tc, classes{i}(index_n_eta(i)).coeffRho(1), classes{i}(index_n_eta(i)).coeffRho(2), classes{i}(index_n_eta(i)).coeffRho(3), classes{i}(index_n_eta(i)).coeffMu(1), classes{i}(index_n_eta(i)).coeffMu(2), classes{i}(index_n_eta(i)).coeffMu(3), classes{i}(index_n_eta(i)).coeffMu(4), classes{i}(index_n_eta(i)).coeffK(1), classes{i}(index_n_eta(i)).coeffK(2), classes{i}(index_n_eta(i)).coeffK(3), classes{i}(index_n_eta(i)).coeffCl(1), classes{i}(index_n_eta(i)).coeffCl(2), classes{i}(index_n_eta(i)).coeffCl(3), classes{i}(index_n_eta(i)).coeffCl(4), classes{i}(index_n_eta(i)).coeffHv(1), classes{i}(index_n_eta(i)).coeffHv(2), classes{i}(index_n_eta(i)).coeffPsat(1), classes{i}(index_n_eta(i)).coeffPsat(2), classes{i}(index_n_eta(i)).coeffPsat(3), classes{i}(index_n_eta(i)).coeffPsat(4), classes{i}(index_n_eta(i)).coeffPsat(5), classes{i}(index_n_eta(i)).coeffSigma(1), classes{i}(index_n_eta(i)).coeffSigma(2), pressure));
    end
    
    phi = exp(sum(mu,2));

elseif strcmp(property,'nu')

    rho = zeros(size(T,1),numel(classes));
    numerator = zeros(size(T,1),numel(classes),1);
    
    for i = 1:numel(classes)
        rho(:,i) = ThermophysicalProperties_SingleLiquid('rho', T, classes{i}(index_n_eta(i)).molWeight, classes{i}(index_n_eta(i)).Tc, classes{i}(index_n_eta(i)).coeffRho(1), classes{i}(index_n_eta(i)).coeffRho(2), classes{i}(index_n_eta(i)).coeffRho(3), classes{i}(index_n_eta(i)).coeffMu(1), classes{i}(index_n_eta(i)).coeffMu(2), classes{i}(index_n_eta(i)).coeffMu(3), classes{i}(index_n_eta(i)).coeffMu(4), classes{i}(index_n_eta(i)).coeffK(1), classes{i}(index_n_eta(i)).coeffK(2), classes{i}(index_n_eta(i)).coeffK(3), classes{i}(index_n_eta(i)).coeffCl(1), classes{i}(index_n_eta(i)).coeffCl(2), classes{i}(index_n_eta(i)).coeffCl(3), classes{i}(index_n_eta(i)).coeffCl(4), classes{i}(index_n_eta(i)).coeffHv(1), classes{i}(index_n_eta(i)).coeffHv(2), classes{i}(index_n_eta(i)).coeffPsat(1), classes{i}(index_n_eta(i)).coeffPsat(2), classes{i}(index_n_eta(i)).coeffPsat(3), classes{i}(index_n_eta(i)).coeffPsat(4), classes{i}(index_n_eta(i)).coeffPsat(5), classes{i}(index_n_eta(i)).coeffSigma(1), classes{i}(index_n_eta(i)).coeffSigma(2), pressure);
        numerator(:,i) = molFrac(i)*classes{i}(index_n_eta(i)).molWeight;
    end
    
    rhoMixture = sum(numerator, 2)./sum(numerator./ rho, 2);

    mu = zeros(size(T,1),numel(classes));
    
    for i = 1:numel(classes)
        mu(:,i) = molFrac(i)*log(ThermophysicalProperties_SingleLiquid('mu', T, classes{i}(index_n_eta(i)).molWeight, classes{i}(index_n_eta(i)).Tc, classes{i}(index_n_eta(i)).coeffRho(1), classes{i}(index_n_eta(i)).coeffRho(2), classes{i}(index_n_eta(i)).coeffRho(3), classes{i}(index_n_eta(i)).coeffMu(1), classes{i}(index_n_eta(i)).coeffMu(2), classes{i}(index_n_eta(i)).coeffMu(3), classes{i}(index_n_eta(i)).coeffMu(4), classes{i}(index_n_eta(i)).coeffK(1), classes{i}(index_n_eta(i)).coeffK(2), classes{i}(index_n_eta(i)).coeffK(3), classes{i}(index_n_eta(i)).coeffCl(1), classes{i}(index_n_eta(i)).coeffCl(2), classes{i}(index_n_eta(i)).coeffCl(3), classes{i}(index_n_eta(i)).coeffCl(4), classes{i}(index_n_eta(i)).coeffHv(1), classes{i}(index_n_eta(i)).coeffHv(2), classes{i}(index_n_eta(i)).coeffPsat(1), classes{i}(index_n_eta(i)).coeffPsat(2), classes{i}(index_n_eta(i)).coeffPsat(3), classes{i}(index_n_eta(i)).coeffPsat(4), classes{i}(index_n_eta(i)).coeffPsat(5), classes{i}(index_n_eta(i)).coeffSigma(1), classes{i}(index_n_eta(i)).coeffSigma(2), pressure));
    end
    
    muMixture = exp(sum(mu,2));

    phi = muMixture./rhoMixture;

elseif strcmp(property,'kappa')

    for i = 1:numel(T)
        % Initialize array of superficial volume fractions
        psi_array = zeros(numel(classes),1);
        for j = 1:numel(classes)
            psi_array(j) = molFrac(j)*classes{j}(index_n_eta(j)).molWeight/ThermophysicalProperties_SingleLiquid('rho', T(i), classes{j}(index_n_eta(j)).molWeight, classes{j}(index_n_eta(j)).Tc, classes{j}(index_n_eta(j)).coeffRho(1), classes{j}(index_n_eta(j)).coeffRho(2), classes{j}(index_n_eta(j)).coeffRho(3), classes{j}(index_n_eta(j)).coeffMu(1), classes{j}(index_n_eta(j)).coeffMu(2), classes{j}(index_n_eta(j)).coeffMu(3), classes{j}(index_n_eta(j)).coeffMu(4), classes{j}(index_n_eta(j)).coeffK(1), classes{j}(index_n_eta(j)).coeffK(2), classes{j}(index_n_eta(j)).coeffK(3), classes{j}(index_n_eta(j)).coeffCl(1), classes{j}(index_n_eta(j)).coeffCl(2), classes{j}(index_n_eta(j)).coeffCl(3), classes{j}(index_n_eta(j)).coeffCl(4), classes{j}(index_n_eta(j)).coeffHv(1), classes{j}(index_n_eta(j)).coeffHv(2), classes{j}(index_n_eta(j)).coeffPsat(1), classes{j}(index_n_eta(j)).coeffPsat(2), classes{j}(index_n_eta(j)).coeffPsat(3), classes{j}(index_n_eta(j)).coeffPsat(4), classes{j}(index_n_eta(j)).coeffPsat(5), classes{j}(index_n_eta(j)).coeffSigma(1), classes{j}(index_n_eta(j)).coeffSigma(2), pressure);
        end

        psi_array = psi_array/sum(psi_array);

        phi = 0.0;

        for j=1:numel(psi_array)
            kappa_j = ThermophysicalProperties_SingleLiquid(property, T(i), classes{j}(index_n_eta(j)).molWeight, classes{j}(index_n_eta(j)).Tc, classes{j}(index_n_eta(j)).coeffRho(1), classes{j}(index_n_eta(j)).coeffRho(2), classes{j}(index_n_eta(j)).coeffRho(3), classes{j}(index_n_eta(j)).coeffMu(1), classes{j}(index_n_eta(j)).coeffMu(2), classes{j}(index_n_eta(j)).coeffMu(3), classes{j}(index_n_eta(j)).coeffMu(4), classes{j}(index_n_eta(j)).coeffK(1), classes{j}(index_n_eta(j)).coeffK(2), classes{j}(index_n_eta(j)).coeffK(3), classes{j}(index_n_eta(j)).coeffCl(1), classes{j}(index_n_eta(j)).coeffCl(2), classes{j}(index_n_eta(j)).coeffCl(3), classes{j}(index_n_eta(j)).coeffCl(4), classes{j}(index_n_eta(j)).coeffHv(1), classes{j}(index_n_eta(j)).coeffHv(2), classes{j}(index_n_eta(j)).coeffPsat(1), classes{j}(index_n_eta(j)).coeffPsat(2), classes{j}(index_n_eta(j)).coeffPsat(3), classes{j}(index_n_eta(j)).coeffPsat(4), classes{j}(index_n_eta(j)).coeffPsat(5), classes{j}(index_n_eta(j)).coeffSigma(1), classes{j}(index_n_eta(j)).coeffSigma(2), pressure);
            for k=1:numel(psi_array)
                kappa_k = ThermophysicalProperties_SingleLiquid(property, T(i), classes{k}(index_n_eta(k)).molWeight, classes{k}(index_n_eta(k)).Tc, classes{k}(index_n_eta(k)).coeffRho(1), classes{k}(index_n_eta(k)).coeffRho(2), classes{k}(index_n_eta(k)).coeffRho(3), classes{k}(index_n_eta(k)).coeffMu(1), classes{k}(index_n_eta(k)).coeffMu(2), classes{k}(index_n_eta(k)).coeffMu(3), classes{k}(index_n_eta(k)).coeffMu(4), classes{k}(index_n_eta(k)).coeffK(1), classes{k}(index_n_eta(k)).coeffK(2), classes{k}(index_n_eta(k)).coeffK(3), classes{k}(index_n_eta(k)).coeffCl(1), classes{k}(index_n_eta(k)).coeffCl(2), classes{k}(index_n_eta(k)).coeffCl(3), classes{k}(index_n_eta(k)).coeffCl(4), classes{k}(index_n_eta(k)).coeffHv(1), classes{k}(index_n_eta(k)).coeffHv(2), classes{k}(index_n_eta(k)).coeffPsat(1), classes{k}(index_n_eta(k)).coeffPsat(2), classes{k}(index_n_eta(k)).coeffPsat(3), classes{k}(index_n_eta(k)).coeffPsat(4), classes{k}(index_n_eta(k)).coeffPsat(5), classes{k}(index_n_eta(k)).coeffSigma(1), classes{k}(index_n_eta(k)).coeffSigma(2), pressure);
                Kij = 2.0/(1.0/kappa_j + 1.0/kappa_k);
                phi = phi + psi_array(j)*psi_array(k)*Kij;
            end
        end
    end

elseif strcmp(property,'specificHeat')

    numerator = zeros(size(T,1),numel(classes));
    denominator = zeros(size(T,1),numel(classes),1);
    
    for i = 1:numel(classes)
        numerator(:,i) = molFrac(i)*classes{i}(index_n_eta(i)).molWeight*ThermophysicalProperties_SingleLiquid(property, T, classes{i}(index_n_eta(i)).molWeight, classes{i}(index_n_eta(i)).Tc, classes{i}(index_n_eta(i)).coeffRho(1), classes{i}(index_n_eta(i)).coeffRho(2), classes{i}(index_n_eta(i)).coeffRho(3), classes{i}(index_n_eta(i)).coeffMu(1), classes{i}(index_n_eta(i)).coeffMu(2), classes{i}(index_n_eta(i)).coeffMu(3), classes{i}(index_n_eta(i)).coeffMu(4), classes{i}(index_n_eta(i)).coeffK(1), classes{i}(index_n_eta(i)).coeffK(2), classes{i}(index_n_eta(i)).coeffK(3), classes{i}(index_n_eta(i)).coeffCl(1), classes{i}(index_n_eta(i)).coeffCl(2), classes{i}(index_n_eta(i)).coeffCl(3), classes{i}(index_n_eta(i)).coeffCl(4), classes{i}(index_n_eta(i)).coeffHv(1), classes{i}(index_n_eta(i)).coeffHv(2), classes{i}(index_n_eta(i)).coeffPsat(1), classes{i}(index_n_eta(i)).coeffPsat(2), classes{i}(index_n_eta(i)).coeffPsat(3), classes{i}(index_n_eta(i)).coeffPsat(4), classes{i}(index_n_eta(i)).coeffPsat(5), classes{i}(index_n_eta(i)).coeffSigma(1), classes{i}(index_n_eta(i)).coeffSigma(2), pressure);
        denominator(:,i) = molFrac(i)*classes{i}(index_n_eta(i)).molWeight;
    end
    
    phi = sum(numerator, 2)./sum(denominator, 2);

elseif strcmp(property,'latentHeat')

    numerator = zeros(size(T,1),numel(classes));
    denominator = zeros(size(T,1),numel(classes),1);
    
    for i = 1:numel(classes)
        numerator(:,i) = molFrac(i)*classes{i}(index_n_eta(i)).molWeight*ThermophysicalProperties_SingleLiquid(property, T, classes{i}(index_n_eta(i)).molWeight, classes{i}(index_n_eta(i)).Tc, classes{i}(index_n_eta(i)).coeffRho(1), classes{i}(index_n_eta(i)).coeffRho(2), classes{i}(index_n_eta(i)).coeffRho(3), classes{i}(index_n_eta(i)).coeffMu(1), classes{i}(index_n_eta(i)).coeffMu(2), classes{i}(index_n_eta(i)).coeffMu(3), classes{i}(index_n_eta(i)).coeffMu(4), classes{i}(index_n_eta(i)).coeffK(1), classes{i}(index_n_eta(i)).coeffK(2), classes{i}(index_n_eta(i)).coeffK(3), classes{i}(index_n_eta(i)).coeffCl(1), classes{i}(index_n_eta(i)).coeffCl(2), classes{i}(index_n_eta(i)).coeffCl(3), classes{i}(index_n_eta(i)).coeffCl(4), classes{i}(index_n_eta(i)).coeffHv(1), classes{i}(index_n_eta(i)).coeffHv(2), classes{i}(index_n_eta(i)).coeffPsat(1), classes{i}(index_n_eta(i)).coeffPsat(2), classes{i}(index_n_eta(i)).coeffPsat(3), classes{i}(index_n_eta(i)).coeffPsat(4), classes{i}(index_n_eta(i)).coeffPsat(5), classes{i}(index_n_eta(i)).coeffSigma(1), classes{i}(index_n_eta(i)).coeffSigma(2), pressure);
        denominator(:,i) = molFrac(i)*classes{i}(index_n_eta(i)).molWeight;
    end
    
    phi = sum(numerator, 2)./sum(denominator, 2);

elseif strcmp(property,'vaporPressure')

    numerator = zeros(size(T,1),numel(classes));
    denominator = zeros(size(T,1),numel(classes),1);
    
    for i = 1:numel(classes)
        numerator(:,i) = molFrac(i)*classes{i}(index_n_eta(i)).molWeight*ThermophysicalProperties_SingleLiquid(property, T, classes{i}(index_n_eta(i)).molWeight, classes{i}(index_n_eta(i)).Tc, classes{i}(index_n_eta(i)).coeffRho(1), classes{i}(index_n_eta(i)).coeffRho(2), classes{i}(index_n_eta(i)).coeffRho(3), classes{i}(index_n_eta(i)).coeffMu(1), classes{i}(index_n_eta(i)).coeffMu(2), classes{i}(index_n_eta(i)).coeffMu(3), classes{i}(index_n_eta(i)).coeffMu(4), classes{i}(index_n_eta(i)).coeffK(1), classes{i}(index_n_eta(i)).coeffK(2), classes{i}(index_n_eta(i)).coeffK(3), classes{i}(index_n_eta(i)).coeffCl(1), classes{i}(index_n_eta(i)).coeffCl(2), classes{i}(index_n_eta(i)).coeffCl(3), classes{i}(index_n_eta(i)).coeffCl(4), classes{i}(index_n_eta(i)).coeffHv(1), classes{i}(index_n_eta(i)).coeffHv(2), classes{i}(index_n_eta(i)).coeffPsat(1), classes{i}(index_n_eta(i)).coeffPsat(2), classes{i}(index_n_eta(i)).coeffPsat(3), classes{i}(index_n_eta(i)).coeffPsat(4), classes{i}(index_n_eta(i)).coeffPsat(5), classes{i}(index_n_eta(i)).coeffSigma(1), classes{i}(index_n_eta(i)).coeffSigma(2), pressure);
        denominator(:,i) = molFrac(i)*classes{i}(index_n_eta(i)).molWeight;
    end
    
    phi = sum(numerator, 2)./sum(denominator, 2);

elseif strcmp(property,'sigma')

    sigma = zeros(size(T,1),numel(classes));
    Xs = zeros(size(T,1),numel(classes)); % surface molar fractions estimated by Raoult's law
    
    for i = 1:numel(classes)
        sigma(:,i) = ThermophysicalProperties_SingleLiquid(property, T, classes{i}(index_n_eta(i)).molWeight, classes{i}(index_n_eta(i)).Tc, classes{i}(index_n_eta(i)).coeffRho(1), classes{i}(index_n_eta(i)).coeffRho(2), classes{i}(index_n_eta(i)).coeffRho(3), classes{i}(index_n_eta(i)).coeffMu(1), classes{i}(index_n_eta(i)).coeffMu(2), classes{i}(index_n_eta(i)).coeffMu(3), classes{i}(index_n_eta(i)).coeffMu(4), classes{i}(index_n_eta(i)).coeffK(1), classes{i}(index_n_eta(i)).coeffK(2), classes{i}(index_n_eta(i)).coeffK(3), classes{i}(index_n_eta(i)).coeffCl(1), classes{i}(index_n_eta(i)).coeffCl(2), classes{i}(index_n_eta(i)).coeffCl(3), classes{i}(index_n_eta(i)).coeffCl(4), classes{i}(index_n_eta(i)).coeffHv(1), classes{i}(index_n_eta(i)).coeffHv(2), classes{i}(index_n_eta(i)).coeffPsat(1), classes{i}(index_n_eta(i)).coeffPsat(2), classes{i}(index_n_eta(i)).coeffPsat(3), classes{i}(index_n_eta(i)).coeffPsat(4), classes{i}(index_n_eta(i)).coeffPsat(5), classes{i}(index_n_eta(i)).coeffSigma(1), classes{i}(index_n_eta(i)).coeffSigma(2), pressure);
        pSat = ThermophysicalProperties_SingleLiquid('vaporPressure', T, classes{i}(index_n_eta(i)).molWeight, classes{i}(index_n_eta(i)).Tc, classes{i}(index_n_eta(i)).coeffRho(1), classes{i}(index_n_eta(i)).coeffRho(2), classes{i}(index_n_eta(i)).coeffRho(3), classes{i}(index_n_eta(i)).coeffMu(1), classes{i}(index_n_eta(i)).coeffMu(2), classes{i}(index_n_eta(i)).coeffMu(3), classes{i}(index_n_eta(i)).coeffMu(4), classes{i}(index_n_eta(i)).coeffK(1), classes{i}(index_n_eta(i)).coeffK(2), classes{i}(index_n_eta(i)).coeffK(3), classes{i}(index_n_eta(i)).coeffCl(1), classes{i}(index_n_eta(i)).coeffCl(2), classes{i}(index_n_eta(i)).coeffCl(3), classes{i}(index_n_eta(i)).coeffCl(4), classes{i}(index_n_eta(i)).coeffHv(1), classes{i}(index_n_eta(i)).coeffHv(2), classes{i}(index_n_eta(i)).coeffPsat(1), classes{i}(index_n_eta(i)).coeffPsat(2), classes{i}(index_n_eta(i)).coeffPsat(3), classes{i}(index_n_eta(i)).coeffPsat(4), classes{i}(index_n_eta(i)).coeffPsat(5), classes{i}(index_n_eta(i)).coeffSigma(1), classes{i}(index_n_eta(i)).coeffSigma(2), pressure);
        Xs(:,i) = molFrac(i).*pSat./pressure;
    end

    Xs = Xs./sum(Xs,2);
    
    phi = sum(Xs.*sigma,2);

end

end