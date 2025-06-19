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

function [volumeFraction, T_distillation]  = DistillationCurve(X0, classes, index_n_eta, pressure)

% DESCRIPTION:
% The DistillationCurve function computes the distillation curve of the given surrogate mixture at a fixed pressure.

% Auxiliary parameters:
% Nc: number of surrogate mixture components  [-]

% Inputs:
% 1) X0         : (1 x Nc-1) array containing the molar fractions of the surrogate
% components in the liquid phase before the distillation process begins  [-]
% 2) classes    : (1 x Nc) cell array, with the i-th cell containing a structure
% array related to the i-th hydrocarbon family, range of number of carbon
% atoms, and range of topochemical atom index
% 3) index_n_eta: (1 x Nc) array, with the i-th element representing the
% index of the species in the i-th cell of "classes" displaying the
% closest topochemical atom index to the index currently investigated for
% the i-th surrogate mixture component
% 4) pressure   : pressure the distillation curve is computed at  [Pa]

% Outputs:
% 1) volumeFraction: (1 x 101) array containing the values of the recovered volume fraction  [-]
% 2) T_distillation: (1 x 101) array containing the temperature values  [K]

mol=0:1:100;

X = X0;

T_distillation = 0:1:100;
volumeFraction = 0:1:100;
Y = zeros(numel(classes)-1,1);

options = optimset('Display', 'off', 'TolFun', 1e-6);

beta = 1e-6;

T0 = 300;

VLE = @(T) Flash(T, beta, classes, index_n_eta, X, pressure);

[temp, ~] = fsolve(VLE,T0,options);

T_distillation(1) = temp;

for i=1:length(mol)-1

    beta = 1/(100.-mol(i));

    VLE = @(T) Flash(T, beta, classes, index_n_eta, X, pressure);

    [temp, ~] = fsolve(VLE,temp,options);

    T_distillation(i+1) = temp;

    if i==1
        rho_i = zeros(numel(classes),1);
        vol_i = zeros(numel(classes),1);
        for j = 1:numel(classes)
            rho_i(j) = ThermophysicalProperties_SingleLiquid('rho', temp, classes{j}(index_n_eta(j)).molWeight, classes{j}(index_n_eta(j)).Tc, classes{j}(index_n_eta(j)).coeffRho(1), classes{j}(index_n_eta(j)).coeffRho(2), classes{j}(index_n_eta(j)).coeffRho(3), classes{j}(index_n_eta(j)).coeffMu(1), classes{j}(index_n_eta(j)).coeffMu(2), classes{j}(index_n_eta(j)).coeffMu(3), classes{j}(index_n_eta(j)).coeffMu(4), classes{j}(index_n_eta(j)).coeffK(1), classes{j}(index_n_eta(j)).coeffK(2), classes{j}(index_n_eta(j)).coeffK(3), classes{j}(index_n_eta(j)).coeffCl(1), classes{j}(index_n_eta(j)).coeffCl(2), classes{j}(index_n_eta(j)).coeffCl(3), classes{j}(index_n_eta(j)).coeffCl(4), classes{j}(index_n_eta(j)).coeffHv(1), classes{j}(index_n_eta(j)).coeffHv(2), classes{j}(index_n_eta(j)).coeffPsat(1), classes{j}(index_n_eta(j)).coeffPsat(2), classes{j}(index_n_eta(j)).coeffPsat(3), classes{j}(index_n_eta(j)).coeffPsat(4), classes{j}(index_n_eta(j)).coeffPsat(5), classes{j}(index_n_eta(j)).coeffSigma(1), classes{j}(index_n_eta(j)).coeffSigma(2), pressure);
            if j < numel(classes)
                vol_i(j) = X0(j)*classes{j}(index_n_eta(j)).molWeight/rho_i(j);
            else
                vol_i(j) = (1-sum(X0))*classes{j}(index_n_eta(j)).molWeight/rho_i(j);
            end
        end
        V0 = sum(vol_i);
    end

    pSat = zeros(numel(classes),1);
    K = zeros(numel(classes),1);
    Vtilde_i = zeros(numel(classes),1);
    for k = 1:numel(classes)
        pSat(k) = ThermophysicalProperties_SingleLiquid('vaporPressure', temp, classes{k}(index_n_eta(k)).molWeight, classes{k}(index_n_eta(k)).Tc, classes{k}(index_n_eta(k)).coeffRho(1), classes{k}(index_n_eta(k)).coeffRho(2), classes{k}(index_n_eta(k)).coeffRho(3), classes{k}(index_n_eta(k)).coeffMu(1), classes{k}(index_n_eta(k)).coeffMu(2), classes{k}(index_n_eta(k)).coeffMu(3), classes{k}(index_n_eta(k)).coeffMu(4), classes{k}(index_n_eta(k)).coeffK(1), classes{k}(index_n_eta(k)).coeffK(2), classes{k}(index_n_eta(k)).coeffK(3), classes{k}(index_n_eta(k)).coeffCl(1), classes{k}(index_n_eta(k)).coeffCl(2), classes{k}(index_n_eta(k)).coeffCl(3), classes{k}(index_n_eta(k)).coeffCl(4), classes{k}(index_n_eta(k)).coeffHv(1), classes{k}(index_n_eta(k)).coeffHv(2), classes{k}(index_n_eta(k)).coeffPsat(1), classes{k}(index_n_eta(k)).coeffPsat(2), classes{k}(index_n_eta(k)).coeffPsat(3), classes{k}(index_n_eta(k)).coeffPsat(4), classes{k}(index_n_eta(k)).coeffPsat(5), classes{k}(index_n_eta(k)).coeffSigma(1), classes{k}(index_n_eta(k)).coeffSigma(2), pressure);
        K(k) = pSat(k)/pressure;
        if k < numel(classes)
            X(k) = X(k)/(1+beta*(K(k)-1));
            Y(k) = K(k)*X(k);
            Vtilde_i(k) = Y(k)*classes{k}(index_n_eta(k)).molWeight/rho_i(k);
        else
            Vtilde_i(k) = (1-sum(Y))*classes{k}(index_n_eta(k)).molWeight/rho_i(k);
        end
    end

    Vtilde = sum(Vtilde_i);

    volumeFraction(i+1) = volumeFraction(i) + Vtilde/V0;

end

end