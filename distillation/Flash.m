function F = Flash(T, beta, classes, index_n_eta, X, pressure)

% DESCRIPTION:
% The Flash function returns the left-hand side of a flash distillation problem for the given
% surrogate mixture at a fixed temperature, recovered mole fraction,
% liquid-phase molar composition, and pressure.

% Auxiliary parameters:
% Nc: number of surrogate mixture components  [-]

% Inputs:
% 1) T          : temperature  [K]
% 2) beta       : recovered mole fraction [-]
% 3) classes    : (1 x Nc) cell array, with the i-th cell containing a structure
% array related to the i-th hydrocarbon family, range of number of carbon
% atoms, and range of topochemical atom index
% 4) index_n_eta: (1 x Nc) array, with the i-th element representing the
% index of the species in the i-th cell of "classes" displaying the
% closest topochemical atom index to the index currently investigated for
% the i-th surrogate mixture component
% 5) X0         : (1 x Nc-1) array containing the molar fractions of the surrogate
% components in the liquid phase
% 6) pressure   : pressure the distillation curve is computed at  [Pa]

% Output:
% 1) F: left-hand side of the flash distillation problem

% ------------------------------------------------------------------------%
% Contributors/Copyright
% 2024 Jacopo Liberatori, jacopo.liberatori@uniroma1.it
% 2024 Davide Cavalieri,  davide.cavalieri@uniroma1.it
% Department of Mechanical and Aerospace Engineering (DIMA)
% Sapienza University of Rome
% ------------------------------------------------------------------------%

pSat = zeros(numel(classes),1);
flash = zeros(numel(classes),1);
for i = 1:numel(classes)
    pSat(i) = ThermophysicalProperties_SingleLiquid('vaporPressure', T, classes{i}(index_n_eta(i)).molWeight, classes{i}(index_n_eta(i)).Tc, classes{i}(index_n_eta(i)).coeffRho(1), classes{i}(index_n_eta(i)).coeffRho(2), classes{i}(index_n_eta(i)).coeffRho(3), classes{i}(index_n_eta(i)).coeffMu(1), classes{i}(index_n_eta(i)).coeffMu(2), classes{i}(index_n_eta(i)).coeffMu(3), classes{i}(index_n_eta(i)).coeffMu(4), classes{i}(index_n_eta(i)).coeffK(1), classes{i}(index_n_eta(i)).coeffK(2), classes{i}(index_n_eta(i)).coeffK(3), classes{i}(index_n_eta(i)).coeffCl(1), classes{i}(index_n_eta(i)).coeffCl(2), classes{i}(index_n_eta(i)).coeffCl(3), classes{i}(index_n_eta(i)).coeffCl(4), classes{i}(index_n_eta(i)).coeffHv(1), classes{i}(index_n_eta(i)).coeffHv(2), classes{i}(index_n_eta(i)).coeffPsat(1), classes{i}(index_n_eta(i)).coeffPsat(2), classes{i}(index_n_eta(i)).coeffPsat(3), classes{i}(index_n_eta(i)).coeffPsat(4), classes{i}(index_n_eta(i)).coeffPsat(5), classes{i}(index_n_eta(i)).coeffSigma(1), classes{i}(index_n_eta(i)).coeffSigma(2), pressure);
    if i < numel(classes)
        flash(i) = X(i)*(pSat(i)/pressure - 1)/(1+beta*(pSat(i)/pressure - 1));
    else
        flash(i) = (1-sum(X))*(pSat(i)/pressure - 1)/(1+beta*(pSat(i)/pressure - 1));
    end
end

F = sum(flash);

end