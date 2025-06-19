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

function  T_threshold = threshold_temperature(families, n_ranges, Tbubble_realFuel, numComponents, pressure_distillation)

% DESCRIPTION:
% The threshold_temperature function calculates the threshold temperature above
% which experimental measurements about thermophysical properties of the real fuel
% (except for the distillation curve data) are discarded.

% Inputs:
% 1) families             : (1 x numComponents) cell array, with the i-th cell containing a string
% that denotes the i-th hydrocarbon family under consideration
% 2) n_ranges             : (1 x numComponents) cell array, with the i-th cell containing the range
% the number of carbon atoms of the i-th surrogate mixture component can vary within
% 3) Tbubble_realFuel     : bubble temperature of the real fuel  [K]
% 4) numComponents        : number of surrogate mixture components  [-]
% 5) pressure_distillation: pressure the distillation curve is computed at  [Pa]

% Outputs:
% 1) T_threshold: threshold temperature calculated as the minimum between Tbubble_realFuel
% and the minimum boiling temperature (calculated at pressure_distillation) among the 
% surrogate mixture candidate components  [K]

% Initialize minimum boiling temperature among all candidate surrogate components
TbMin = 1e+18;

% The minimum boiling temperature among all candidate surrogate components is calculated
% Initialize empty cell array to store the classes
for j = 1:numComponents
    class = Hydrocarbons(families{j}, n_ranges{j});
    for k = 1:numel(n_ranges{j})
        TbMin = min(TbMin, ThermophysicalProperties_SingleLiquid('boilingTemperature', 0, class(k).molWeight, class(k).Tc, class(k).coeffRho(1), class(k).coeffRho(2), class(k).coeffRho(3), class(k).coeffMu(1), class(k).coeffMu(2), class(k).coeffMu(3), class(k).coeffMu(4), class(k).coeffK(1), class(k).coeffK(2), class(k).coeffK(3), class(k).coeffCl(1), class(k).coeffCl(2), class(k).coeffCl(3), class(k).coeffCl(4), class(k).coeffHv(1), class(k).coeffHv(2), class(k).coeffPsat(1), class(k).coeffPsat(2), class(k).coeffPsat(3), class(k).coeffPsat(4), class(k).coeffPsat(5), class(k).coeffSigma(1), class(k).coeffSigma(2), pressure_distillation));
    end
end

T_threshold = min(TbMin, Tbubble_realFuel);

end