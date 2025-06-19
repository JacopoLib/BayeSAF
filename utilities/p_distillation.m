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

function  pressure_distillation = p_distillation(Dataset)

% DESCRIPTION:
% The p_distillation function calculates the pressure the distillation
% curve is computed at.

% Auxiliary parameters:
% Np: number of thermophysical properties targeted during the formulation
% of the surrogate mixture  [-]

% Inputs:
% 1) Dataset: (1 x Np) cell array, with the i-th cell containing a
% string denoting the name of the experimental dataset about the i-th thermophysical property

% Outputs:
% 1) pressure_distillation: pressure the distillation curve is computed at  [Pa]

if any(contains(Dataset,'distillation'))

    index = find(contains(Dataset, 'distillation'));
    data_import = table2array(readtable(Dataset{index}, 'VariableNamingRule', 'preserve'));
    pressure_distillation = data_import(1, 4);

else

    pressure_distillation = 101325;

end

end