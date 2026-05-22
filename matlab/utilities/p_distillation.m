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
% ------------------------------------------------------------------------

function  pressure_distillation = p_distillation(Dataset)

if any(contains(Dataset,'distillation'))

    index = find(contains(Dataset, 'distillation'));
    data_import = table2array(readtable(Dataset{index}, 'VariableNamingRule', 'preserve'));
    pressure_distillation = data_import(1, 4);

elseif any(contains(Dataset,'deltaT_dist'))

    index = find(contains(Dataset, 'deltaT_dist'));
    data_import = table2array(readtable(Dataset{index}, 'VariableNamingRule', 'preserve'));
    pressure_distillation = data_import(1, 4);

else

    pressure_distillation = 101325;

end

end