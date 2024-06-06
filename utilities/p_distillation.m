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

% ------------------------------------------------------------------------%
% Contributors/Copyright
% 2024 Jacopo Liberatori, jacopo.liberatori@uniroma1.it
% 2024 Davide Cavalieri,  davide.cavalieri@uniroma1.it
% Department of Mechanical and Aerospace Engineering (DIMA)
% Sapienza University of Rome
% ------------------------------------------------------------------------%

if any(contains(Dataset,'distillation'))

    index = find(contains(Dataset, 'distillation'));
    data_import = table2array(readtable(Dataset{index}, 'VariableNamingRule', 'preserve'));
    pressure_distillation = data_import(1, 4);

else

    pressure_distillation = 101325;

end

end