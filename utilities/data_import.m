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

function  [fullData, minT_array, maxT_array, variable_names] = data_import(Dataset, families, n_ranges, Tbubble_realFuel, numComponents, pressure_distillation)

% DESCRIPTION:
% The data_import function provides tools to import information about
% the experimental measurements concerning the real fuel.

% Auxiliary parameters:
% Np: number of thermophysical properties targeted during the formulation
% of the surrogate mixture  [-]

% Inputs:
% 1) Dataset                 : (1 x Np) cell array, with the i-th cell containing a
% string denoting the name of the experimental dataset about the i-th thermophysical property
% 2) families                : (1 x numComponents) cell array, with the i-th cell containing a string
% that denotes the i-th hydrocarbon family under consideration
% 3) n_ranges                : (1 x numComponents) cell array, with the i-th cell containing the range
% the number of carbon atoms of the i-th surrogate mixture component can vary within
% 4) Tbubble_realFuel        : bubble temperature of the real fuel  [K]
% 5) numComponents           : number of surrogate mixture components  [-]
% 6) pressure_distillation   : pressure the distillation curve is computed at  [Pa]

% Outputs:
% 1) fullData      : (1 x Np) cell array, with the i-th cell being a (Nd x 3) array
% that contains the experimental measurements about the thermophysical property under
% consideration (independent variable, dependent thermophysical property, standard deviation). 
% Note that the i-th cell is a (Nd x 4) array in case the distillation curve is the i-th
% thermophysical property under consideration.
% 2) minT_array    : (1 x Np) array, with the i-th element denoting the minimum
% temperature in the dataset of the i-th thermophysical property
% 3) maxT_array    : (1 x Np) array, with the i-th element denoting the maximum
% temperature in the dataset of the i-th thermophysical property
% 4) variable_names: (1 x Np) cell array, with the i-th cell containing a string
% that denotes the i-th thermophysical property under consideration

fullData = cell(1,numel(Dataset));
minT_array = zeros(1,numel(Dataset));
maxT_array = zeros(1,numel(Dataset));
variable_names = cell(1,numel(Dataset));

for i = 1:numel(Dataset)

    if contains(Dataset{i},'rho')
        variable_name = 'rho';
    end

    if contains(Dataset{i},'mu')
        variable_name = 'mu';
    end

    if contains(Dataset{i},'nu')
        variable_name = 'nu';
    end

    if contains(Dataset{i},'kappa')
        variable_name = 'kappa';
    end

    if contains(Dataset{i},'specificHeat')
        variable_name = 'specificHeat';
    end

    if contains(Dataset{i},'latentHeat')
        variable_name = 'latentHeat';
    end

    if contains(Dataset{i},'vaporPressure')
        variable_name = 'vaporPressure';
    end

    if contains(Dataset{i},'sigma')
        variable_name = 'sigma';
    end

    if contains(Dataset{i},'distillation')
        variable_name = 'distillation';
    end

    if contains(Dataset{i},'molWeight')
        variable_name = 'molWeight';
    end

    if contains(Dataset{i},'HC')
        variable_name = 'HC';
    end

    if contains(Dataset{i},'DCN')
        variable_name = 'DCN';
    end

    if contains(Dataset{i},'flash')
        variable_name = 'flash';
    end

    if contains(Dataset{i},'freezing')
        variable_name = 'freezing';
    end

    if contains(Dataset{i},'LHV')
        variable_name = 'LHV';
    end

    variable_names{i} = variable_name;

    % Define a threshold temperature as the minimum between the minimum boiling
    % temperature of candidate surrogate components and the bubble temperature
    % of the real fuel. Every experimental measurement above this threshold
    % temperature will be discarded during Bayesian inference
    T_threshold = threshold_temperature(families, n_ranges, Tbubble_realFuel, numComponents, pressure_distillation);

    % Import experimental data and calculate the standard deviation associated with experimental results
    data_import = table2array(readtable(Dataset{i}, 'VariableNamingRule', 'preserve'));

    if contains(Dataset{i},'distillation') || contains(Dataset{i},'molWeight') || contains(Dataset{i},'HC') || contains(Dataset{i},'DCN') || contains(Dataset{i},'flash') || contains(Dataset{i},'freezing') || contains(Dataset{i},'LHV')
        data = data_import;
    else
        data = data_import(0.999*T_threshold > data_import(:, 1), :);
    end

    if ~contains(Dataset{i},'distillation')
        minT_array(i) = min(data(:,1));
        maxT_array(i) = max(data(:,1));
    else
        minT_array(i) = min(data(:,2));
        maxT_array(i) = max(data(:,2));
    end

    fullData{i} = data;

end

end