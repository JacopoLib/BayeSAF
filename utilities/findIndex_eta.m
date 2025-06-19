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

function  row_idx = findIndex_eta(class, nC, eta_B_star)

% DESCRIPTION:
% The findIndex_eta function returns the row index in the given hydrocarbon family structure array
% corresponding to the chemical species under consideration.

% Auxiliary parameters:
% Np: number of thermophysical properties targeted during the formulation
% of the surrogate mixture  [-]

% Inputs:
% 1) class  : structure array for the hydrocarbon family the chemical species under
% consideration belongs to, containing detailed information about the whole
% set of chemical species (number of carbon atoms, normalized topochemical
% atom index, molecular weight [g/mol], critical temperature [K], set of
% coefficients for the Yaws' polynomials concerning liquid-phase dynamic viscosity,
% liquid-phase density, vapor pressure, liquid-phase thermal conductivity, 
% liquid-phase specific heat capacity, heat of vaporization and liquid-phase surface tension)
% 2) nC        : number of carbon atoms for the chemical species under consideration
% 3) eta_B_star: topochemical atom index for the chemical species under consideration

% Outputs:
% 1) row_idx: row index in the given hydrocarbon family structure array
% corresponding to the chemical species under consideration

% Find elements in the structure array displaying a number of carbon atoms equal to nC
filtered_indices = find([class.nC] == nC);
% Extract the eta_B_star_norm values for these elements
eta_B_star_values = [class.eta_B_star_norm];
filtered_eta_values = eta_B_star_values(filtered_indices);

% Compute the absolute difference between these eta values and the target eta
differences = abs(filtered_eta_values - eta_B_star);

% Find the index of the minimum difference
[~, min_index_relative] = min(differences);

% Get the original index in the data structure array
row_idx = filtered_indices(min_index_relative);

end