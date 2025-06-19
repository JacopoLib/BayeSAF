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

function  Yi = mol2mass(Xi, Wi)

% DESCRIPTION:
% The mol2mass function calculates the mass fractions of the surrogate mixture components
% given their molar fractions and molecular weights.

% Auxiliary parameters:
% Nc: number of surrogate mixture components  [-]

% Inputs:
% 1) Xi: (1 x Nc) array containing the molar fractions of the surrogate mixture components
% 2) Wi: (1 x Nc) array containing the molecular weights of the surrogate mixture components

% Outputs:
% 1) Yi: (1 x Nc) array containing the mass fractions of the surrogate mixture components

% Mean molecular weight of the surrogate mixture [g/mol]
Wavg  = sum(Xi.*Wi);

% Mass fractions of the surrogate mixture components
Yi = (Xi.*Wi)/Wavg;

end