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
% The mol2mass function calculates the mass fractions of the surrogate mixture components
% given their molar fractions and molecular weights.

% Auxiliary parameters:
% Nc: number of surrogate mixture components  [-]

% Inputs:
% 1) Xi: (1 x Nc) array containing the molar fractions of the surrogate mixture components
% 2) Wi: (1 x Nc) array containing the molecular weights of the surrogate mixture components

% Outputs:
% 1) Yi: (1 x Nc) array containing the mass fractions of the surrogate mixture components
% ------------------------------------------------------------------------

function  Yi = mol2mass(Xi, Wi)

% Mean molecular weight of the surrogate mixture [g/mol]
Wavg  = sum(Xi.*Wi);

% Mass fractions of the surrogate mixture components
Yi = (Xi.*Wi)/Wavg;

end