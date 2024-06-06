function  Yi = mol2mass(Xi, Wi)

% DESCRIPTION:
% The mol2mass function calculates the mass fractions of the surrogate mixture components
% given their molar fractions and molecular weights.

% Auxiliary parameters:
% Nc: number of surrogate mixture components  [-]

% Inputs:
% 1) Xi: (1 x Nc-1) array containing the molar fractions of the surrogate mixture components
% 2) Wi: (1 x Nc) array containing the molecular weights of the surrogate mixture components

% Outputs:
% 1) Yi: (1 x Nc) array containing the mass fractions of the surrogate mixture components

% ------------------------------------------------------------------------%
% Contributors/Copyright
% 2024 Jacopo Liberatori, jacopo.liberatori@uniroma1.it
% 2024 Davide Cavalieri,  davide.cavalieri@uniroma1.it
% Department of Mechanical and Aerospace Engineering (DIMA)
% Sapienza University of Rome
% ------------------------------------------------------------------------%

% Mean molecular weight of the surrogate mixture [g/mol]
Wavg  = sum(Xi.*Wi);

% Mass fractions of the surrogate mixture components
Yi = (Xi.*Wi)/Wavg;

end