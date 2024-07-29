function  Wavg = molWeight(Xi, Wi)

% DESCRIPTION:
% The mol2mass function calculates the mean molecular weight of the surrogate mixture
% given the molar fractions and molecular weights of its components.

% Auxiliary parameters:
% Nc: number of surrogate mixture components  [-]

% Inputs:
% 1) Xi: (1 x Nc-1) array containing the molar fractions of the surrogate mixture components
% 2) Wi: (1 x Nc) array containing the molecular weights of the surrogate mixture components

% Outputs:
% 1) Wavg: mean molecular weight of the surrogate mixture [g/mol]

% ------------------------------------------------------------------------%
% Contributors/Copyright
% 2024 Jacopo Liberatori, jacopo.liberatori@uniroma1.it
% 2024 Davide Cavalieri,  davide.cavalieri@uniroma1.it
% Department of Mechanical and Aerospace Engineering (DIMA)
% Sapienza University of Rome
% ------------------------------------------------------------------------%

Wavg  = sum(Xi.*Wi);

end