function  row_idx = findIndex_eta(class, eta_B_star)

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
% liquid-phase specific heat capacity, heat of vaporization and liquid-phase surface tension).
% 2) eta_B_star: topochemical atom index for the chemical species
% under consideration

% Outputs:
% 1) row_idx: row index in the given hydrocarbon family structure array
% corresponding to the chemical species under consideration

% ------------------------------------------------------------------------%
% Contributors/Copyright
% 2024 Jacopo Liberatori, jacopo.liberatori@uniroma1.it
% 2024 Davide Cavalieri,  davide.cavalieri@uniroma1.it
% Department of Mechanical and Aerospace Engineering (DIMA)
% Sapienza University of Rome
% ------------------------------------------------------------------------%

eta_B_star_values = [class.eta_B_star_norm];
eta_B_star_diff = abs(eta_B_star_values - eta_B_star);
min_eta_B_star_diff = min(abs(eta_B_star_values - eta_B_star));
row_idx = find(eta_B_star_diff == min_eta_B_star_diff, 1);

end