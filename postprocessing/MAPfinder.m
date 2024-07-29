function [x_MAP, nc_MAP, eta_B_star_MAP] = MAPfinder(posterior, chain, numComponents)

% DESCRIPTION:
% The MAPfinder function returns the molar fractions, numbers of carbon atoms,
% and topochemical atom indices characterizing the maximum a posteriori (MAP) surrogate.

% Auxiliary parameters:
% Nc: number of surrogate mixture components  [-]
% numIterations: number of iterations for each chain in the DE-MC algorithm  [-]
% t_convergence: number of iterations the convergence of the DE-MC run is
% reached at, checking the R-hat statistic  [-]
% N_chains: number of chains in the DE-MC algorithm  [-]

% Inputs:
% 1) posterior    : (numIterations-t_convergence x N_chains) array containing the posterior density
% values from each chain in the DE-MC algorithm discarding the burn-in period
% 2) chain    : (numIterations-t_convergence x 3*Nc x N_chains) array containing the
% samples from the posterior PDF from each chain in the DE-MC algorithm discarding the burn-in period,
% concerning the molar fractions, numbers of carbon atoms, and topochemical atom indices
% 3) numComponents: number of surrogate mixture components  [-]


% Outputs:
% 1) x_MAP: molar fractions characterizing the MAP surrogate components
% 2) nc_MAP: numbers of carbon atoms characterizing the MAP surrogate components
% 3) eta_B_star_MAP: topochemical atom indices characterizing the MAP surrogate components

% ------------------------------------------------------------------------%
% Contributors/Copyright
% 2024 Jacopo Liberatori, jacopo.liberatori@uniroma1.it
% 2024 Davide Cavalieri,  davide.cavalieri@uniroma1.it
% Department of Mechanical and Aerospace Engineering (DIMA)
% Sapienza University of Rome
% ------------------------------------------------------------------------%

[~, max_posterior_index] = maxk(posterior(:,1), 1);
x_MAP = chain(max_posterior_index,1:numComponents-1);
x_MAP = [x_MAP, 1-sum(x_MAP)];
nc_MAP = chain(max_posterior_index,numComponents:2*numComponents-1);
eta_B_star_MAP = chain(max_posterior_index,2*numComponents:3*numComponents-1);

end