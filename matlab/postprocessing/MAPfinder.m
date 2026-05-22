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
% The MAPfinder function returns the molar fractions, numbers of carbon atoms,
% and topochemical atom indices characterizing the maximum a posteriori (MAP) surrogate.

% Auxiliary parameters:
% Nc: number of surrogate mixture components  [-]
% numIterations: number of iterations for each chain in the DE-MC algorithm  [-]
% t_convergence: number of iterations the convergence of the DE-MC run is
% reached at, checking the R-hat statistic  [-]
% t_burnin: number of iterations in the DE-MC run to be discarded as burn-in [-]
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
% ------------------------------------------------------------------------

function [x_MAP, nc_MAP, eta_B_star_MAP] = MAPfinder(posterior, chain, numComponents)

[~, max_posterior_index] = maxk(posterior(:,1), 1);
x_MAP = chain(max_posterior_index,1:numComponents-1);
x_MAP = [x_MAP, 1-sum(x_MAP)];
nc_MAP = chain(max_posterior_index,numComponents:2*numComponents-1);

eta_B_star_MAP = chain(max_posterior_index,2*numComponents:3*numComponents-1);

end
