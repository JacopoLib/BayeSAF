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

function prob = dirichlet_pdf(x, alpha, LowerBound, UpperBound)

% DESCRIPTION:
% The dirichlet_pdf function returns the logarithmic probability density
% of the Dirichlet distribution for the set of molar fractions.

% Auxiliary parameters:
% Nc: number of surrogate mixture components  [-]

% Inputs:
% 1) x         : (1 x Nc-1) array containing the molar fractions of the surrogate components  [-]
% 2) alpha     : (1 x Nc) vector of concentration parameters for the Dirichlet distribution  [-]
% 3) LowerBound: (1 x Nc) array, with the i-th element representing the
% lower bound for the range the mole fraction of the i-th surrogate mixture
% 4) UpperBound: (1 x Nc) array, with the i-th element representing the
% upper bound for the range the mole fraction of the i-th surrogate mixture

% Outputs:
% 1) prob: logarithmic probability density of the Dirichlet distribution

% Reconstruct full vector of molar fractions
x_full = [x 1-sum(x)];

% Check that the inputs are valid
if any(x_full < LowerBound) || any(x_full > UpperBound)
    prob = -Inf;
    return;
end

% Compute the logarithmic probability density
if any(alpha <= 0)
    error('Alpha values must be strictly positive.');
else
    log_const = gammaln(sum(alpha)) - sum(gammaln(alpha));
    prob = log_const + sum((alpha - 1) .* log(x_full));
end
end
