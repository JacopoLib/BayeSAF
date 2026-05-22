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
% dirichlet_pdf returns the logarithmic probability density of the
% Dirichlet distribution for the given molar-fraction vector. Returns
% -Inf if any component violates its bound constraints.
%
% Inputs:
% 1) x         : (1 x Nc-1) array of molar fractions  [-]
% 2) alpha     : (1 x Nc) concentration parameters  [-]
% 3) LowerBound: (1 x Nc) lower bounds for molar fractions  [-]
% 4) UpperBound: (1 x Nc) upper bounds for molar fractions  [-]
%
% Outputs:
% 1) prob: log probability density of the Dirichlet distribution  [-]
% ------------------------------------------------------------------------

function prob = dirichlet_pdf(x, alpha, LowerBound, UpperBound)

x_full = [x, 1 - sum(x)];

if any(x_full < LowerBound) || any(x_full > UpperBound)
    prob = -Inf;
    return;
end

if any(alpha <= 0)
    error('Alpha values must be strictly positive.');
end

log_const = gammaln(sum(alpha)) - sum(gammaln(alpha));
prob = log_const + sum((alpha - 1) .* log(x_full));

end
