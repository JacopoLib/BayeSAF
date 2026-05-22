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
% uniform_pdf returns the logarithmic probability density of a continuous
% uniform distribution over the topochemical atom index of each surrogate
% mixture component. Returns -Inf if any component lies outside its bounds.
%
% Inputs:
% 1) x         : (1 x Nc) array of topochemical atom indices  [-]
% 2) LowerBound: (1 x Nc) lower bounds  [-]
% 3) UpperBound: (1 x Nc) upper bounds  [-]
%
% Outputs:
% 1) prob: log probability density of the continuous uniform distribution  [-]
% ------------------------------------------------------------------------

function prob = uniform_pdf(x, LowerBound, UpperBound)

prob = 0;
for j = 1:numel(x)
    if x(j) < LowerBound(j) || x(j) > UpperBound(j)
        prob = -Inf;
        return;
    end
    prob = prob - log(UpperBound(j) - LowerBound(j));
end

end
