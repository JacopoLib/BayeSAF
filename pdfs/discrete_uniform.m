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

function prob = discrete_uniform(n, n_ranges)

% DESCRIPTION:
% The discrete_uniform function calculates the logarithmic probability density of a discrete
% uniform distribution of the number of carbon atoms of each surrogate mixture component.

% Auxiliary parameters:
% Nc: number of surrogate mixture components  [-]

% Inputs:
% 1) n       : (1 x Nc) array containing the numbers of carbon atoms of the
% surrogate components  [-]
% 2) n_ranges: (1 x Nc) cell array, with the i-th cell containing the range
% the number of carbon atoms of the i-th surrogate mixture component can vary within

% Outputs:
% 1) prob: logarithmic probability density of the discrete uniform distribution

% Calculate the total number of possible outcomes
num_outcomes = prod(cellfun(@numel, n_ranges));

% Check if the input values for n are within their respective ranges of existence and are integers
for i = 1:numel(n)
    if ~ismember(n(i), n_ranges{i})
        prob = -Inf;
        return
    end
end

% Compute the logarithmic probability density
prob = -log(num_outcomes);
end
