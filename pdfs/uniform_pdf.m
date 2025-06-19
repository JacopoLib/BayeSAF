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

function prob = uniform_pdf(x, LowerBound, UpperBound)

% DESCRIPTION:
% The uniform_pdf function calculates the logarithmic probability density of a continuous
% uniform distribution of the topochemical atom index of each surrogate mixture component.

% Auxiliary parameters:
% Nc: number of surrogate mixture components  [-]

% Inputs:
% 1) x         : (1 x Nc-1) array containing the topochemical atom indices
% of the surrogate components  [-]
% 2) LowerBound: (1 x Nc) array, with the i-th element representing the
% lower bound for the range the topochemical atom index of the i-th 
% surrogate mixture component can vary within
% 3) UpperBound: (1 x Nc) array, with the i-th element representing the
% upper bound for the range the topochemical atom index of the i-th 
% surrogate mixture component can vary within

% Outputs:
% 1) prob: logarithmic probability density of the continuous uniform distribution

% Calculate the logarithmic probability density
prob = 0;
for j = 1:numel(x)
    if x(j) < LowerBound(j) || x(j) > UpperBound(j)
        prob = prob - Inf;
    else
        prob = prob - log(UpperBound(j) - LowerBound(j));
    end
end
end