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

% ------------------------------------------------------------------------%
% Contributors/Copyright
% 2024 Jacopo Liberatori, jacopo.liberatori@uniroma1.it
% 2024 Davide Cavalieri,  davide.cavalieri@uniroma1.it
% Department of Mechanical and Aerospace Engineering (DIMA)
% Sapienza University of Rome
% ------------------------------------------------------------------------%

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