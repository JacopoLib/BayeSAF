function prob = dirichlet_pdf(x, alpha, LowerBound, UpperBound)

% DESCRIPTION:
% The dirichlet_pdf function returns the logarithmic probability density
% of the symmetric Dirichlet distribution for the set of molar fractions.

% Auxiliary parameters:
% Nc: number of surrogate mixture components  [-]

% Inputs:
% 1) x         : (1 x Nc-1) array containing the molar fractions of the surrogate components  [-]
% 2) alpha     : concentration parameter of the symmetric Dirichlet distribution  [-]
% 3) LowerBound: (1 x Nc) array, with the i-th element representing the
% lower bound for the range the mole fraction of the i-th surrogate mixture
% 3) UpperBound: (1 x Nc) array, with the i-th element representing the
% upper bound for the range the mole fraction of the i-th surrogate mixture

% Outputs:
% 1) prob: logarithmic probability density of the symmetric Dirichlet distribution

% ------------------------------------------------------------------------%
% Contributors/Copyright
% 2024 Jacopo Liberatori, jacopo.liberatori@uniroma1.it
% 2024 Davide Cavalieri,  davide.cavalieri@uniroma1.it
% Department of Mechanical and Aerospace Engineering (DIMA)
% Sapienza University of Rome
% ------------------------------------------------------------------------%

% Check that the inputs are valid
if sum(x) > 1 - LowerBound(end) || sum(x) < 1 - UpperBound(end)
    prob = -Inf;
    return
end

% Compute the logarithmic probability density
if alpha <= 0
    error('alpha must be positive.')
elseif alpha == 1
    prob = 1;
elseif alpha == 2
    prob = log(6) + sum(log(x)) + log(1 - sum(x));
else
    c = alpha*(alpha+1)*(alpha+2)/6;
    prob = log(c) + (alpha - 1) * sum(log(x)) + (alpha - 1) * log(1 - sum(x));
end
end