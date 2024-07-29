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

% ------------------------------------------------------------------------%
% Contributors/Copyright
% 2024 Jacopo Liberatori, jacopo.liberatori@uniroma1.it
% 2024 Davide Cavalieri,  davide.cavalieri@uniroma1.it
% Department of Mechanical and Aerospace Engineering (DIMA)
% Sapienza University of Rome
% ------------------------------------------------------------------------%

% Calculate the total number of possible outcomes
num_outcomes = prod(cellfun(@(r) (max(r) - min(r) + 1), n_ranges));

% Check if the input values for n are within their respective ranges of existence and are integers
for i = 1:numel(n)
    if n(i) < min(n_ranges{i}) || n(i) > max(n_ranges{i}) || n(i) ~= floor(n(i))
        prob = -Inf;
        return
    end
end

% Compute the logarithmic probability density
prob = -log(num_outcomes);
end