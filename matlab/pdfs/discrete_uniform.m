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
% discrete_uniform returns the logarithmic probability density of a
% discrete uniform distribution over the number of carbon atoms of each
% surrogate mixture component. Returns -Inf if any component's carbon
% number falls outside its admissible range.
%
% Inputs:
% 1) n       : (1 x Nc) array of carbon-atom numbers  [-]
% 2) n_ranges: (1 x Nc) cell array of admissible carbon-number ranges  [-]
%
% Outputs:
% 1) prob: log probability density of the discrete uniform distribution  [-]
% ------------------------------------------------------------------------

function prob = discrete_uniform(n, n_ranges)

for i = 1:numel(n)
    if ~ismember(n(i), n_ranges{i})
        prob = -Inf;
        return;
    end
end

num_outcomes = prod(cellfun(@numel, n_ranges));
prob = -log(num_outcomes);

end
