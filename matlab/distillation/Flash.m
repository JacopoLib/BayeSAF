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
% Flash returns the left-hand side of the Rachford-Rice flash distillation
% equation for the given surrogate mixture at a fixed temperature, recovered
% mole fraction, liquid-phase molar composition, and pressure. Vapour
% pressures of all components are evaluated via Yaws' polynomial
% correlations and used to compute the equilibrium ratios K_i = p_sat,i/p.
%
% Inputs:
% 1) T          : temperature  [K]
% 2) beta       : recovered mole fraction (vapour fraction)  [-]
% 3) classes    : (1 x Nc) cell array of struct arrays — candidate
%                 species per surrogate component
% 4) index_n_eta: (1 x Nc) array of species indices in classes  [-]
% 5) X          : (1 x Nc-1) array of liquid molar fractions  [-]
% 6) pressure   : pressure  [Pa]
%
% Outputs:
% 1) F: left-hand side of the Rachford-Rice flash equation  [-]
% ------------------------------------------------------------------------

function F = Flash(T, beta, classes, index_n_eta, X, pressure)

Nc    = numel(classes);
Apsat = zeros(1, Nc);
Bpsat = zeros(1, Nc);
Cpsat = zeros(1, Nc);
Dpsat = zeros(1, Nc);
Epsat = zeros(1, Nc);
W     = zeros(1, Nc);
Tc    = zeros(1, Nc);
Asigma = zeros(1, Nc);
Bsigma = zeros(1, Nc);

for i = 1:Nc
    c = classes{i}(index_n_eta(i));
    W(i)      = c.molWeight;
    Tc(i)     = c.Tc;
    Apsat(i)  = c.coeffPsat(1);
    Bpsat(i)  = c.coeffPsat(2);
    Cpsat(i)  = c.coeffPsat(3);
    Dpsat(i)  = c.coeffPsat(4);
    Epsat(i)  = c.coeffPsat(5);
    Asigma(i) = c.coeffSigma(1);
    Bsigma(i) = c.coeffSigma(2);
end

% Vectorized vapor pressures for all components
pSat = ThermophysicalProperties_SingleLiquid('vaporPressure', T, W, Tc, ...
    [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],Apsat,Bpsat,Cpsat,Dpsat,Epsat,Asigma,Bsigma, pressure);

K = pSat / pressure;

X_full = [X, 1 - sum(X)];
flash  = X_full .* (K - 1) ./ (1 + beta * (K - 1));
F = sum(flash);

end
