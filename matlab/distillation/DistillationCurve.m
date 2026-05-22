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
% DistillationCurve computes the distillation curve of a surrogate
% fuel mixture. It integrates the Rachford-Rice flash equation step by step
% from 0 to 100% recovered mole fraction, finding the bubble-point
% temperature at each step via fzero with adaptive bracket widening.
% Species thermophysical data (Psat, rho) are pre-extracted once before the
% main loop for maximum efficiency.
%
% Inputs:
% 1) X0          : (1 x Nc-1) array of initial liquid molar fractions  [-]
% 2) classes     : (1 x Nc) cell array of struct arrays — candidate
%                  species per surrogate component
% 3) index_n_eta : (1 x Nc) array of species indices in classes  [-]
% 4) pressure    : distillation pressure  [Pa]
%
% Outputs:
% 1) volumeFraction  : (1 x 101) recovered volume fraction  [-]
% 2) T_distillation  : (1 x 101) bubble-point temperatures  [K]
% ------------------------------------------------------------------------

function [volumeFraction, T_distillation] = DistillationCurve(X0, classes, index_n_eta, pressure)

% -----------------------
% Setup and preallocation
% -----------------------
Nc = numel(classes);
assert(numel(index_n_eta) == Nc, 'index_n_eta must be 1xNc.');

X = X0(:).';                   % row vector (1 x Nc-1)
mol = 0:100;
nSteps = numel(mol);

T_distillation = zeros(1, nSteps);
volumeFraction = zeros(1, nSteps);
pSat           = zeros(Nc,1);
K              = zeros(Nc,1);
Vtilde_i       = zeros(Nc,1);
rho_i          = zeros(Nc,1);

% ---------------------------------------
% Cache all species data once (big win!)
% ---------------------------------------
molWeight  = zeros(Nc,1);
Tc         = zeros(Nc,1);
coeffRho   = zeros(Nc,3);
coeffMu    = zeros(Nc,4);
coeffK     = zeros(Nc,3);
coeffCl    = zeros(Nc,4);
coeffHv    = zeros(Nc,2);
coeffPsat  = zeros(Nc,5);
coeffSigma = zeros(Nc,2);

for j = 1:Nc
    s = classes{j}(index_n_eta(j));
    molWeight(j)    = s.molWeight;
    Tc(j)           = s.Tc;
    coeffRho(j,:)   = s.coeffRho(:).';
    coeffMu(j,:)    = s.coeffMu(:).';
    coeffK(j,:)     = s.coeffK(:).';
    coeffCl(j,:)    = s.coeffCl(:).';
    coeffHv(j,:)    = s.coeffHv(:).';
    coeffPsat(j,:)  = s.coeffPsat(:).';
    coeffSigma(j,:) = s.coeffSigma(:).';
end

% -------------------------------
% Solver setup & initial guesses
% -------------------------------
beta = 1e-6;
T0   = 300.0;

VLE_fun = @(T,beta_local,X_local) Flash(T, beta_local, classes, index_n_eta, X_local, pressure);

% -----------------------------------
% Step 0: find first temperature (LLE)
% -----------------------------------
temp = find_root_fzero(VLE_fun, T0, beta, X);

T_distillation(1) = temp;

% ------------------------------------------
% First-step volume (V0) using densities at T
% ------------------------------------------
rho_i = eval_props('rho', temp, ...
    molWeight, Tc, coeffRho, coeffMu, coeffK, coeffCl, coeffHv, coeffPsat, coeffSigma, pressure);

vol_i = zeros(Nc,1);
for j = 1:Nc
    if j < Nc
        vol_i(j) = X0(j) * molWeight(j) / rho_i(j);
    else
        vol_i(j) = (1 - sum(X0)) * molWeight(j) / rho_i(j);
    end
end
V0 = sum(vol_i);

% ----------------
% Distillation run
% ----------------
for ii = 1:(nSteps-1)
    beta = 1.0 / (100.0 - mol(ii));

    temp = find_root_fzero(VLE_fun, temp, beta, X);
    T_distillation(ii+1) = temp;

    rho_i = eval_props('rho', temp, ...
        molWeight, Tc, coeffRho, coeffMu, coeffK, coeffCl, coeffHv, coeffPsat, coeffSigma, pressure);

    pSat = eval_props('vaporPressure', temp, ...
        molWeight, Tc, coeffRho, coeffMu, coeffK, coeffCl, coeffHv, coeffPsat, coeffSigma, pressure);

    K = pSat / pressure;

    den = 1 + beta * (K(1:end-1) - 1);
    X   = X ./ den';                         % update liquid composition
    Y   = (K(1:end-1) .* X(:)).';            % vapor composition
    Y_last = 1 - sum(Y);

    % Vectorized vapor-phase volume contribution
    Vtilde_i(1:Nc-1) = (Y(:) .* molWeight(1:Nc-1)) ./ rho_i(1:Nc-1);
    Vtilde_i(Nc)     = Y_last * molWeight(Nc) / rho_i(Nc);

    Vtilde = sum(Vtilde_i);
    volumeFraction(ii+1) = volumeFraction(ii) + Vtilde / V0;
end

end

% ============================================================
% Helper: scalar root finder using fzero with adaptive bracket
% ============================================================
function T = find_root_fzero(VLE_fun, T_guess, beta, X)

T_low  = max(200, T_guess - 100);
T_high = min(1200, T_guess + 100);
F = @(T) VLE_fun(T, beta, X);

try
    f_low  = F(T_low);
    f_high = F(T_high);
    expand = 0;

    while sign(f_low) == sign(f_high) && expand < 3
        expand = expand + 1;
        T_low  = max(150, T_low  - 100);
        T_high = min(1500, T_high + 100);
        f_low  = F(T_low);
        f_high = F(T_high);
    end

    T = fzero(F, [T_low, T_high]);
catch
    T = T_guess;
end

end

% ============================================================
% Local utility: evaluate properties for all components at T
% ============================================================
function out = eval_props(what, T, ...
    molWeight, Tc, coeffRho, coeffMu, coeffK, coeffCl, coeffHv, coeffPsat, coeffSigma, pressure)

Nc  = numel(molWeight);
out = zeros(Nc,1);

switch what
    case 'rho'
        for j = 1:Nc
            out(j) = ThermophysicalProperties_SingleLiquid('rho', T, ...
                molWeight(j), Tc(j), ...
                coeffRho(j,1),  coeffRho(j,2),  coeffRho(j,3), ...
                coeffMu(j,1),   coeffMu(j,2),   coeffMu(j,3),  coeffMu(j,4), ...
                coeffK(j,1),    coeffK(j,2),    coeffK(j,3), ...
                coeffCl(j,1),   coeffCl(j,2),   coeffCl(j,3),  coeffCl(j,4), ...
                coeffHv(j,1),   coeffHv(j,2), ...
                coeffPsat(j,1), coeffPsat(j,2), coeffPsat(j,3), coeffPsat(j,4), coeffPsat(j,5), ...
                coeffSigma(j,1), coeffSigma(j,2), pressure);
        end

    case 'vaporPressure'
        for j = 1:Nc
            out(j) = ThermophysicalProperties_SingleLiquid('vaporPressure', T, ...
                molWeight(j), Tc(j), ...
                coeffRho(j,1),  coeffRho(j,2),  coeffRho(j,3), ...
                coeffMu(j,1),   coeffMu(j,2),   coeffMu(j,3),  coeffMu(j,4), ...
                coeffK(j,1),    coeffK(j,2),    coeffK(j,3), ...
                coeffCl(j,1),   coeffCl(j,2),   coeffCl(j,3),  coeffCl(j,4), ...
                coeffHv(j,1),   coeffHv(j,2), ...
                coeffPsat(j,1), coeffPsat(j,2), coeffPsat(j,3), coeffPsat(j,4), coeffPsat(j,5), ...
                coeffSigma(j,1), coeffSigma(j,2), pressure);
        end

    otherwise
        error('Unknown property request: %s', what);
end

end
