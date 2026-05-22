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
% ThermophysicalProperties_LiquidMixture_speed evaluates a liquid-mixture
% thermophysical property for multiple MCMC samples in a vectorized
% fashion. Mixing rules: volume-weighted (rho), log-linear (mu),
% mass-weighted (kappa, specificHeat, latentHeat, vaporPressure),
% molar-weighted (sigma).
%
% Inputs:
% 1) property        : string denoting the thermophysical property
% 2) T               : (N_exp x 1) experimental temperatures  [K]
% 3) X_all           : (N_samples x Nc-1) molar fractions  [-]
% 4) classes         : (1 x Nc) cell array of struct arrays
% 5) index_n_eta_all : (N_samples x Nc) species indices  [-]
% 6) pressure        : pressure  [Pa]
%
% Outputs:
% 1) phi: (N_samples x N_exp) property values in SI units
% ------------------------------------------------------------------------

function phi = ThermophysicalProperties_LiquidMixture_speed(property, T, X_all, classes, index_n_eta_all, pressure)

N_samples = size(X_all, 1);
N_exp     = length(T);
Nc        = numel(classes);

% Full molar fraction matrix [N_samples x Nc]
X_full = [X_all, 1 - sum(X_all, 2)];

switch property
    case 'rho'
        rho       = zeros(N_samples, N_exp, Nc);
        numerator = zeros(N_samples, N_exp, Nc);

        for i = 1:Nc
            idx_i = index_n_eta_all(:, i);
            for j = 1:N_samples
                sp = classes{i}(idx_i(j));
                rho(j,:,i) = ThermophysicalProperties_SingleLiquid('rho', T, ...
                    sp.molWeight, sp.Tc, ...
                    sp.coeffRho(1), sp.coeffRho(2), sp.coeffRho(3), ...
                    sp.coeffMu(1),  sp.coeffMu(2),  sp.coeffMu(3),  sp.coeffMu(4), ...
                    sp.coeffK(1),   sp.coeffK(2),   sp.coeffK(3), ...
                    sp.coeffCl(1),  sp.coeffCl(2),  sp.coeffCl(3),  sp.coeffCl(4), ...
                    sp.coeffHv(1),  sp.coeffHv(2), ...
                    sp.coeffPsat(1), sp.coeffPsat(2), sp.coeffPsat(3), sp.coeffPsat(4), sp.coeffPsat(5), ...
                    sp.coeffSigma(1), sp.coeffSigma(2), pressure);
                numerator(j,:,i) = X_full(j,i) * sp.molWeight;
            end
        end

        phi = sum(numerator, 3) ./ sum(numerator ./ rho, 3);

    case 'mu'
        mu_log = zeros(N_samples, N_exp, Nc);

        for i = 1:Nc
            idx_i = index_n_eta_all(:, i);
            for j = 1:N_samples
                sp = classes{i}(idx_i(j));
                mu_log(j,:,i) = X_full(j,i) * log(ThermophysicalProperties_SingleLiquid('mu', T, ...
                    sp.molWeight, sp.Tc, ...
                    sp.coeffRho(1), sp.coeffRho(2), sp.coeffRho(3), ...
                    sp.coeffMu(1),  sp.coeffMu(2),  sp.coeffMu(3),  sp.coeffMu(4), ...
                    sp.coeffK(1),   sp.coeffK(2),   sp.coeffK(3), ...
                    sp.coeffCl(1),  sp.coeffCl(2),  sp.coeffCl(3),  sp.coeffCl(4), ...
                    sp.coeffHv(1),  sp.coeffHv(2), ...
                    sp.coeffPsat(1), sp.coeffPsat(2), sp.coeffPsat(3), sp.coeffPsat(4), sp.coeffPsat(5), ...
                    sp.coeffSigma(1), sp.coeffSigma(2), pressure));
            end
        end

        phi = exp(sum(mu_log, 3));

    case 'nu'
        rho = ThermophysicalProperties_LiquidMixture_speed('rho', T, X_all, classes, index_n_eta_all, pressure);
        mu  = ThermophysicalProperties_LiquidMixture_speed('mu',  T, X_all, classes, index_n_eta_all, pressure);
        phi = mu ./ rho;

    case {'kappa', 'specificHeat', 'latentHeat', 'vaporPressure'}
        numerator   = zeros(N_samples, N_exp, Nc);
        denominator = zeros(N_samples, Nc);

        for i = 1:Nc
            idx_i = index_n_eta_all(:, i);
            for j = 1:N_samples
                sp = classes{i}(idx_i(j));
                numerator(j,:,i) = X_full(j,i) * sp.molWeight * ...
                    ThermophysicalProperties_SingleLiquid(property, T, ...
                    sp.molWeight, sp.Tc, ...
                    sp.coeffRho(1), sp.coeffRho(2), sp.coeffRho(3), ...
                    sp.coeffMu(1),  sp.coeffMu(2),  sp.coeffMu(3),  sp.coeffMu(4), ...
                    sp.coeffK(1),   sp.coeffK(2),   sp.coeffK(3), ...
                    sp.coeffCl(1),  sp.coeffCl(2),  sp.coeffCl(3),  sp.coeffCl(4), ...
                    sp.coeffHv(1),  sp.coeffHv(2), ...
                    sp.coeffPsat(1), sp.coeffPsat(2), sp.coeffPsat(3), sp.coeffPsat(4), sp.coeffPsat(5), ...
                    sp.coeffSigma(1), sp.coeffSigma(2), pressure);
                denominator(j,i) = X_full(j,i) * sp.molWeight;
            end
        end

        phi = sum(numerator, 3) ./ sum(denominator, 2);

    case 'sigma'
        sigma_mat = zeros(N_samples, N_exp, Nc);
        Xs        = zeros(N_samples, N_exp, Nc);

        for i = 1:Nc
            idx_i = index_n_eta_all(:, i);
            for j = 1:N_samples
                sp = classes{i}(idx_i(j));
                sigma_mat(j,:,i) = ThermophysicalProperties_SingleLiquid('sigma', T, ...
                    sp.molWeight, sp.Tc, ...
                    sp.coeffRho(1), sp.coeffRho(2), sp.coeffRho(3), ...
                    sp.coeffMu(1),  sp.coeffMu(2),  sp.coeffMu(3),  sp.coeffMu(4), ...
                    sp.coeffK(1),   sp.coeffK(2),   sp.coeffK(3), ...
                    sp.coeffCl(1),  sp.coeffCl(2),  sp.coeffCl(3),  sp.coeffCl(4), ...
                    sp.coeffHv(1),  sp.coeffHv(2), ...
                    sp.coeffPsat(1), sp.coeffPsat(2), sp.coeffPsat(3), sp.coeffPsat(4), sp.coeffPsat(5), ...
                    sp.coeffSigma(1), sp.coeffSigma(2), pressure);
                Xs(j,:,i) = X_full(j,i);
            end
        end

        phi = sum(sigma_mat .* Xs, 3);

    otherwise
        error('Unknown thermophysical property: %s', property);
end

end
