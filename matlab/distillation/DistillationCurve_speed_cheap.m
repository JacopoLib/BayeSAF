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
% DistillationCurve_speed_cheap computes partial distillation
% curves for a batch of surrogate mixture samples in parallel using parfor.
% Each sample calls DistillationCurve_cheap, which terminates early once
% the recovered volume fraction exceeds the first experimental data point
% vol1, reducing computational cost during the cheap MCMC stage.
%
% Inputs:
% 1) X_all           : (nSamples x Nc-1) molar fractions  [-]
% 2) classes         : (1 x Nc) cell array of struct arrays — candidate
%                      species per surrogate component
% 3) index_n_eta_all : (nSamples x Nc) species indices in classes  [-]
% 4) pressure        : distillation pressure  [Pa]
% 5) vol1            : first experimental volume fraction target  [-]
%
% Outputs:
% 1) volumeFraction_all  : (nSamples x 101) recovered volume fractions  [-]
% 2) T_distillation_all  : (nSamples x 101) bubble-point temperatures  [K]
% ------------------------------------------------------------------------

function [volumeFraction_all, T_distillation_all] = DistillationCurve_speed_cheap(X_all, classes, index_n_eta_all, pressure, vol1)

nSamples = size(X_all, 1);
nSteps   = 101;

volumeFraction_all = zeros(nSamples, nSteps);
T_distillation_all = zeros(nSamples, nSteps);

parfor s = 1:nSamples
    X       = X_all(s,:);
    idx_eta = index_n_eta_all(s,:);
    [vfrac, Tcurve] = DistillationCurve_cheap(X, classes, idx_eta, pressure, vol1);
    volumeFraction_all(s,:) = vfrac;
    T_distillation_all(s,:) = Tcurve;
end

end
