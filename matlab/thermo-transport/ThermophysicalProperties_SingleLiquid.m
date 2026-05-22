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
% ThermophysicalProperties_SingleLiquid calculates the value of a
% liquid-phase thermophysical property of a pure chemical species at a
% given temperature and pressure using Yaws' polynomial correlations.
% Supported properties: rho, mu, nu, kappa, specificHeat, latentHeat,
% vaporPressure, boilingTemperature, sigma.
%
% Inputs:
% 1) property : string denoting the thermophysical property
% 2) T        : temperature  [K]
% 3) W        : molecular weight  [g/mol]
% 4) Tcr      : critical temperature  [K]
% 5) Arho     : A-coefficient for liquid density (Yaws')
% 6) Brho     : B-coefficient for liquid density (Yaws')
% 7) Crho     : C-coefficient for liquid density (Yaws')
% 8) Amu      : A-coefficient for liquid dynamic viscosity (Yaws')
% 9) Bmu      : B-coefficient for liquid dynamic viscosity (Yaws')
% 10) Cmu     : C-coefficient for liquid dynamic viscosity (Yaws')
% 11) Dmu     : D-coefficient for liquid dynamic viscosity (Yaws')
% 12) Ak      : A-coefficient for liquid thermal conductivity (Yaws')
% 13) Bk      : B-coefficient for liquid thermal conductivity (Yaws')
% 14) Ck      : C-coefficient for liquid thermal conductivity (Yaws')
% 15) Ac      : A-coefficient for liquid specific heat capacity (Yaws')
% 16) Bc      : B-coefficient for liquid specific heat capacity (Yaws')
% 17) Cc      : C-coefficient for liquid specific heat capacity (Yaws')
% 18) Dc      : D-coefficient for liquid specific heat capacity (Yaws')
% 19) Ah      : A-coefficient for latent heat of vaporization (Yaws')
% 20) Bh      : B-coefficient for latent heat of vaporization (Yaws')
% 21) Ap      : A-coefficient for vapour pressure (Yaws')
% 22) Bp      : B-coefficient for vapour pressure (Yaws')
% 23) Cp      : C-coefficient for vapour pressure (Yaws')
% 24) Dp      : D-coefficient for vapour pressure (Yaws')
% 25) Ep      : E-coefficient for vapour pressure (Yaws')
% 26) Asigma  : A-coefficient for liquid surface tension (Yaws')
% 27) Bsigma  : B-coefficient for liquid surface tension (Yaws')
% 28) pressure: pressure  [Pa]
%
% Outputs:
% 1) phi: value of the thermophysical property in SI units
% ------------------------------------------------------------------------

function phi = ThermophysicalProperties_SingleLiquid(property, T, W, Tcr, Arho, Brho, Crho, Amu, Bmu, Cmu, Dmu, Ak, Bk, Ck, Ac, Bc, Cc, Dc, Ah, Bh, Ap, Bp, Cp, Dp, Ep, Asigma, Bsigma, pressure)

if strcmp(property, 'rho')
    phi = 1000 * Arho * Brho.^(-(1 - T./Tcr).^Crho);           % [kg/m^3]
elseif strcmp(property, 'mu')
    phi = 0.001 * 10.^(Amu + Bmu./T + Cmu.*T + Dmu.*T.^2);     % [Pa s]
elseif strcmp(property, 'nu')
    rho = 1000 * Arho * Brho.^(-(1 - T./Tcr).^Crho);
    mu  = 0.001 * 10.^(Amu + Bmu./T + Cmu.*T + Dmu.*T.^2);
    phi = mu ./ rho;                                              % [m^2/s]
elseif strcmp(property, 'kappa')
    phi = Ak + Bk.*T + Ck*T.^2;                                  % [W/m/K]
elseif strcmp(property, 'specificHeat')
    phi = 1000./W * (Ac + Bc.*T + Cc.*T.^2 + Dc.*T.^3);         % [J/kg/K]
elseif strcmp(property, 'latentHeat')
    phi = (Ah.*(1 - T./Tcr).^Bh) ./ W * 1e6;                    % [J/kg]
elseif strcmp(property, 'vaporPressure')
    phi = 10.^(Ap + Bp./T + Cp.*log10(T) + Dp.*T + Ep*T.^2) .* 133.322387415;  % [Pa]
elseif strcmp(property, 'boilingTemperature')
    pSat = @(T) 10.^(Ap + Bp./T + Cp.*log10(T) + Dp.*T + Ep*T.^2) .* 133.322387415 - pressure;
    options = optimset('Display', 'off');
    phi = fsolve(pSat, 500, options);                             % [K]
elseif strcmp(property, 'sigma')
    phi = (Asigma.*(1 - T./Tcr).^Bsigma) ./ 1e3;                % [N/m]
else
    error('Unknown thermophysical property: %s', property);
end

end
