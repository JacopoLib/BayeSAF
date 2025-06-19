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

function phi = ThermophysicalProperties_SingleLiquid(property, T, W, Tcr, Arho, Brho, Crho, Amu, Bmu, Cmu, Dmu, Ak, Bk, Ck, Ac, Bc, Cc, Dc, Ah, Bh, Ap, Bp, Cp, Dp, Ep, Asigma, Bsigma, pressure)

% DESCRIPTION:
% The ThermophysicalProperties_SingleLiquid function calculates the value of a liquid-phase thermophysical
% property of a pure chemical species at given temperature and pressure.

% Inputs:
% 1) property : string denoting the thermophysical property under consideration
% 2) T        : temperature  [K]
% 3) W        : molecular weight  [g/mol]
% 4) Tcr      : critical temperature  [K]
% 5) Arho     : A-coefficient for liquid density computation by Yaws' polynomials
% 6) Brho     : B-coefficient for liquid density computation by Yaws' polynomials
% 7) Crho     : C-coefficient for liquid density computation by Yaws' polynomials
% 8) Amu      : A-coefficient for liquid dynamic viscosity computation by Yaws' polynomials
% 9) Bmu      : B-coefficient for liquid dynamic viscosity computation by Yaws' polynomials
% 10) Cmu     : C-coefficient for liquid dynamic viscosity computation by Yaws' polynomials
% 11) Dmu     : D-coefficient for liquid dynamic viscosity computation by Yaws' polynomials
% 12) Ak      : A-coefficient for liquid thermal conductivity computation by Yaws' polynomials
% 13) Bk      : B-coefficient for liquid thermal conductivity computation by Yaws' polynomials
% 14) Ck      : C-coefficient for liquid thermal conductivity computation by Yaws' polynomials
% 15) Ac      : A-coefficient for liquid specific heat capacity computation by Yaws' polynomials
% 16) Bc      : B-coefficient for liquid specific heat capacity computation by Yaws' polynomials
% 17) Cc      : C-coefficient for liquid specific heat capacity computation by Yaws' polynomials
% 18) Dc      : D-coefficient for liquid specific heat capacity computation by Yaws' polynomials
% 19) Ah      : A-coefficient for latent heat of vaporization computation by Yaws' polynomials
% 20) Bh      : B-coefficient for latent heat of vaporization computation by Yaws' polynomials
% 21) Ap      : A-coefficient for vapor pressure computation by Yaws' polynomials
% 22) Bp      : B-coefficient for vapor pressure computation by Yaws' polynomials
% 23) Cp      : C-coefficient for vapor pressure computation by Yaws' polynomials
% 24) Dp      : D-coefficient for vapor pressure computation by Yaws' polynomials
% 25) Ep      : E-coefficient for vapor pressure computation by Yaws' polynomials
% 26) Asigma  : A-coefficient for liquid surface tension computation by Yaws' polynomials
% 27) Bsigma  : B-coefficient for liquid surface tension computation by Yaws' polynomials
% 28) pressure: pressure  [Pa]

% Outputs:
% 1) phi: value of the thermophysical property under consideration

if strcmp(property,'rho')    % [kg/m^3]
    phi = 1000 * Arho * Brho.^(-(1 - T./Tcr).^Crho);
elseif strcmp(property,'mu')
    phi = 0.001 * 10.^(Amu + Bmu./T + Cmu.*T + Dmu.*T.^2);    % [Pa s]
elseif strcmp(property,'nu')
    rho = 1000 * Arho * Brho.^(-(1 - T./Tcr).^Crho);
    mu = 0.001 * 10.^(Amu + Bmu./T + Cmu.*T + Dmu.*T.^2);
    phi = mu./rho;    % [m^2/s]
elseif strcmp(property,'kappa')
    phi = Ak + Bk.*T + Ck*T.^2;    % [W/m/K]
elseif strcmp(property,'specificHeat')
    phi = 1000./W*(Ac + Bc.*T + Cc.*T.^2 + Dc.*T.^3);    % [J/kg/K]
elseif strcmp(property,'latentHeat')
    phi = (Ah.*(1 - T./Tcr).^Bh)./W*1e+6;    % [J/kg]
elseif strcmp(property,'vaporPressure')
    phi = 10.^(Ap + Bp./T + Cp.*log10(T) + Dp.*T + Ep*T.^2).*133.322387415;    % [Pa]
elseif strcmp(property,'boilingTemperature')
    pSat = @(T) 10.^(Ap + Bp./T + Cp.*log10(T) + Dp.*T + Ep*T.^2).*133.322387415 - pressure;
    options = optimset('Display', 'off');
    phi = fsolve(pSat, 500, options);    % [K]
elseif strcmp(property,'sigma')
    phi = (Asigma.*(1 - T./Tcr).^Bsigma)./1e+3;    % [N/m]
else
    error('Unknown thermophysical property.');
end
end