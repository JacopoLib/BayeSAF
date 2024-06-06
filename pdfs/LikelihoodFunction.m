function L = LikelihoodFunction(families, data, std, X, n, eta_B_star, variable, pressure, n_ranges)

% DESCRIPTION:
% The LikelihoodFunction function computes a Gaussian approximation for marginalized log-likelihoods.

% Auxiliary parameters:
% Nc: number of surrogate mixture components  [-]
% Nd: number of experimental measurements for the thermophysical property under consideration  [-]

% Inputs:
% 1) families: (1 x Nc) cell array, with the i-th cell containing a string
% that denotes the i-th hydrocarbon family under consideration
% 2) data    : (Nd x 3) array containing the experimental measurements about the
% thermophysical property under consideration (independent variable,
% dependent thermophysical property, standard deviation). Note that "data"
% is a (Nd x 4) array in case the distillation curve is addressed in the
% thermophysical characterization of the surrogate mixture, with the 4th
% column denoting the operating pressure the distillation curve is computed at.
% 3) std     : (Nd x 1) array containing the standard deviation associated with the
% experimental measurements about the thermophysical property under consideration
% 4) X       : (1 x Nc-1) array containing the molar fractions of the surrogate
% components in the liquid phase before the distillation process begins  [-]
% 5) n       : (1 x Nc) array containing the numbers of carbon atoms of the
% surrogate components  [-]
% 6) eta_B_star : (1 x Nc) array containing the topochemical atom indices
% of the surrogate components  [-]
% 7) variable: string denoting the thermophysical property under
% consideration
% 8) pressure: operating pressure  [Pa]
% 9) n_ranges: (1 x Nc) cell array, with the i-th cell containing the range
% the number of carbon atoms of the i-th surrogate mixture component can vary within

% Outputs:
% 1) L: Gaussian log-likelihood

% ------------------------------------------------------------------------%
% Contributors/Copyright
% 2024 Jacopo Liberatori, jacopo.liberatori@uniroma1.it
% 2024 Davide Cavalieri,  davide.cavalieri@uniroma1.it
% Department of Mechanical and Aerospace Engineering (DIMA)
% Sapienza University of Rome
% ------------------------------------------------------------------------%

index_n_eta = zeros(1, numel(families)); % indices of the components
% Initialize empty cell array to store the classes
classes = cell(1, numel(families));

for i = 1:numel(families)
    classes{i} = Hydrocarbons(families{i}, n_ranges{i});
    index_n_eta(i) = findIndex_eta(classes{i}, eta_B_star(i));
end

arrayLikelihood = zeros(size(data, 1), 1);
N = length(data);

if numel(unique(n)) == 1

    L = 0;

else

    if ~strcmp(variable,'distillation')
        for i = 1:size(data, 1)
            phiModel = ThermophysicalProperties_LiquidMixture(variable, data(i,1), X, classes, index_n_eta, pressure);
            arrayLikelihood(i) = log(2.*pi.*std(i)^2) + (phiModel - data(i, 2)).^2./(std(i).^2);
        end
        L = - 1./2*sum(arrayLikelihood);
    elseif strcmp(variable,'distillation')
        [volumeFraction, T_distillation]  = DistillationCurve(X, classes, index_n_eta, pressure);
        % Check if all the temperature values along the distillation curve are real
        if all(isreal(T_distillation)) && all(isreal(volumeFraction))
            for i = 1:size(data, 1)
                T_interpolated = interp1(volumeFraction, T_distillation, data(i, 1), 'linear');
                arrayLikelihood(i) = log(2.*pi.*std(i)^2) + (T_interpolated - data(i, 2)).^2./(std(i).^2);
            end
        end
        L = - 1./2*sum(arrayLikelihood);
    end

end

end