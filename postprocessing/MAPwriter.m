function species_MAP = MAPwriter(families, numComponents, x_MAP, nc_MAP, eta_B_star_MAP, fuel_name)

% DESCRIPTION:
% The MAPwriter function writes a text file and display information about
% the maximum a posteriori (MAP) surrogate.

% Inputs:
% 1) families          : (1 x numComponents) cell array, with the i-th cell containing a string
% that denotes the i-th hydrocarbon family under consideration
% 2) numComponents: number of surrogate mixture components  [-]
% 3) x_MAP: (1 x numComponents) array containing the molar fractions characterizing the MAP surrogate components
% 4) nc_MAP: (1 x numComponents) array containing the numbers of carbon atoms characterizing the MAP surrogate components
% 5) eta_B_star_MAP: (1 x numComponents) array containing the topochemical atom indices characterizing the MAP surrogate components
% 6) fuel_name     : string denoting the real fuel name

% Outputs:
% 1) species_MAP: (1 x numComponents) cell array, with the i-th cell containing a
% string that denotes the name of the i-th species in the MAP surrogate

% ------------------------------------------------------------------------%
% Contributors/Copyright
% 2024 Jacopo Liberatori, jacopo.liberatori@uniroma1.it
% 2024 Davide Cavalieri,  davide.cavalieri@uniroma1.it
% Department of Mechanical and Aerospace Engineering (DIMA)
% Sapienza University of Rome
% ------------------------------------------------------------------------%

species_MAP = cell(1, numComponents); % name of MAP species
Wi_MAP = zeros(1, numComponents); % initializing the array containing the molecular weights of MAP components
nc_idx = zeros(1000, numComponents);

for l = 1:numComponents
    if strcmp(families(l), 'nparaffins')
        data = 'nparaffins/nparaffins_MatlabData.csv';
        nparaffins = readtable(data, 'VariableNamingRule', 'preserve');
        nc_idx(1:numel(find(nparaffins.nC == nc_MAP(l))),l) = find(nparaffins.nC == nc_MAP(l));
        eta_B_star_diff = abs(nparaffins.eta_B_star_norm(nc_idx(1:numel(find(nparaffins.nC == nc_MAP(l))),l)) - eta_B_star_MAP(l));
        min_eta_B_star_diff = min(abs(nparaffins.eta_B_star_norm(nc_idx(1:numel(find(nparaffins.nC == nc_MAP(l))),l)) - eta_B_star_MAP(l)));
        row_MAP_idx(l) = find(eta_B_star_diff == min_eta_B_star_diff, 1);

        names = 'nparaffins/nparaffins_Names.csv';
        nparaffins_names = readtable(names, 'VariableNamingRule', 'preserve');

        species_MAP(l) = nparaffins_names.Name(nc_idx(1,l) + row_MAP_idx(l) - 1);
        Wi_MAP(l) = nparaffins.W(nc_idx(1,l) + row_MAP_idx(l) - 1);
    elseif strcmp(families(l), 'isoparaffins')
        data = 'isoparaffins/isoparaffins_MatlabData.csv';
        isoparaffins = readtable(data, 'VariableNamingRule', 'preserve');
        nc_idx(1:numel(find(isoparaffins.nC == nc_MAP(l))),l) = find(isoparaffins.nC == nc_MAP(l));
        eta_B_star_diff = abs(isoparaffins.eta_B_star_norm(nc_idx(1:numel(find(isoparaffins.nC == nc_MAP(l))),l)) - eta_B_star_MAP(l));
        min_eta_B_star_diff = min(abs(isoparaffins.eta_B_star_norm(nc_idx(1:numel(find(isoparaffins.nC == nc_MAP(l))),l)) - eta_B_star_MAP(l)));
        row_MAP_idx(l) = find(eta_B_star_diff == min_eta_B_star_diff, 1);

        names = 'isoparaffins/isoparaffins_Names.csv';
        isoparaffins_names = readtable(names, 'VariableNamingRule', 'preserve');

        species_MAP(l) = isoparaffins_names.Name(nc_idx(1,l) + row_MAP_idx(l) - 1);
        Wi_MAP(l) = isoparaffins.W(nc_idx(1,l) + row_MAP_idx(l) - 1);
    elseif strcmp(families(l), 'cycloparaffins')
        data = 'cycloparaffins/cycloparaffins_MatlabData.csv';
        cycloparaffins = readtable(data, 'VariableNamingRule', 'preserve');
        nc_idx(1:numel(find(cycloparaffins.nC == nc_MAP(l))),l) = find(cycloparaffins.nC == nc_MAP(l));
        eta_B_star_diff = abs(cycloparaffins.eta_B_star_norm(nc_idx(1:numel(find(cycloparaffins.nC == nc_MAP(l))),l)) - eta_B_star_MAP(l));
        min_eta_B_star_diff = min(abs(cycloparaffins.eta_B_star_norm(nc_idx(1:numel(find(cycloparaffins.nC == nc_MAP(l))),l)) - eta_B_star_MAP(l)));
        row_MAP_idx(l) = find(eta_B_star_diff == min_eta_B_star_diff, 1);

        names = 'cycloparaffins/cycloparaffins_Names.csv';
        cycloparaffins_names = readtable(names, 'VariableNamingRule', 'preserve');

        species_MAP(l) = cycloparaffins_names.Name(nc_idx(1,l) + row_MAP_idx(l) - 1);
        Wi_MAP(l) = cycloparaffins.W(nc_idx(1,l) + row_MAP_idx(l) - 1);
    elseif strcmp(families(l), 'dicycloparaffins')
        data = 'dicycloparaffins/dicycloparaffins_MatlabData.csv';
        dicycloparaffins = readtable(data, 'VariableNamingRule', 'preserve');
        nc_idx(1:numel(find(dicycloparaffins.nC == nc_MAP(l))),l) = find(dicycloparaffins.nC == nc_MAP(l));
        eta_B_star_diff = abs(dicycloparaffins.eta_B_star_norm(nc_idx(1:numel(find(dicycloparaffins.nC == nc_MAP(l))),l)) - eta_B_star_MAP(l));
        min_eta_B_star_diff = min(abs(dicycloparaffins.eta_B_star_norm(nc_idx(1:numel(find(dicycloparaffins.nC == nc_MAP(l))),l)) - eta_B_star_MAP(l)));
        row_MAP_idx(l) = find(eta_B_star_diff == min_eta_B_star_diff, 1);

        names = 'dicycloparaffins/dicycloparaffins_Names.csv';
        dicycloparaffins_names = readtable(names, 'VariableNamingRule', 'preserve');

        species_MAP(l) = dicycloparaffins_names.Name(nc_idx(1,l) + row_MAP_idx(l) - 1);
        Wi_MAP(l) = dicycloparaffins.W(nc_idx(1,l) + row_MAP_idx(l) - 1);
    elseif strcmp(families(l), 'alkylbenzenes')
        data = 'alkylbenzenes/alkylbenzenes_MatlabData.csv';
        alkylbenzenes = readtable(data, 'VariableNamingRule', 'preserve');
        nc_idx(1:numel(find(alkylbenzenes.nC == nc_MAP(l))),l) = find(alkylbenzenes.nC == nc_MAP(l));
        eta_B_star_diff = abs(alkylbenzenes.eta_B_star_norm(nc_idx(1:numel(find(alkylbenzenes.nC == nc_MAP(l))),l)) - eta_B_star_MAP(l));
        min_eta_B_star_diff = min(abs(alkylbenzenes.eta_B_star_norm(nc_idx(1:numel(find(alkylbenzenes.nC == nc_MAP(l))),l)) - eta_B_star_MAP(l)));
        row_MAP_idx(l) = find(eta_B_star_diff == min_eta_B_star_diff, 1);

        names = 'alkylbenzenes/alkylbenzenes_Names.csv';
        alkylbenzenes_names = readtable(names, 'VariableNamingRule', 'preserve');

        species_MAP(l) = alkylbenzenes_names.Name(nc_idx(1,l) + row_MAP_idx(l) - 1);
        Wi_MAP(l) = alkylbenzenes.W(nc_idx(1,l) + row_MAP_idx(l) - 1);
    elseif strcmp(families(l), 'alkylnaphtalenes')
        data = 'alkylnaphtalenes/alkylnaphtalenes_MatlabData.csv';
        alkylnaphtalenes = readtable(data, 'VariableNamingRule', 'preserve');
        nc_idx(1:numel(find(alkylnaphtalenes.nC == nc_MAP(l))),l) = find(alkylnaphtalenes.nC == nc_MAP(l));
        eta_B_star_diff = abs(alkylnaphtalenes.eta_B_star_norm(nc_idx(1:numel(find(alkylnaphtalenes.nC == nc_MAP(l))),l)) - eta_B_star_MAP(l));
        min_eta_B_star_diff = min(abs(alkylnaphtalenes.eta_B_star_norm(nc_idx(1:numel(find(alkylnaphtalenes.nC == nc_MAP(l))),l)) - eta_B_star_MAP(l)));
        row_MAP_idx(l) = find(eta_B_star_diff == min_eta_B_star_diff, 1);

        names = 'alkylnaphtalenes/alkylnaphtalenes_Names.csv';
        alkylnaphtalenes_names = readtable(names, 'VariableNamingRule', 'preserve');

        species_MAP(l) = alkylnaphtalenes_names.Name(nc_idx(1,l) + row_MAP_idx(l) - 1);
        Wi_MAP(l) = alkylnaphtalenes.W(nc_idx(1,l) + row_MAP_idx(l) - 1);
    elseif strcmp(families(l), 'cycloaromatics')
        data = 'cycloaromatics/cycloaromatics_MatlabData.csv';
        cycloaromatics = readtable(data, 'VariableNamingRule', 'preserve');
        nc_idx(1:numel(find(cycloaromatics.nC == nc_MAP(l))),l) = find(cycloaromatics.nC == nc_MAP(l));
        eta_B_star_diff = abs(cycloaromatics.eta_B_star_norm(nc_idx(1:numel(find(cycloaromatics.nC == nc_MAP(l))),l)) - eta_B_star_MAP(l));
        min_eta_B_star_diff = min(abs(cycloaromatics.eta_B_star_norm(nc_idx(1:numel(find(cycloaromatics.nC == nc_MAP(l))),l)) - eta_B_star_MAP(l)));
        row_MAP_idx(l) = find(eta_B_star_diff == min_eta_B_star_diff, 1);

        names = 'cycloaromatics/cycloaromatics_Names.csv';
        cycloaromatics_names = readtable(names, 'VariableNamingRule', 'preserve');

        species_MAP(l) = cycloaromatics_names.Name(nc_idx(1,l) + row_MAP_idx(l) - 1);
        Wi_MAP(l) = cycloaromatics.W(nc_idx(1,l) + row_MAP_idx(l) - 1);
    end
end

% Calculate mass fractions of MAP components and mean molecular weight of the MAP surrogate mixture
y_MAP = mol2mass(x_MAP, Wi_MAP);
W_MAP = molWeight(x_MAP, Wi_MAP);

% Compute the molecular formula of the MAP surrogate mixture
nCarbon = round(sum(x_MAP.*nc_MAP),2); % number of carbon atoms denoting the molecular formula
nHydrogen = 0.0; % number of hydrogen atoms denoting the molecular formula
for l = 1:numComponents
    if strcmp(families(l), 'nparaffins') || strcmp(families(l), 'isoparaffins')
        nHydrogen = nHydrogen + x_MAP(l)*(2*nc_MAP(l)+2);
    elseif strcmp(families(l), 'cycloparaffins')
        nHydrogen = nHydrogen + x_MAP(l)*2*nc_MAP(l);
    elseif strcmp(families(l), 'dicycloparaffins')
        nHydrogen = nHydrogen + x_MAP(l)*(2*nc_MAP(l)-2);
    elseif strcmp(families(l), 'alkylbenzenes')
        nHydrogen = nHydrogen + x_MAP(l)*(2*nc_MAP(l)-6);
    elseif strcmp(families(l), 'alkylnaphtalenes')
        nHydrogen = nHydrogen + x_MAP(l)*(2*nc_MAP(l)-12);
    elseif strcmp(families(l), 'cycloaromatics')
        nHydrogen = nHydrogen + x_MAP(l)*(2*nc_MAP(l)-8);
    end
end

nHydrogen = round(nHydrogen, 2);

% Hydrogen-to-carbon ratio (H/C)
HC = round(nHydrogen/nCarbon, 3);

% +++ Calculate the pseudocritical temperature (Tc_m) and pressure (Pc_m) of the MAP surrogate mixture according to Lee-Kesler mixing rule +++ %

Tc_m = 0.0; % initialize critical temperature
Vc_m = 0.0; % initialize critical volume
omega_m = 0.0; % initialize acentric factor

for i = 1:numComponents

    for j = 1:numComponents
        
        % Compute critical volume and temperature of the i-th and j-th MAP surrogate mixture components
        if strcmp(families(i), 'nparaffins')
            Vc_i = nparaffins.Vc(nc_idx(i,1) + row_MAP_idx(i) - 1);
            Tc_i = nparaffins.Tc(nc_idx(i,1) + row_MAP_idx(i) - 1);
        elseif strcmp(families(i), 'isoparaffins')
            Vc_i = isoparaffins.Vc(nc_idx(i,1) + row_MAP_idx(i) - 1);
            Tc_i = isoparaffins.Tc(nc_idx(i,1) + row_MAP_idx(i) - 1);
        elseif strcmp(families(i), 'cycloparaffins')
            Vc_i = cycloparaffins.Vc(nc_idx(i,1) + row_MAP_idx(i) - 1);
            Tc_i = cycloparaffins.Tc(nc_idx(i,1) + row_MAP_idx(i) - 1);
        elseif strcmp(families(i), 'dicycloparaffins')
            Vc_i = dicycloparaffins.Vc(nc_idx(i,1) + row_MAP_idx(i) - 1);
            Tc_i = dicycloparaffins.Tc(nc_idx(i,1) + row_MAP_idx(i) - 1);
        elseif strcmp(families(i), 'alkylbenzenes')
            Vc_i = alkylbenzenes.Vc(nc_idx(i,1) + row_MAP_idx(i) - 1);
            Tc_i = alkylbenzenes.Tc(nc_idx(i,1) + row_MAP_idx(i) - 1);
        elseif strcmp(families(i), 'alkylnaphtalenes')
            Vc_i = alkylnaphtalenes.Vc(nc_idx(i,1) + row_MAP_idx(i) - 1);
            Tc_i = alkylnaphtalenes.Tc(nc_idx(i,1) + row_MAP_idx(i) - 1);
        elseif strcmp(families(i), 'cycloaromatics')
            Vc_i = cycloaromatics.Vc(nc_idx(i,1) + row_MAP_idx(i) - 1);
            Tc_i = cycloaromatics.Tc(nc_idx(i,1) + row_MAP_idx(i) - 1);
        end

        if strcmp(families(j), 'nparaffins')
            Vc_j = nparaffins.Vc(nc_idx(j,1) + row_MAP_idx(j) - 1);
            Tc_j = nparaffins.Tc(nc_idx(j,1) + row_MAP_idx(j) - 1);
        elseif strcmp(families(j), 'isoparaffins')
            Vc_j = isoparaffins.Vc(nc_idx(j,1) + row_MAP_idx(j) - 1);
            Tc_j = isoparaffins.Tc(nc_idx(j,1) + row_MAP_idx(j) - 1);
        elseif strcmp(families(j), 'cycloparaffins')
            Vc_j = cycloparaffins.Vc(nc_idx(j,1) + row_MAP_idx(j) - 1);
            Tc_j = cycloparaffins.Tc(nc_idx(j,1) + row_MAP_idx(j) - 1);
        elseif strcmp(families(j), 'dicycloparaffins')
            Vc_j = dicycloparaffins.Vc(nc_idx(j,1) + row_MAP_idx(j) - 1);
            Tc_j = dicycloparaffins.Tc(nc_idx(j,1) + row_MAP_idx(j) - 1);
        elseif strcmp(families(j), 'alkylbenzenes')
            Vc_j = alkylbenzenes.Vc(nc_idx(j,1) + row_MAP_idx(j) - 1);
            Tc_j = alkylbenzenes.Tc(nc_idx(j,1) + row_MAP_idx(j) - 1);
        elseif strcmp(families(j), 'alkylnaphtalenes')
            Vc_j = alkylnaphtalenes.Vc(nc_idx(j,1) + row_MAP_idx(j) - 1);
            Tc_j = alkylnaphtalenes.Tc(nc_idx(j,1) + row_MAP_idx(j) - 1);
        elseif strcmp(families(j), 'cycloaromatics')
            Vc_j = cycloaromatics.Vc(nc_idx(j,1) + row_MAP_idx(j) - 1);
            Tc_j = cycloaromatics.Tc(nc_idx(j,1) + row_MAP_idx(j) - 1);
        end

        Vc_ij = 1./8*(Vc_i.^(1./3) + Vc_j.^(1./3)).^3;
        Tc_ij = (Tc_i.*Tc_j).^(1./2);

        Vc_m = Vc_m + x_MAP(i).*x_MAP(j).*Vc_ij;
        Tc_m = Tc_m + x_MAP(i).*x_MAP(j).*Vc_ij.^(1./4).*Tc_ij;

    end
    
    % Compute acentric factor of the i-th MAP surrogate mixture component
    if strcmp(families(i), 'nparaffins')
            omega_i = nparaffins.omega(nc_idx(i,1) + row_MAP_idx(i) - 1);
        elseif strcmp(families(i), 'isoparaffins')
            omega_i = isoparaffins.omega(nc_idx(i,1) + row_MAP_idx(i) - 1);
        elseif strcmp(families(i), 'cycloparaffins')
            omega_i = cycloparaffins.omega(nc_idx(i,1) + row_MAP_idx(i) - 1);
        elseif strcmp(families(i), 'dicycloparaffins')
            omega_i = dicycloparaffins.omega(nc_idx(i,1) + row_MAP_idx(i) - 1);
        elseif strcmp(families(i), 'alkylbenzenes')
            omega_i = alkylbenzenes.omega(nc_idx(i,1) + row_MAP_idx(i) - 1);
        elseif strcmp(families(i), 'alkylnaphtalenes')
            omega_i = alkylnaphtalenes.omega(nc_idx(i,1) + row_MAP_idx(i) - 1);
        elseif strcmp(families(i), 'cycloaromatics')
            omega_i = cycloaromatics.omega(nc_idx(i,1) + row_MAP_idx(i) - 1);
    end

    omega_m = omega_m + x_MAP(i).*omega_i;

end

Tc_m = Tc_m./(Vc_m.^(1./4)); % pseudocritical temperature of the MAP surrogate mixture [K]

R_univ = 8.314; % universal gas constant [J/mol/K]

Pc_m = (0.2905 - 0.085.*omega_m).*(R_univ.*Tc_m)./Vc_m./1e+05; % pseudocritical pressure of the MAP surrogate mixture [bar]


fid = fopen('MAP.txt', 'w'); % Open the file for writing

% Write the ASCII art header
fprintf(fid, '    ____                  _____ ___    ______                            \n');
fprintf(fid, '   / __ )____ ___  _____ / ___//   |  / ____/                            \n');
fprintf(fid, '  / __  / __ `/ / / / _ \\\\__ \\/ /| | / /_                             \n');
fprintf(fid, ' / /_/ / /_/ / /_/ /  __/__/ / ___ |/ __/                                \n');
fprintf(fid, '/_____/\\__,_/\\__, /\\___/____/_/  |_/_/                                \n');
fprintf(fid, '            /____/                                                       \n');
fprintf(fid, '\n'); % Add an extra newline for separation
fprintf(fid, '-----------------------------------------------------------------------  \n');
fprintf(fid, 'Contributors/Copyright                                                   \n');
fprintf(fid, '2024 Jacopo Liberatori, jacopo.liberatori@uniroma1.it                    \n');
fprintf(fid, '2024 Davide Cavalieri,  davide.cavalieri@uniroma1.it                     \n');
fprintf(fid, 'Department of Mechanical and Aerospace Engineering (DIMA)                \n');
fprintf(fid, 'Sapienza University of Rome                                              \n');
fprintf(fid, '-----------------------------------------------------------------------  \n');
fprintf(fid, '\n'); % Add an extra newline for separation

% Write information to the file
fprintf(fid, ['+++++ MAP surrogate for ', fuel_name, ' +++++\n']);
fprintf(fid, 'Number of MAP components: %s\n', num2str(numComponents));
fprintf(fid, 'MAP hydrocarbon classes: %s\n', strjoin(families, ', '));
fprintf(fid, 'MAP hydrocarbons: %s\n', strjoin(species_MAP, ', '));
fprintf(fid, ['MAP molecular weight: ', num2str(W_MAP), ' g/mol\n']);
fprintf(fid, ['MAP molecular formula: C', num2str(nCarbon), ' H', num2str(nHydrogen), '\n']);
fprintf(fid, 'MAP hydrogen-to-carbon ratio: %s\n', num2str(HC));
fprintf(fid, 'MAP molar fractions: %s\n', num2str(x_MAP));
fprintf(fid, 'MAP mass fractions: %s\n', num2str(y_MAP));
fprintf(fid, 'MAP carbon atoms: %s\n', num2str(nc_MAP));
fprintf(fid, 'MAP topochemical atom indices: %s\n', num2str(eta_B_star_MAP));
fprintf(fid, '\n'); % Add an extra newline for separation
fprintf(fid, ['MAP pseudocritical temperature according to Lee-Kesler mixing rule: ', num2str(Tc_m), ' K\n']);
fprintf(fid, ['MAP pseudocritical pressure according to Lee-Kesler mixing rule: ', num2str(Pc_m), ' bar\n']);

fclose(fid); % Close the file

end