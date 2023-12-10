% Clear workspace and command window
clc;
clear;

E = 25e9;              % Young's modulus of the matrix material (Pa)
nu = 0.28;             % Poisson's ratio of the matrix material
k_H2O = 2.3e9;         % Stiffness of the pore fluid (Pa)
one = [1 1 1 0 0 0];   % identitiy in Kelvin Mandel Notation

% Define loading scenarios
F_load = linspace(0, 10e6, 20);   % Load factor in the range [0, 10] MPa
phi = linspace(0.25, 0.75, 20);   % Porosity in the range [0.25, 0.75]

% Compute the stiffness tensor for the matrix phase
C_matrix = (1 / (1 + nu)) * [...
    (1-nu)*E/(1-2*nu), nu*E/(1-2*nu), nu*E/(1-2*nu), 0, 0, 0;
    nu*E/(1-2*nu), (1-nu)*E/(1-2*nu), nu*E/(1-2*nu), 0, 0, 0;
    nu*E/(1-2*nu), nu*E/(1-2*nu), (1-nu)*E/(1-2*nu), 0, 0, 0;
    0, 0, 0, E, 0, 0;
    0, 0, 0, 0, E, 0;
    0, 0, 0, 0, 0, E];

K = [1/3, 1/3, 1/3, 0, 0, 0;
    1/3, 1/3, 1/3, 0, 0, 0;
    1/3, 1/3, 1/3, 0, 0, 0;
    0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0];

K_matrix = E / (3 * (1 - 2 * nu));
G_matrix = E / (2 * (1 + nu));

% Compute the Hill tensor
P_sph_iso = inv(C_matrix) * (3 * K_matrix / (3 * K_matrix + 4 * G_matrix) * K + ...
    6 * (K_matrix + 2 * G_matrix) / (5 * (3 * K_matrix + 4 * G_matrix)) * (eye(6) - K));

% Compute the stiffness tensor for the pore phase
C_pore = 3 * k_H2O * K;

% Initialize arrays to store results
strain_energy_density = zeros(length(F_load), length(phi), 4);
hydrostatic_pressure = zeros(length(F_load), length(phi), 4);

% Perform the computations for each loading scenario and porosity value
for i = 1:length(F_load)
    for j = 1:length(phi)

        % Compute strain concentration tensor for the matrix phase
        A_matrix = calculate_A_matrix(C_pore, C_matrix, P_sph_iso, phi(j));
        
        % Compute strain concentration tensor for the pore phase
        A_pore = calculate_A_pore(C_pore, C_matrix, P_sph_iso, phi(j));
        C_RVE = calculate_C_RVE(C_pore, C_matrix, A_matrix, A_pore, phi(j));
        
        % Iterate over all four loading scenarios
        for k = 1:4
            % Compute the macroscopic stress tensor (Kelvin-Mandel notation)
            Sigma_RVE = zeros(1, 6);
            switch k
                case 1
                    Sigma_RVE(1) = -2 * F_load(i);
                case 2
                    Sigma_RVE(1:3) = -F_load(i) * [1, 1, 1];
                case 3
                    Sigma_RVE(6) = sqrt(2) * F_load(i);
                case 4
                    Sigma_RVE(1:6) = F_load(i) * [-3, -1, -1.5, sqrt(2)*1.8, sqrt(2)*0.1, sqrt(2)*0.7];
            end
            
            % Compute the macroscopic strain tensor (Kelvin-Mandel notation)
            E_RVE = inv(C_RVE) * Sigma_RVE';

            % Compute the strain energy density for the matrix material
            epsilon_matrix = A_matrix * E_RVE;
            strain_energy_density(i, j, k) = 0.5 * epsilon_matrix' * C_matrix * epsilon_matrix;

            % Compute the hydrostatic pressure in the pore fluid

            % Compute strain concentration tensor for the matrix phase
            A_matrix_drained = calculate_A_matrix(zeros(6,6), C_matrix, P_sph_iso, phi(j));
        
            % Compute strain concentration tensor for the pore phase
            A_pore_drained = calculate_A_pore(zeros(6,6), C_matrix, P_sph_iso, phi(j));
            C_RVE_drained = calculate_C_RVE(zeros(6,6), C_matrix, A_matrix_drained, A_pore, phi(j));

            b_pore = (phi(j) * one * A_pore_drained)';
            N = 1 / (one * inv(C_matrix) * (b_pore - (phi(j) * one')));

            C_RVE_undrained = C_RVE_drained + inv(phi(j) / k_H2O + 1 / N) * (b_pore * b_pore');
            E_RVE = inv(C_RVE_undrained) * Sigma_RVE';

            p = - inv(phi(j) / k_H2O + 1/N) * b_pore' * E_RVE;
            hydrostatic_pressure(i, j, k) = p;
        end
    end
end

% Plot strain energy density as a function of load factor and porosity
figure;
for k = 1:4
    subplot(2, 2, k);
    surf(phi, F_load / 1e6, squeeze(strain_energy_density(:, :, k)));
    ylabel('Load Factor (MPa)');
    xlabel('Porosity');
    zlabel('Strain Energy Density (J/m^3)');
    title(sprintf('Loading Scenario %d', k));
end

% Plot hydrostatic pressure as a function of load factor and porosity
figure;
for k = 1:4
    subplot(2, 2, k);
    surf(phi, F_load / 1e6, squeeze(hydrostatic_pressure(:, :, k)));
    ylabel('Load Factor (MPa)');
    xlabel('Porosity');
    zlabel('Hydrostatic Pressure (Pa)');
    title(sprintf('Loading Scenario %d', k));
end

% Plot strain energy density and hydrostatic pressure in one figure
figure;

% Plot strain energy density for different load factors
subplot(1, 2, 1);
hold on;
for k = 1:4
    plot(phi, squeeze(strain_energy_density(20, :, k)), 'DisplayName', sprintf('Load Scenario %d', k));
end
hold off;
xlabel('Porosity');
ylabel('Strain Energy Density (J/m^3)');
title('Strain Energy Density vs. Porosity');
legend('Location', 'best');

% Plot hydrostatic pressure for different load factors
subplot(1, 2, 2);
hold on;
for k = 1:4
    plot(phi, squeeze(hydrostatic_pressure(20, :, k)), 'DisplayName', sprintf('Load Scenario %d', k));
end
hold off;
xlabel('Porosity');
ylabel('Hydrostatic Pressure (Pa)');
title('Hydrostatic Pressure vs. Porosity');
legend('Location', 'best');

function A_pore = calculate_A_pore(C_pore, C_matrix, P_sph_iso, phi)
    A_pore = inv(eye(6) + P_sph_iso * (C_pore - C_matrix)) * inv(phi * inv(eye(6) + P_sph_iso * (C_pore - C_matrix)) + (1 - phi) * eye(6));
end

function A_matrix = calculate_A_matrix(C_pore, C_matrix, P_sph_iso, phi)
       A_matrix = eye(6) * inv(phi * inv(eye(6) + P_sph_iso * (C_pore - C_matrix)) + (1 - phi) * eye(6));
end

function C_RVE = calculate_C_RVE(C_pore, C_matrix, A_matrix, A_pore, phi)
    C_RVE = phi * C_pore * A_pore + (1 - phi) * C_matrix * A_matrix;
end
