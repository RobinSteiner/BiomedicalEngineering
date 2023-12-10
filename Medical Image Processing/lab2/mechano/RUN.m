clear all;

% Material properties
E = 25e9;          % Young's modulus of the matrix material (GPa)
nu = 0.28;         % Poisson's ratio of the matrix material
nu2 = 2 * nu;      % Poisson's ratio of the matrix material
k_H2O = 2.3e9;     % Stiffness of the pore fluid (GPa)
F_load = 0:10e5:10e6; % Load factor (Pa)
phi = 0.25:0.1:0.75; % Porosity

% Stiffness tensors
C_matrix = 1/(1+nu) * [((1-nu)*E)/(1-nu2), (nu*E)/(1-nu2), (nu*E)/(1-nu2), 0, 0, 0;
                      (nu*E)/(1-nu2), ((1-nu)*E)/(1-nu2), (nu*E)/(1-nu2), 0, 0, 0;
                      (nu*E)/(1-nu2), (nu*E)/(1-nu2), ((1-nu)*E)/(1-nu2), 0, 0, 0;
                      0, 0, 0, E, 0, 0;
                      0, 0, 0, 0, E, 0;
                      0, 0, 0, 0, 0, E];

K = [1/3, 1/3, 1/3, 0, 0, 0;
     1/3, 1/3, 1/3, 0, 0, 0;
     1/3, 1/3, 1/3, 0, 0, 0;
     0, 0, 0, 0, 0, 0;
     0, 0, 0, 0, 0, 0;
     0, 0, 0, 0, 0, 0];

C_pore = 3 * k_H2O * K; % Assuming isotropic pore phase

% Initialize result arrays
psi_results = zeros(4, length(F_load), length(phi), length(F_load));

K_matrix = E / (3 * (1 - nu2));
G_matrix = E / (2 * (1 + nu));
P_iso = inv(C_matrix) * (3 * K_matrix / (3 * K_matrix + 4 * G_matrix) * K + 6 * (K_matrix + 2 * G_matrix) / (5 * (3 * K_matrix + 4 * G_matrix)) * (eye(6) - K));

% Task 1: Compute strain energy density
for i = 1:length(F_load)
    for j = 1:length(phi)
        A_pore = inv(eye(6) + P_iso * (C_pore - C_matrix)) * inv(phi(j) * inv(eye(6) + P_iso * (C_pore - C_matrix)) + (1 - phi(j)) * eye(6));
        A_matrix = eye(6) * inv(phi(j) * inv(eye(6) + P_iso * (C_pore - C_matrix)) + (1 - phi(j)) * eye(6));
        C_RVE = phi(j) * C_pore * A_pore + (1 - phi(j)) * C_matrix * A_matrix;
        
        for k = 1:length(F_load)
            stresstensors = [-2 * F_load(k), 0, 0, 0, 0, 0;
                            -F_load(k), -F_load(k), -F_load(k), 0, 0, 0;  
                            0, 0, 0, 0, 0, sqrt(2) * F_load(k);  
                            -3 * F_load(k), -F_load(k), -1.5 * F_load(k), sqrt(2) * 1.8 * F_load(k), sqrt(2) * 0.1 * F_load(k), sqrt(2) * 0.7 * F_load(k)];
            for l = 1:length(stresstensors(:, 1))
                ERVE = inv(C_RVE) * stresstensors(l, :)';
                epsilon_matrix = A_matrix * ERVE;
                psi_matrix = 0.5 * epsilon_matrix' * C_matrix * epsilon_matrix;
            
                psi_results(l, i, j, k) = psi_matrix;
            end
        end
    end
end


% Initialize result arrays
p_results = zeros(4, length(F_load), length(phi), length(F_load));
% Task 2: Compute pore pressure
for i = 1:length(F_load)
    for j = 1:length(phi)        
        A_pore_drained = inv(eye(6) + P_iso*C_matrix) * inv(phi(j) * inv(eye(6) + P_iso*C_matrix)+(1-phi(j))*eye(6)); % Assuming drained pore phase
        C_RVE_drained = (phi(j) * inv(eye(6) + P_iso *C_matrix) + (1 - phi(j)) * C_matrix * eye(6))*inv(phi(j) * inv(eye(6) + P_iso * C_matrix) + (1 - phi(j)) * eye(6));
        b_pore = ((phi(1) * ones(1,6)) * A_pore_drained)';
        one_over_N = ones(1,6) * inv(C_matrix) * (b_pore-(phi(1)*ones(6,1)));
        C_RVE_undrained = C_RVE_drained + inv(phi(1)/k_H2O + one_over_N) * (b_pore' * b_pore);

        for k = 1:length(F_load)
                    stresstensors = [-2*F_load(k), 0, 0, 0, 0, 0;
                         -F_load(k), -F_load(k), -F_load(k), 0, 0, 0;  
                         0, 0, 0, 0, 0, sqrt(2)*F_load(k);  
                         -3*F_load(k), -F_load(k), -1.5*F_load(k), sqrt(2)*1.8*F_load(k), sqrt(2)*0.1*F_load(k), sqrt(2)*0.7*F_load(k)];
            for l = 1:length(stresstensors(:,1))
                ERVE = inv(C_RVE_undrained) * stresstensors(l,:)';
                pore_pressure = (-inv(phi(1)/k_H2O + one_over_N) * b_pore') * ERVE;
                                
                p_results(l, i, j, k) = pore_pressure;
            end
        end
    end
end

plotResults(psi_results, p_results, F_load, phi);


function plotResults(psi_results, p_results, F_load_values, phi_values)
    % psi_results: Strain energy density results
    % p_results: Pore pressure results
    % phi_values: Array of porosity values
    % F_load_values: Array of load factor values
    
    num_phi = length(phi_values);
    
    % Plot strain energy density
    for l = 1:4
        figure;
        hold on;
        line_handles = gobjects(size(phi_values));
        for j = 1:num_phi
            for i = 1:size(psi_results, 1)
                line_handles(j) = plot(F_load_values, squeeze(psi_results(l, i, j, :)), '-o', 'LineWidth', 1.5);
            end
        end
        hold off;
        xlabel('Load Factor (Pa)');
        ylabel('Strain Energy Density');
        title(['Loading Configuration ', num2str(l)]);
        legend(line_handles, num2str(phi_values', 'Porosity %0.2f'), 'Location', 'best');
        grid on;
    end

     for l = 1:4
        figure;
        hold on;
        line_handles = gobjects(size(phi_values));
        for j = 1:num_phi
            for i = 1:size(p_results, 1)
                line_handles(j) = plot(F_load_values, squeeze(p_results(l, i, j, :)), '-o', 'LineWidth', 1.5);
            end
        end
        hold off;
        xlabel('Load Factor (Pa)');
        ylabel('Pore Preasure (Pa)');
        title(['Loading Configuration ', num2str(l)]);
        legend(line_handles, num2str(phi_values', 'Porosity %0.2f'), 'Location', 'best');
        grid on;
    end
     
end