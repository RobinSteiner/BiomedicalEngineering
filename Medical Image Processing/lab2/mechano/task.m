% task for 202.674 Mechanobiology - Fundamentals and Modeling Concepts
% Hannah Katharina Fabro, 11719424

clear all
close all
clc

makeplots = 1;

% -------------------- defintion of variables -----------------------
E = 25; %E modulus in GPa
nue = 0.28; %Poisson's ratio
K_matrix = E/(3*(1-2*nue)); %bulk modulus
G_matrix = E/(2*(1+nue)); %shear modulus
C_matrix = 1/((1+nue)*(1-2*nue))*...
    [(1-nue)*E nue*E nue*E 0 0 0;
    nue*E (1-nue)*E nue*E 0 0 0;
    nue*E nue*E (1-nue)*E 0 0 0;
    0 0 0 E*(1-2*nue) 0 0;
    0 0 0 0 E*(1-2*nue) 0;
    0 0 0 0 0 E*(1-2*nue)];

k_H2O = 2.3; %GPa
K = [1/3 1/3 1/3 0 0 0;
    1/3 1/3 1/3 0 0 0;
    1/3 1/3 1/3 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0];
C_pore = 3*k_H2O*K;

% -------------------------- task 1 ---------------------------
% compute strain energy density for all loading scenarios
F_load = linspace(0,10,50); % load factor in MPa
phi = linspace(0.25,0.75,50); %porosity

% ----------------------- calculations task 1 ---------------------------

for i = 1:4 %amount of stress tensors
    for j = 1:length(phi)
        for k = 1:length(F_load)

            % stress tensors in Kelvin-Mandel-Notation
            if i == 1
                Sigma_RVE = [-2*F_load(k) 0 0 0 0 0]';
            end
            if i == 2
                Sigma_RVE = [-F_load(k) -F_load(k) -F_load(k) 0 0 0]';
            end
            if i == 3
                Sigma_RVE = [0 0 0 0 0 sqrt(2)*F_load(k)]';
            end
            if i == 4
                Sigma_RVE = [-3*F_load(k) -F_load(k) -1.5*F_load(k) sqrt(2)*1.8*F_load(k) sqrt(2)*0.1*F_load(k) sqrt(2)*0.7*F_load(k)]';
            end

            P_isosphere = inv(C_matrix).*(3*K_matrix*K/(3*K_matrix+4*G_matrix)+6*(K_matrix+2*G_matrix)*(eye(6,6)-K)/(15*K_matrix+20*G_matrix));
            A_pore = inv(eye(6,6)+P_isosphere.*(C_pore-C_matrix)).*(inv(phi(j)*(inv(eye(6,6)+P_isosphere.*(C_pore-C_matrix)))));
            A_matrix = eye(6,6).*(inv(phi(j)*(inv(eye(6,6)+P_isosphere.*(C_pore-C_matrix)))+(1-phi(j))*eye(6,6)));

            %homogenized stiffness tensor
            C_RVE = phi(j)*C_pore*A_pore + (1-phi(j))*C_matrix*A_matrix;

            Epsilon_RVE = inv(C_RVE)*Sigma_RVE;
            epsilon_matrix = A_matrix*Epsilon_RVE;
            psi_matrix(j,k) = 1/2*epsilon_matrix'*C_matrix*epsilon_matrix; %strain energy density
        end
    end
    if i == 1
        psi_matrix_Sigma_RVE1 = psi_matrix;
    end
    if i == 2
        psi_matrix_Sigma_RVE2 = psi_matrix;
    end
    if i == 3
        psi_matrix_Sigma_RVE3 = psi_matrix;
    end
    if i == 4
        psi_matrix_Sigma_RVE4 = psi_matrix;
    end
end

clear Epsilon_RVE i j k

% -------------------------- task 2 ---------------------------
% compute hydrostatic pressure for all loading scenarios

%get b_pore in Kelvin-Mandel-Notation by making \fat{1} tensor 6x1
allone = [1 1 1 0 0 0]';

phi_0 = 0.25; %initial porosity of the RVE; i.e., the volume fraction of
% the pore phase before applying any mechanical loading

for i = 1:4 %amount of stress tensors
    for j = 1:length(phi)
        for k = 1:length(F_load)

            % stress tensors in Kelvin-Mandel-Notation
            if i == 1
                Sigma_RVE = [-2*F_load(k) 0 0 0 0 0]';
            end
            if i == 2
                Sigma_RVE = [-F_load(k) -F_load(k) -F_load(k) 0 0 0]';
            end
            if i == 3
                Sigma_RVE = [0 0 0 0 0 sqrt(2)*F_load(k)]';
            end
            if i == 4
                Sigma_RVE = [-3*F_load(k) -F_load(k) -1.5*F_load(k) sqrt(2)*1.8*F_load(k) sqrt(2)*0.1*F_load(k) sqrt(2)*0.7*F_load(k)]';
            end

            A_drainedpore = inv(eye(6,6)+P_isosphere.*(zeros(6,6)-C_matrix)).*(inv(phi(j)*(inv(eye(6,6)+P_isosphere.*(zeros(6,6)-C_matrix))))); %C_pore = 0
            A_drainedmatrix = eye(6,6).*(inv(phi(j).*(inv(eye(6,6)+P_isosphere.*(zeros(6,6)-C_matrix)))+(1-phi(j))*eye(6,6))); %C_pore = 0
            b_pore = (phi_0*allone'*A_drainedpore)';

            C_drainedRVE = (1-phi(j))*C_matrix*A_drainedmatrix;
            C_undrainedRVE = C_drainedRVE+inv(phi_0/k_H2O+(allone'*inv(C_matrix)*(b_pore-phi_0*allone)))*(b_pore*b_pore');

            Epsilon_RVE = inv(C_undrainedRVE)*Sigma_RVE;

            p(j,k) = -inv(phi_0/k_H2O + (allone'*inv(C_matrix)*(b_pore-phi_0*allone)))*b_pore'*Epsilon_RVE;
        end
    end
    if i == 1
        p_Sigma_RVE1 = p;
    end
    if i == 2
        p_Sigma_RVE2 = p;
    end
    if i == 3
        p_Sigma_RVE3 = p;
    end
    if i == 4
        p_Sigma_RVE4 = p;
    end
end

% clear A_drainedmatrix A_drainedpore A_matrix A_pore allone b_pore C_drainedRVE
% clear C_matrix C_pore C_RVE C_undrainedRVE E epsilon_matrix Epsilon_RVE G_matrix
% clear P_isosphere Sigma_RVE

% ------------------------------ plots ------------------------------

if makeplots == 1

    % ----------------------- figures for task 1 ------------------------
    subplot(2,2,1)
    surf(phi,F_load,psi_matrix_Sigma_RVE1)
    xlabel('phi')
    ylabel('F_{load}')
    zlabel('psi')
    set(gca,'FontSize',20,'LineWidth',2)
    grid on
    title('Strain energy density of stress tensor \Sigma_{RVE1}')
    subplot(2,2,2)
    surf(phi,F_load,psi_matrix_Sigma_RVE2)
    xlabel('phi')
    ylabel('F_{load}')
    zlabel('psi')
    set(gca,'FontSize',20,'LineWidth',2)
    grid on
    title('Strain energy density of stress tensor \Sigma_{RVE2}')
    subplot(2,2,3)
    surf(phi,F_load,psi_matrix_Sigma_RVE3)
    xlabel('phi')
    ylabel('F_{load}')
    zlabel('psi')
    set(gca,'FontSize',20,'LineWidth',2)
    grid on
    title('Strain energy density of stress tensor \Sigma_{RVE3}')
    subplot(2,2,4)
    surf(phi,F_load,psi_matrix_Sigma_RVE4)
    xlabel('phi')
    ylabel('F_{load}')
    zlabel('psi')
    set(gca,'FontSize',20,'LineWidth',2)
    grid on
    title('Strain energy density of stress tensor \Sigma_{RVE4}')

    figure
    surf(phi,F_load,psi_matrix_Sigma_RVE1)
    hold on
    surf(phi,F_load,psi_matrix_Sigma_RVE2)
    surf(phi,F_load,psi_matrix_Sigma_RVE3)
    surf(phi,F_load,psi_matrix_Sigma_RVE4)
    set(gca,'zscale','log')

    % ----------------------- figures for task 1 ------------------------
    subplot(2,2,1)
    surf(phi,F_load,p_Sigma_RVE1)
    xlabel('phi')
    ylabel('F_{load}')
    zlabel('p')
    set(gca,'FontSize',20,'LineWidth',2)
    grid on
    title('Pore pressure of stress tensor \Sigma_{RVE1}')
    subplot(2,2,2)
    surf(phi,F_load,p_Sigma_RVE2)
    xlabel('phi')
    ylabel('F_{load}')
    zlabel('p')
    set(gca,'FontSize',20,'LineWidth',2)
    grid on
    title('Pore pressure of stress tensor \Sigma_{RVE2}')
    subplot(2,2,3)
    surf(phi,F_load,p_Sigma_RVE3)
    xlabel('phi')
    ylabel('F_{load}')
    zlabel('p')
    set(gca,'FontSize',20,'LineWidth',2)
    grid on
    title('Pore pressure of stress tensor \Sigma_{RVE3}')
    subplot(2,2,4)
    surf(phi,F_load,p_Sigma_RVE4)
    xlabel('phi')
    ylabel('F_{load}')
    zlabel('p')
    set(gca,'FontSize',20,'LineWidth',2)
    grid on
    title('Pore pressure of stress tensor \Sigma_{RVE4}')
end
