% Homogenization from EXTRAVASCULAR MATRIX to MACRO BONE
function [Chom_macro]=hom_exvas_to_macro(phi_vas) 
% Edited: IMWS Pircher 2023-01
% Input: phi_vas ... vascular porosity = V_vas / V_macro
% Output: Chom_macro ... homogenized stiffness matrix of macro bone
% material in [GPa]
% Chom_exvas (taken from Blanchard et al 2016) =
% homogenized stiffness matrix of extravascular matrix in [GPa] 


%% 1.0 Specification of Parameters
%% 1.1 General - Tensordefinition
% 4th-order identity tensor
I=[1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1];
% volumetric part of the 4th-order identity tensor
J=[1/3 1/3 1/3 0 0 0; 1/3 1/3 1/3 0 0 0; 1/3 1/3 1/3 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
% deviatoric part of the 4th-order identity tensor
K=I-J; 
% second-order unit tensor
one2=eye(3);

%% 1.2 Material Properties General

%%%%%%%%%% STIFFNESS PARAMETERS START %%%%%%%%%%
% Stiffness Water
K_H2O = 2.3; %GPa
mu_H2O = 0.0; %GPa

% Stiffness vas = vascular pores (water filled pores)
K_vas = K_H2O; %GPa
mu_vas = mu_H2O; %GPa
Cmat_vas = 3*K_vas*J + 2*mu_vas*K ;%;%GPa

% Stiffness of extravascular matrix (taken from Blanchard et al 2016)
Chom_exvas = [ 10.56 4.82 5.30 0 0 0;
    4.82 10.56 5.3 0 0 0;
    5.3 5.3 13.06 0 0 0;
    0 0 0 5.99 0 0;
    0 0 0 0 5.99 0;
    0 0 0 0 0 5.75]; %GPa
%%%%%%%%%% STIFFNESS PARAMETERS END %%%%%%%%%%

%% 1.3 Volume Fractions
% RVE macro - Macro bone material
f_vas_macro = phi_vas;
f_exvas_macro = 1-f_vas_macro;

%% 2.0 RVE MACRO - Macro bone material
% Macro bone material - Mori Tanaka Scheme
% Phase 1 = Matrix: exvas
% Phase 2 = vas = cylindrical inclusions in transversal isotropic matrix

% Morphology Tensors = Hill Tensors
P_cyl_in_Chom_exvas = P_isotrans_cyl(Chom_exvas);

% Infinite Strain Concentration Tensors
A_inf_exvas = I;
A_inf_vas = inv(I + P_cyl_in_Chom_exvas*(Cmat_vas - Chom_exvas));

% Actual Strain Concentration Tensors
A_exvas = A_inf_exvas* inv( f_exvas_macro * A_inf_exvas + f_vas_macro * A_inf_vas );
A_vas = A_inf_vas* inv( f_exvas_macro * A_inf_exvas + f_vas_macro * A_inf_vas );

% Estimated Homogenized Stiffness
Chom_macro = f_exvas_macro * Chom_exvas * A_exvas + f_vas_macro * Cmat_vas * A_vas;
