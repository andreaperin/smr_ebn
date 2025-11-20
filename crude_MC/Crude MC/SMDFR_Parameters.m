%% Steam generator data
%STEAM ENTHALPY 4.4*10^6
Tc_in = 973;
Tin_c = Tc_in;
d_c = (10.678-13.174*10^(-4)*(Tc_in-600.6))*1000; %coolant density [kg/m^3]
k_sg = 500; %trial value [W/mK]
A_sg = 370; %trial value [m^2]
C_c = 140.2; % [J/kgK]
G_s = 27000; %trial value [kg/s]
V_s = 100; %trial value [m^3]
d_st = 6.123; %trial value [kg/m^3]
T_fw = 300; %trail value [°C]
T_st = 900; %trial value [°C]
h_st_in = 4.38645*10^6; % trial value [J/K]
h_fw_st = h_st_in; % trial value [J/K]
V_sg = 100; %trial value [m^3]

%%Rankine data
eta=0.4;
alpha = 0.9926; %percent of steam that goes to rankine cycle

%% Cu-Cl cycle data
V_hst = 100; %trial value [m^3]
d_hst = 3.14; %trial value [kg/m^3]
C_hst = 1996000; % trial value [J/kgK]
h_hpo = 7159680; % trial value [kJ/K]
K_h = k_sg;
A_h = A_sg/1.7;
T_hpi = 100; %trial value [°C]
A_h20 = 18.015; %[kg/kmol]
A_hcl = 36.458; %[kg/kmol]
A_cucl = 134.446; %[kg/kmol]
A_cuocl = 213.991; %[kg/kmol]
p_h20 = 40; %trial value [bar]
C_h20 = 4.48; %trial value [kJ/kgK]
C_hcl = 0.8; %trial value [kJ/kgK]
C_cucl = 0.71; %trial value [kJ/kgK]
C_cuocl = 1; %trial value [kJ/kgK]
G_h20i = 10; %v_h20 stable at 3!!!
G_h20o = 7;
G_hcl = 12.1426; %v_hcl = 12.142
G_cucl = 44.778; %v_cucl = 44.778
G_cuocl = 35.6355; %v_cuocl = 35.6355
h_in = 7*10^3; % trial value [kJ/K]
h_out = 10^3; % trial value [kJ/K]
P_hydro = 130+37; % trial value [kW]
k_h = 5*10^(-3);
E_h = 5; % trial value [J/mol]
V_h = 160; %trial value [m^3]
V_dst = V_hst;
d_dst = d_hst;
C_dst = C_hst;
h_dpo = 6.8010e+06;
K_d = K_h;
A_d = A_h;
A_o = 32;
C_o = 0.918;
G_o = 2.66;
P_decomposition = 160;
k_d = 2.7187;
E_d = E_h;
V_d = V_h;
k_e = k_d;
E_e = E_d;
M_hcl0 = G_hcl;
M_cucl0 = G_cucl;
V_h2 = V_d;
A_h2 = A_d;
M_cucl20 = G_cucl;
M_h20 = 8;
G_h2 = 32;
A_cucl2 = A_cucl;
G_cucl2 = G_cucl;
k_hte = k_d;
E_hte = E_d;
V_c = V_d;
G_h20 = G_h20i;


%% Neutronics data for the point kinetic model

Lambda = 1.05E-6;
beta_i = [7.91E-5, 7.03E-4, 5.04E-4, 1.17E-3, 4.57E-4, 1.1E-4];
beta = 3.02E-3;
lambda_i = [1.27E-2, 3E-2, 1.1E-1, 3.19E-1, 1.18, 7.02];
phi_tot = 3.61E21;
v_rec = 4.16E-9;
N_0 = 1.5E13;
Pi_0 = [1.58E-1, sum([6.42E-2, 6.98E-2, 7.1E-2, 7.14E-2, 7.12E-2, 7.07E-2, 7E-2, 6.91E-2, 6.76E-2, 6.49E-2]), 1.52E-1];
P_0 = 1.E8;
tau_c = 5.43;
tau_e = 10;
lambda_c = 0.184;
lambda_e = 0.1;
alfa_f = -2.08E-4;
alfa_c = -8.28E-6;
rho_0 = 1.11E-3; %reactivity, in cents, reported in Liu et al 2022
%rho_0 = 204.2E-6; % Considering rho_0 reported in Wang, 2017
C_i0 = [4.29E14, 4.15E15, 3.83E15, 1.16E16, 5.7E15, 1.53E15]';
Ce_i0 = [7.01E14, 5.88E15, 3.36E15, 5.09E15, 8.20E14, 3.96E13]';

%% Parameters
L = 6;          %delayed neutronic groups
Nodes = 3;     %nodes for energy balance

%% Thermodynamic parameters - low resolution model
fuelpipes = 1027;
coolpipes = [2166, 2166,0];
A_fw = [13.6, 90, 13.6];
A_wc = [10.9, 103, 10.9];
mdot_f = [327.7, 327.7, 327.7];
mdot_c = [5550.5, 5550.5, 5550.5];
M_f = [307.7, 1164, 307.7];
M_w = [39.3, 311, 39.3];
M_c = [219.9, 7705, 219.9];
c_pf = [400, 400, 400];
c_pw = [690, 690, 690];
c_pc = [140.2, 140.2, 140.2];
Heat = [18744000, 5970000, 16956000];
h_fw = [6546.0, 4305.1, 7629.5];
h_wc = [27405.4, 13833.4, 32177.3];

%% Assumptions (vact. Sep 27. 2023)
%0. Total reactivity formula: modified 

%1. Tin_x and Tout_x reported in Liu et al. (2021)

%Defined by design
Tin_f = 1300;
Tout_f = 1300;

Tout_c = 1100;

%2. Tbar_x: key for stability at first seconds of simulation

%They do not have to be equal to Tout_x!

%Defined by inspection (in order to P_N(t=tsim) = 1)
Tbar_f0 = 1288;
Tbar_c0 = 1060;

% HYPERPARAMETERS - SIMULATION
%nsim: taken from the HRM
nsim = 1;
step = 0.1;

% BIAS AND DRIFT - INITIAL VALUES FOR THE SIMULATION
mdot_f_value = mdot_f(1);    
mdot_c_value = mdot_c(1);
leakage = 5;
% external reactivity
rhoex_value = 0;    

ACS1_flow = 15; %[kg/s]
ACS2_flow = 15; %[kg/s]
ACS3_flow = 15; %[kg/s]
ACS4_flow = 15; %[kg/s]

LOPC_time = 10^5;
EPB = 1;
RPB = 1;
IPB = 1;
DPB = 1;