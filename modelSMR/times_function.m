function [LOCA_time, LOOP_time, LHS_time, MSLB_time, MSLBH2_time, LOOPH2_time,ACS_time, EDG_time,pdp_time] = times_function (PGA,age)

load('Fragility_interp.mat','sw_fun','turbine_fun','tr_fun')
load('US_hourly_electricity_demand','demand')
load('Joint_Failure_Functon.mat','jointFailure')
load('Seismic_fragility_electrical.mat','cb_fun','transformer_fun')

maint_effectiveness = 0.01; %it is actually 1-effectiveness

%% INPUTS
tsim = 1200; %simulation time
distance = 250;

%Coolant pipes corrosion defect [m]
mu_l_in_CP = 0.00159; %mean of the initial depth of the corrosion defect
sigma_l_in_CP = 0.000619; %standard deviation of the initial depth of the corrosion defect

%ACS corrosion defect [m]
mu_l_in_ACS = 0.00159; %mean of the initial depth of the corrosion defect
sigma_l_in_ACS = 0.000619; %standard deviation of the initial depth of the corrosion defect

%Seismic capacity of coolant pipes
sigma_C_C = 0.14; %standard deviation of the capacity

%Seismic structural response of coolant pipes
sigma_R_C = 0.07; %standard deviation of the structural response

%Seismic capacity of ACS
sigma_C_ACS=0.14; %standard deviation of the capacity

%Seismic structural response of ACS
sigma_R_ACS = 0.07; %standard deviation of the structural response

%Corrosion depth of EDG cooling pipes [m]
mu_l_in_EDG = 0.00159; %mean of the initial depth of the corrosion defect
sigma_l_in_EDG = 0.000619; %standard deviation of the initial depth of the corrosion defect

%Failure time of EDG [s]
t_min_EDG = 15; %minimum value of the failure time
t_max_EDG = tsim; %maximum value of the failure time

%Failure time of ACS [s]
t_min_ACS = 15; %minimum value of the failure time
t_max_ACS = tsim; %maximum value of the failure time

PGA2=PGA;
design_improvements=0;

[ACS_time,~,~,~,LOCA_time,~,~,...
    ~,EDG_time,~,LOOP_time,~,~,~,~,...
    ~,~,~,~,~,...
    ~,~,~,~,pdp_time,~,...
    ~,~,~,~,~,~,MSLB_time,~,~,~,...
    ~,LOOPH2_time,MSLBH2_time,~,~,...
    ~, ~,LHS_time,~,~,...
    ~] = Failure_sim (age,PGA,PGA2,sigma_C_C,sigma_R_C,sigma_C_ACS,sigma_R_ACS, ...
    t_min_EDG,t_max_EDG,t_min_ACS,t_max_ACS, ...
    maint_effectiveness,mu_l_in_ACS,mu_l_in_CP,sigma_l_in_ACS,sigma_l_in_CP,mu_l_in_EDG, ...
    sigma_l_in_EDG,distance,sw_fun,cb_fun,transformer_fun,tr_fun,turbine_fun,demand,jointFailure,design_improvements);

end