function [ACS_1,ACS_2,ACS_3,ACS_4,LOCA1_time,LOCA2_time,LOCA3_time,LOCA4_time,diesel_1,...
    diesel_2,Power,alpha_1,alpha_2,alpha_3,alpha_4,ACS_initial,pipes_initial,EDG_initial,...
    diesel_1_corrosion,diesel_1_earthquake,ACS_1_corrosion,ACS_1_earthquake,CP_corrosion,...
    CP_earthquake,PDP11,PDP12,PDP21,PDP22,PDP31,PDP32,PDP41,PDP42,MSLB1,MSLB2,MSLB3,...
    MSLB4,LOP_earthquake,LOP_hydrogen,MSLB1_hydrogen,MSLB2_hydrogen,MSLB3_hydrogen,MSLB4_hydrogen, MSLB_earthquake,...
    LHS, LHS_explosion,LHS_overcurrent, thermal_failure] = Failure_sim (age,PGA,PGA2,sigma_C_C, ...
    sigma_R_C,sigma_C_ACS,sigma_R_ACS,t_min_EDG,t_max_EDG,t_min_ACS,t_max_ACS, ...
    maint_effectiveness,mu_l_in_ACS,mu_l_in_CP,sigma_l_in_ACS,sigma_l_in_CP,mu_l_in_EDG, ...
    sigma_l_in_EDG,distance,sw_fun,cb_fun,tra_fun,tr_fun,turbine_fun,demand,jointFailure,design_improvements)
%% Auxiliary cooling data
%SMDFR
if design_improvements
    capacity_factor = 1.5747;
else
    capacity_factor = 1;
end
params = default_params(design_improvements);
LHS=0;
LHS_explosion=0;
LHS_overcurrent=0;
% closed-form failure probabilities
p_ACS = acs_failure_prob(age, PGA, params);
p_EDG = edg_failure_prob(age, PGA, params);
p_LOCA1 = loca_prob(age, PGA, params);
p_MSLB1 = mslb_prob(age, PGA2, params);
p_LOOP = loop_prob(PGA);
diesel_1_corrosion=0;
diesel_1_earthquake=0;
ACS_1_corrosion=0;
ACS_1_earthquake=0;
CP_corrosion=0;
CP_earthquake=0;
Cr_m = 8.769*0.375; %[kg/m^2/y]
Cr_m_red = Cr_m/2;
rho = 7750; %kg/m^3
Cr = Cr_m/rho; %[m/y]
Cr_red = Cr_m_red/rho;
Cr_steam = 30*10^(-6)*8.76; %[m/y] from SANDIA "Design Considerations for Concentrating Solar Power Tower Systems Employing Molten Salt"
t = 0.02/2.5; %m
%IMPERFECT MAINTENANCE schedule
maint_interval = [0.5,1.5,3];
%initial length of corrosion defect
t_ACS = 0.015/2.5; %m
if age>0
    corrosion_ACS = normrnd(mu_l_in_ACS,sigma_l_in_ACS,[1000,1]);
else
    corrosion_ACS=zeros(1000,1);
end
for corr=1:floor(age/maint_interval(1))
    corrosion_ACS = corrosion_ACS + Cr_red*maint_interval(1);
    corrosion_ACS = corrosion_ACS*maint_effectiveness;
end
corrosion_ACS = corrosion_ACS+Cr_red*mod(age,maint_interval(1));
mu = mean(corrosion_ACS);
sigma_ACS = std(corrosion_ACS);
mu_f = t_ACS*0.8-mu;
sigma_f = sigma_ACS;
ACS_lambda = 1-normcdf(0,mu_f,sigma_f);

if age>0
    corrosion_ACS = abs(normrnd(mu_l_in_ACS,sigma_l_in_ACS,[4,1]));
    ACS_initial = corrosion_ACS;
    corrosion_pipes = abs(normrnd(mu_l_in_CP,sigma_l_in_CP,[4,1]));
    corrosion_steam_pipes = abs(normrnd(mu_l_in_CP,sigma_l_in_CP,[4,1]));
    pipes_initial=corrosion_pipes;
else
    corrosion_ACS=zeros(4,1);
    corrosion_pipes=zeros(4,1);
    corrosion_steam_pipes=zeros(4,1);
end
corrosion_ACS = corrosion_ACS+corrosion_fun(age,Cr_red,maint_interval(1),1-maint_effectiveness);
corrosion_pipes = corrosion_pipes+corrosion_fun(age,Cr,maint_interval(3),1-maint_effectiveness);
corrosion_steam_pipes = corrosion_steam_pipes+corrosion_fun(age,Cr_steam,maint_interval(3),1-maint_effectiveness);
%ccf combinations
% n = 4;
% binaryCombinations = dec2bin(0:(2^n-1), n);
% combinations = zeros(size(binaryCombinations));
% % Convert each binary string to a numeric vector
% for w = 1:size(binaryCombinations, 1)
%     combinations(w, :) = arrayfun(@(x) str2double(x), binaryCombinations(w, :));
% end
% alpha values from corrosion state (kept for flow scaling)
c_p = min(90,corrosion_ACS/t*100);
alpha_1 = c_p(1)/100;
alpha_2 = c_p(2)/100;
alpha_3 = c_p(3)/100;
alpha_4 = c_p(4)/100;

% ACS failures: independent draws using closed-form probability
ACS_times = t_min_ACS + randi(t_max_ACS - t_min_ACS, [4,1]);
fail_mask = rand(4,1) < p_ACS;
ACS_times(~fail_mask) = 1200;
ACS_1 = ACS_times(1); ACS_2 = ACS_times(2); ACS_3 = ACS_times(3); ACS_4 = ACS_times(4);

%% Other failure data
%LOP
sigma = 1;
mu = log(10) + sigma^2;
% Generating data
x_loca = linspace(0.01, 500, 4000);
lognorm_pdf = lognpdf(x_loca, mu, sigma);
Power = 1200;
if rand() < p_LOOP
    Power = fix(randsample(x_loca, 1, true, lognorm_pdf));
    LOP_earthquake = 1;
else
    LOP_earthquake = 0;
end

% LOCA times (independent draws using closed-form probability)
loca_draws = rand(4,1) < p_LOCA1;
loca_times = fix(randsample(x_loca, 4, true, lognorm_pdf));
loca_times(~loca_draws) = 1200;
LOCA1_time = loca_times(1);
LOCA2_time = loca_times(2);
LOCA3_time = loca_times(3);
LOCA4_time = loca_times(4);

%EDG
EDG_rate = [ones(1,2).*1.22*10^(-3),ones(1,8).*1.22*10^(-3)/4];
corrosion_EDG = abs(normrnd(mu_l_in_EDG,sigma_l_in_EDG,[1,10]));
EDG_initial = corrosion_EDG;
corrosion_EDG = corrosion_EDG+corrosion_fun(age,EDG_rate,maint_interval(2),1-maint_effectiveness);
diesel_1 = 1200; diesel_2 = 1200;
PDP11 = 1200; PDP12 = 1200; PDP21 = 1200; PDP22 = 1200;
PDP31 = 1200; PDP32 = 1200; PDP41 = 1200; PDP42 = 1200;

% independent Bernoulli draws for each EDG component using closed-form probability
diesel_times = t_min_EDG + randi(t_max_EDG - t_min_EDG, [2,1]);
fail_mask = rand(2,1) < p_EDG;
diesel_times(~fail_mask) = 1200;
diesel_1 = diesel_times(1); diesel_2 = diesel_times(2);

pdp_times = t_min_EDG + randi(t_max_EDG - t_min_EDG, [8,1]);
fail_pdp = rand(8,1) < p_EDG;
pdp_times(~fail_pdp) = 1200;
PDP11 = pdp_times(1); PDP12 = pdp_times(2); PDP21 = pdp_times(3); PDP22 = pdp_times(4);
PDP31 = pdp_times(5); PDP32 = pdp_times(6); PDP41 = pdp_times(7); PDP42 = pdp_times(8);

% LOCA and MSLB seismic/corrosion blocks removed in favor of closed-form sampling

fragility_fun = @(PGA,alpha,beta_r,beta_u) normcdf(log(PGA/alpha)/sqrt(beta_r^2+beta_u^2));

%overcurrent
relay = fragility_fun(PGA2/9.81,0.9*capacity_factor,0.35,0.37);
busbar = fragility_fun(PGA2/9.81,1.476*capacity_factor,0.35,0.37);
transformer = tra_fun(PGA2/9.81);
circuit_breaker = cb_fun(PGA2/9.81);
r1 = rand();r2 = rand();r3 = rand();r4 = rand();r5 = rand();r6 = rand();r7 = rand(); r8 = rand(); r9=rand();
if ((r1<transformer && r2>transformer && r3<transformer) || (r4<busbar && r5<busbar && r6<busbar)) && (r7<circuit_breaker && r8<relay && r9<relay)
    LHS = 1;
    LHS_overcurrent = 1;
end

%main steam line break probability
valve_f = fragility_fun(PGA2/9.81,3.8*capacity_factor,0.35,0.5); %from "Separation Requirements for a Hydrogen Production Plant and High-Temperature Nuclear Reactor"
%valve_f = fragility_fun(PGA,22.13,0.27,0.37); %from "NuScale - Chapter Nineteen Probabilistic Risk Assessment" pag. 214
heat_f = fragility_fun(PGA2/9.81,3.65*capacity_factor,0.12,0.51);
%heat_f = fragility_fun(PGA/9.81,6.81,0.32,0.51); %from "NuScale - Chapter Nineteen Probabilistic Risk Assessment" pag. 214

%fault tree from "Expansion of Hazards and Probabilistic Risk Assessments of a Light-Water Reactor Coupled with Electrolysis Hydrogen Production Plants"

%valve = 1 means "fails open" valve == 2 means "ruptures"
%the percentage of fail open and ruptures is taken from the fault tree reference
MSLB1_hydrogen = 0;
MSLB2_hydrogen = 0;
MSLB3_hydrogen = 0;
MSLB4_hydrogen = 0;
MSLB_earthquake = 0;

%valve = 1 means "fails open" valve == 2 means "ruptures"
%the percentage of fail open and ruptures is taken from the fault tree reference
v22 = rand()<valve_f; if v22, v22 = 2*(rand()<0.02) + 1*(rand()>=0.02); end
v23 = rand()<valve_f; if v23, v23 = 2*(rand()<0.02) + 1*(rand()>=0.02); end
v24 = rand()<valve_f; if v24, v24 = 2*(rand()<0.02) + 1*(rand()>=0.02); end
v7 = rand()<valve_f; if v7, v7 = 2*(rand()<0.02) + 1*(rand()>=0.02); end
v2 = rand()<valve_f; if v2, v2 = 2*(rand()<0.02) + 1*(rand()>=0.02); end
v27 = rand()<valve_f; if v27, v27 = 2*(rand()<0.02) + 1*(rand()>=0.02); end
v3 = rand()<valve_f; if v3, v3 = 2*(rand()<0.02) + 1*(rand()>=0.02); end
v4 = rand()<valve_f; if v4, v4 = 2*(rand()<0.02) + 1*(rand()>=0.02); end
h1 = rand()<heat_f;
h2 = rand()<heat_f;

%check for CCF
if v22>0 && rand()<jointFailure(PGA2)
    v7 = v22;
else
    if v7>0 && rand()<jointFailure(PGA2)
        v22 = v7;
    end
end
if v22>0 && rand()<jointFailure(PGA2)
    v27 = v22;
else
    if v27>0 && rand()<jointFailure(PGA2)
        v22 = v27;
    end
end
if v2>0 && rand()<jointFailure(PGA2)
    v7 = v2;
else
    if v7>0 && rand()<jointFailure(PGA2)
        v2 = v7;
    end
end
if v2>0 && rand()<jointFailure(PGA2)
    v27 = v2;
else
    if v27>0 && rand()<jointFailure(PGA2)
        v2 = v27;
    end
end

%if there are reboilers
MSLB1_hydrogen = 0;
MSLB2_hydrogen = 0;
MSLB3_hydrogen = 0;
MSLB4_hydrogen = 0;
MSLB_earthquake = 0;
if  LHS || (v22 && (v23 == 2 || v24 == 2)) || (h1 && (v7 && (v2 || v22))) || (h2 && (v27 && (v2 || v22))) || ((v2 || v22) && (v7 == 2 || v27 == 2)) || (v2 && (v3 == 2 || v4 == 2)) || v2 == 2 || v22 == 2
    MSLB1 = fix(randsample(x_loca, 1, true, lognorm_pdf));
else
    MSLB1=1200;
end

v22 = rand()<valve_f; if v22, v22 = 2*(rand()<0.02) + 1*(rand()>=0.02); end
v23 = rand()<valve_f; if v23, v23 = 2*(rand()<0.02) + 1*(rand()>=0.02); end
v24 = rand()<valve_f; if v24, v24 = 2*(rand()<0.02) + 1*(rand()>=0.02); end
v7 = rand()<valve_f; if v7, v7 = 2*(rand()<0.02) + 1*(rand()>=0.02); end
v2 = rand()<valve_f; if v2, v2 = 2*(rand()<0.02) + 1*(rand()>=0.02); end
v27 = rand()<valve_f; if v27, v27 = 2*(rand()<0.02) + 1*(rand()>=0.02); end
v3 = rand()<valve_f; if v3, v3 = 2*(rand()<0.02) + 1*(rand()>=0.02); end
v4 = rand()<valve_f; if v4, v4 = 2*(rand()<0.02) + 1*(rand()>=0.02); end
h1 = rand()<heat_f;
h2 = rand()<heat_f;

%check for CCF
if v22>0 && rand()<jointFailure(PGA2)
    v7 = v22;
else
    if v7>0 && rand()<jointFailure(PGA2)
        v22 = v7;
    end
end
if v22>0 && rand()<jointFailure(PGA2)
    v27 = v22;
else
    if v27>0 && rand()<jointFailure(PGA2)
        v22 = v27;
    end
end
if v2>0 && rand()<jointFailure(PGA2)
    v7 = v2;
else
    if v7>0 && rand()<jointFailure(PGA2)
        v2 = v7;
    end
end
if v2>0 && rand()<jointFailure(PGA2)
    v27 = v2;
else
    if v27>0 && rand()<jointFailure(PGA2)
        v2 = v27;
    end
end

if LHS || (v22 && (v23 == 2 || v24 == 2)) || (h1 && (v7 && (v2 || v22))) || (h2 && (v27 && (v2 || v22))) || ((v2 || v22) && (v7 == 2 || v27 == 2)) || (v2 && (v3 == 2 || v4 == 2)) || v2 == 2 || v22 == 2
    MSLB2 = fix(randsample(x_loca, 1, true, lognorm_pdf));
else
    MSLB2=1200;
end

v22 = rand()<valve_f; if v22, v22 = 2*(rand()<0.02) + 1*(rand()>=0.02); end
v23 = rand()<valve_f; if v23, v23 = 2*(rand()<0.02) + 1*(rand()>=0.02); end
v24 = rand()<valve_f; if v24, v24 = 2*(rand()<0.02) + 1*(rand()>=0.02); end
v7 = rand()<valve_f; if v7, v7 = 2*(rand()<0.02) + 1*(rand()>=0.02); end
v2 = rand()<valve_f; if v2, v2 = 2*(rand()<0.02) + 1*(rand()>=0.02); end
v27 = rand()<valve_f; if v27, v27 = 2*(rand()<0.02) + 1*(rand()>=0.02); end
v3 = rand()<valve_f; if v3, v3 = 2*(rand()<0.02) + 1*(rand()>=0.02); end
v4 = rand()<valve_f; if v4, v4 = 2*(rand()<0.02) + 1*(rand()>=0.02); end
h1 = rand()<heat_f;
h2 = rand()<heat_f;

%check for CCF
if v22>0 && rand()<jointFailure(PGA2)
    v7 = v22;
else
    if v7>0 && rand()<jointFailure(PGA2)
        v22 = v7;
    end
end
if v22>0 && rand()<jointFailure(PGA2)
    v27 = v22;
else
    if v27>0 && rand()<jointFailure(PGA2)
        v22 = v27;
    end
end
if v2>0 && rand()<jointFailure(PGA2)
    v7 = v2;
else
    if v7>0 && rand()<jointFailure(PGA2)
        v2 = v7;
    end
end
if v2>0 && rand()<jointFailure(PGA2)
    v27 = v2;
else
    if v27>0 && rand()<jointFailure(PGA2)
        v2 = v27;
    end
end

if  LHS || (v22 && (v23 == 2 || v24 == 2)) || (h1 && (v7 && (v2 || v22))) || (h2 && (v27 && (v2 || v22))) || ((v2 || v22) && (v7 == 2 || v27 == 2)) || (v2 && (v3 == 2 || v4 == 2)) || v2 == 2 || v22 == 2
    MSLB3 = fix(randsample(x_loca, 1, true, lognorm_pdf));
else
    MSLB3=1200;
end

v22 = rand()<valve_f; if v22, v22 = 2*(rand()<0.02) + 1*(rand()>=0.02); end
v23 = rand()<valve_f; if v23, v23 = 2*(rand()<0.02) + 1*(rand()>=0.02); end
v24 = rand()<valve_f; if v24, v24 = 2*(rand()<0.02) + 1*(rand()>=0.02); end
v7 = rand()<valve_f; if v7, v7 = 2*(rand()<0.02) + 1*(rand()>=0.02); end
v2 = rand()<valve_f; if v2, v2 = 2*(rand()<0.02) + 1*(rand()>=0.02); end
v27 = rand()<valve_f; if v27, v27 = 2*(rand()<0.02) + 1*(rand()>=0.02); end
v3 = rand()<valve_f; if v3, v3 = 2*(rand()<0.02) + 1*(rand()>=0.02); end
v4 = rand()<valve_f; if v4, v4 = 2*(rand()<0.02) + 1*(rand()>=0.02); end
h1 = rand()<heat_f;
h2 = rand()<heat_f;

%check for CCF
if v22>0 && rand()<jointFailure(PGA2)
    v7 = v22;
else
    if v7>0 && rand()<jointFailure(PGA2)
        v22 = v7;
    end
end
if v22>0 && rand()<jointFailure(PGA2)
    v27 = v22;
else
    if v27>0 && rand()<jointFailure(PGA2)
        v22 = v27;
    end
end
if v2>0 && rand()<jointFailure(PGA2)
    v7 = v2;
else
    if v7>0 && rand()<jointFailure(PGA2)
        v2 = v7;
    end
end
if v2>0 && rand()<jointFailure(PGA2)
    v27 = v2;
else
    if v27>0 && rand()<jointFailure(PGA2)
        v2 = v27;
    end
end

if  LHS || (v22 && (v23 == 2 || v24 == 2)) || (h1 && (v7 && (v2 || v22))) || (h2 && (v27 && (v2 || v22))) || ((v2 || v22) && (v7 == 2 || v27 == 2)) || (v2 && (v3 == 2 || v4 == 2)) || v2 == 2 || v22 == 2
    MSLB4 = fix(randsample(x_loca, 1, true, lognorm_pdf));
else
    MSLB4=1200;
end

% overwrite with closed-form MSLB probabilities (independent per train)
mslb_times = fix(randsample(x_loca, 4, true, lognorm_pdf));
fail_mslb = rand(4,1) < p_MSLB1;
mslb_times(~fail_mslb) = 1200;
MSLB1 = min(MSLB1,mslb_times(1)); MSLB2 = min(MSLB2,mslb_times(2)); MSLB3 = min(MSLB3,mslb_times(3)); MSLB4 = min(MSLB4,mslb_times(4));
if any(fail_mslb)
    MSLB_earthquake = 1;
end

%thermal stress failure in HTE
thermal_failure = 0;
pipe_failure_prob = 0.646;
if (MSLB1<1200 || MSLB2<1200 || MSLB3<1200 || MSLB4<1200) && rand()<heat_f
    if rand()<pipe_failure_prob
        thermal_failure = 1;
        thermal_failure_time = 60;
    end
end

%hydrogen release probability

%p_isolate = 0.9; %from HyRAM ESD
%alternative approach from "PERFORMANCE ASSESSMENT OF SAFETY BARRIERS IN LIQUID HYDROGEN BUNKERING OPERATIONS USING BAYESIAN NETWORK"
%all barriers have to fail to have an explosion
p_26 = 2.4e-2; p_27 = 1.17e-1; p_28 = 2.5e-1; p_29 = 4.77e-2; p_30 = 2.11E-2;
p_31 = 3.36e-1; p_32 = 4.42e-3; p_33 = 1.42e-3; p_34 = 2.52e-3; p_35 = 1.7e-3;
p_36 = 2.36e-3; p_37 = 1.54e-4; p_7 = 9.82e-3;
p_ESD = (1-(1-p_26)*(1-p_27)*(1-p_28))*(1-(1-p_29)*(1-p_30)*(1-p_31));
p_alarm_activation = (1-(1-p_34)*(1-(p_32+p_33-p_32*p_33)));
p_alarm = (1-(1-p_alarm_activation)*(1-(p_35+p_36-p_35*p_36)));
p_detection = p_alarm*(1-(1-p_29)*(1-p_7)*(1-p_37));
p_isolate = 1-(p_ESD+p_detection-p_ESD*p_detection); %Release Detection Barrier
hydrogen_barriers_age = mod(age,1.5)*12;
RDP_deg = (3.19e-6*hydrogen_barriers_age^2-5.34e-5*hydrogen_barriers_age+3.89e-4)/4.27e-4;
p_isolate = p_isolate*RDP_deg;

%piping_f = fragility_fun(PGA2,4.8,0.35,0.5);
piping_f = fragility_fun(PGA2,2.8*capacity_factor,0.35,0.5);
LOP=0;
demand_p = demand(randi(8760));
LOP_thermal = 0;
LHS_thermal = 0;
if (rand()<piping_f && rand()<(1-p_isolate)) || thermal_failure  
    t = randi(120*60);
    d = 0.001+rand()*0.66;
    [LOP,m_dot_release,~,~,overpressure] = Loss_of_Power(d,t,sw_fun,tr_fun,demand_p,distance);
    %from HyRAM ESD
    EPB_deg = (3.72e-5*hydrogen_barriers_age+4.5817e-4)/5.89e-4;
    IPB_deg = (1.51e-6*hydrogen_barriers_age^2-2.1e-5*hydrogen_barriers_age+7.19e-5)/1.42e-4;
    if m_dot_release < 0.125
        p_delayed_ign = 0.004; %Ignition Prevention Barrier + Escalation Prevention Barrier
    else
        if m_dot_release < 6.25
            p_delayed_ign = 0.027; %Ignition Prevention Barrier + Escalation Prevention Barrier
        else
            p_delayed_ign = 0.120; %Ignition Prevention Barrier + Escalation Prevention Barrier
        end
    end
    p_delayed_ign = p_delayed_ign*IPB_deg*EPB_deg;
    if rand()>p_delayed_ign
        LOP = 0;
    else
        if rand()<turbine_fun(overpressure)
            LHS = 1;
            LHS_explosion = 1;

            if thermal_failure
                LHS_thermal=1;
            end
        end

        if LOP
            if thermal_failure
                LOP_thermal=1;
            end
        end
    end
end

if LOP
    if LOP_thermal
        Power = thermal_failure_time;
    else
        Power = fix(randsample(x_loca, 1, true, lognorm_pdf));
    end
    LOP_hydrogen = 1;
else
    LOP_hydrogen = 0;
end

if LHS
    if LHS_thermal  
         MSLB1 = thermal_failure_time;
         MSLB2 = thermal_failure_time;
         MSLB3 = thermal_failure_time;
         MSLB4 = thermal_failure_time;
    else
        MSLB1 = fix(randsample(x_loca, 1, true, lognorm_pdf));
        MSLB2 = fix(randsample(x_loca, 1, true, lognorm_pdf));
        MSLB3 = fix(randsample(x_loca, 1, true, lognorm_pdf));
        MSLB4 = fix(randsample(x_loca, 1, true, lognorm_pdf));
    end
end

if LHS == 0
    LHS = 1200;
else
    LHS = MSLB1;
end

if MSLB1_hydrogen == 0
    MSLB1_hydrogen=1200;
else
    MSLB1_hydrogen = MSBL1;
end

end
