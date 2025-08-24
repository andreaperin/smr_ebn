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
LHS=0;
LHS_explosion=0;
LHS_overcurrent=0;
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
for corr=1:floor(age/maint_interval(1))
    corrosion_ACS = corrosion_ACS + Cr_red*maint_interval(1);
    corrosion_ACS = max(ACS_initial,corrosion_ACS*maint_effectiveness);
end
for corr=1:floor(age/maint_interval(3))
    corrosion_pipes = corrosion_pipes + Cr*maint_interval(3);
    corrosion_pipes = max(pipes_initial,corrosion_pipes*maint_effectiveness);
    corrosion_steam_pipes = corrosion_steam_pipes + Cr_steam*maint_interval(3);
    corrosion_steam_pipes = max(pipes_initial,corrosion_steam_pipes*maint_effectiveness);
end
corrosion_ACS = corrosion_ACS+Cr_red*mod(age,maint_interval(1));
corrosion_pipes = corrosion_pipes+Cr*mod(age,maint_interval(3));
corrosion_steam_pipes = corrosion_steam_pipes+Cr_steam*mod(age,maint_interval(3));
%ccf combinations
% n = 4;
% binaryCombinations = dec2bin(0:(2^n-1), n);
% combinations = zeros(size(binaryCombinations));
% % Convert each binary string to a numeric vector
% for w = 1:size(binaryCombinations, 1)
%     combinations(w, :) = arrayfun(@(x) str2double(x), binaryCombinations(w, :));
% end
combinations = [0	0	0	0;
0	0	0	1;
0	0	1	0;
0	0	1	1;
0	1	0	0;
0	1	0	1;
0	1	1	0;
0	1	1	1;
1	0	0	0;
1	0	0	1;
1	0	1	0;
1	0	1	1;
1	1	0	0;
1	1	0	1;
1	1	1	0;
1	1	1	1];
ACS_comb = (1+ACS_lambda)/2;
ccf_check = 0;
for w=1:size(combinations,1)-1
    P_comb = ACS_lambda*ACS_comb^(sum(1-combinations(w,:))-1)*(1-ACS_comb)^sum(combinations(w,:));
    if rand()<P_comb
        if combinations(w,1) == 1
            ACS_1 = 1200;
        else
            ACS_1 = t_min_ACS+randi(t_max_ACS-t_min_ACS);
        end
        if combinations(w,2) == 1
            ACS_2 = 1200;
        else
            ACS_2 = t_min_ACS+randi(t_max_ACS-t_min_ACS);
        end
        if combinations(w,3) == 1
            ACS_3 = 1200;
        else
            ACS_3 = t_min_ACS+randi(t_max_ACS-t_min_ACS);
        end
        if combinations(w,4) == 1
            ACS_4 = 1200;
        else
            ACS_4 = t_min_ACS+randi(t_max_ACS-t_min_ACS);
        end
        ccf_check = 1;
        break
    end
end

%PGA_fit = [1,1.5,2,2.5,3].*9.81;
%R = [523.05,655.63,806.04,950.15,1060];
C = 687*1.2;
c_p = min(90,corrosion_ACS/t*100);
alpha_1 = c_p(1)/100;
alpha_2 = c_p(2)/100;
alpha_3 = c_p(3)/100;
alpha_4 = c_p(4)/100;
C = C*((-0.009611)*c_p+1)*capacity_factor;
%R_f = fit(PGA_fit',R','poly2');
R_f = (-0.1537).*PGA.^2+33.9285.*PGA+199.8460;

F = 1-normcdf((log(C)-log(R_f))/sqrt(sigma_C_ACS^2+sigma_R_ACS^2));

%seismic failure with CCF
if ccf_check == 0
    if corrosion_ACS(1)>t_ACS*0.8 || rand()<F(1)
        ACS_1 = t_min_ACS+randi(t_max_ACS-t_min_ACS);
        if Pb_acs(1)<P_acs(1)
            ACS_1_corrosion=1;
        else
            ACS_1_earthquake=1;
        end
    else
        ACS_1 = 1200;
    end
    if corrosion_ACS(2)>t_ACS*0.8 || rand()<F(2)
        ACS_2 = t_min_ACS+randi(t_max_ACS-t_min_ACS);
    else
        ACS_2 = 1200;
    end
    if corrosion_ACS(3)>t_ACS*0.8 || rand()<F(3)
        ACS_3 = t_min_ACS+randi(t_max_ACS-t_min_ACS);
    else
        ACS_3 = 1200;
    end
    if corrosion_ACS(4)>t_ACS*0.8 || rand()<F(4)
        ACS_4 = t_min_ACS+randi(t_max_ACS-t_min_ACS);
    else
        ACS_4 = 1200;
    end
end

%% Other failure data
%LOP
sigma = 1;
mu = log(10) + sigma^2;
% Generating data
x_loca = linspace(0.01, 500, 4000);
lognorm_pdf = lognpdf(x_loca, mu, sigma);
F_power = normcdf((log(PGA/(9.81*0.2))+0.25*normcdf(0.95))/0.2);
Power = 1200;
if rand() < F_power
    Power = fix(randsample(x_loca, 1, true, lognorm_pdf));
    LOP_earthquake = 1;
else
    LOP_earthquake = 0;
end

%EDG
%EDG_rate = [1.22*10^(-3),1.22*10^(-3),1.22*10^(-3)/2,1.22*10^(-3)/2];
EDG_rate = [ones(1,2).*1.22*10^(-3),ones(1,8).*1.22*10^(-3)/4];
corrosion_EDG = abs(normrnd(mu_l_in_EDG,sigma_l_in_EDG,[1,10]));
EDG_initial = corrosion_EDG;
%IMPERFECT MAINTENANCE
for corr=1:floor(age/maint_interval(2))
    corrosion_EDG = corrosion_EDG+EDG_rate.*maint_interval(2);
    corrosion_EDG = max(EDG_initial,corrosion_EDG*maint_effectiveness);
end
corrosion_EDG = corrosion_EDG+EDG_rate.*mod(age,maint_interval(2));
t_diesel = [0.015/2.5, 0.045/2.5]; %m
diesel_1 = 1200;
diesel_2 = 1200;
PDP11 = 1200;
PDP12 = 1200;
PDP21 = 1200;
PDP22 = 1200;
PDP31 = 1200;
PDP32 = 1200;
PDP41 = 1200;
PDP42 = 1200;
if corrosion_EDG(1) > t_diesel(1)*0.8
    diesel_1 = t_min_EDG+randi(t_max_EDG-t_min_EDG);
    diesel_1_corrosion=1;
    %check ccf
    if rand() < 1.19e-2
        diesel_2 = t_min_EDG+randi(t_max_EDG-t_min_EDG);
    end
end
if corrosion_EDG(2) > t_diesel(1)*0.8
    diesel_2 = t_min_EDG+randi(t_max_EDG-t_min_EDG);
    %check ccf
    if rand() < 1.19e-2
        diesel_1 = t_min_EDG+randi(t_max_EDG-t_min_EDG);
    end
end
if corrosion_EDG(3) > t_diesel(2)*0.8
    PDP11 = t_min_EDG+randi(t_max_EDG-t_min_EDG);
    %check ccf
    if rand() < 1.19e-2
        PDP12 = t_min_EDG+randi(t_max_EDG-t_min_EDG);
    end
end
if corrosion_EDG(4) > t_diesel(2)*0.8
    PDP12 = t_min_EDG+randi(t_max_EDG-t_min_EDG);
    %check ccf
    if rand() < 1.19e-2
        PDP11 = t_min_EDG+randi(t_max_EDG-t_min_EDG);
    end
end
if corrosion_EDG(5) > t_diesel(2)*0.8
    PDP21 = t_min_EDG+randi(t_max_EDG-t_min_EDG);
    %check ccf
    if rand() < 1.19e-2
        PDP22 = t_min_EDG+randi(t_max_EDG-t_min_EDG);
    end
end
if corrosion_EDG(6) > t_diesel(2)*0.8
    PDP22 = t_min_EDG+randi(t_max_EDG-t_min_EDG);
    %check ccf
    if rand() < 1.19e-2
        PDP21 = t_min_EDG+randi(t_max_EDG-t_min_EDG);
    end
end
if corrosion_EDG(7) > t_diesel(2)*0.8
    PDP31 = t_min_EDG+randi(t_max_EDG-t_min_EDG);
    %check ccf
    if rand() < 1.19e-2
        PDP32 = t_min_EDG+randi(t_max_EDG-t_min_EDG);
    end
end
if corrosion_EDG(8) > t_diesel(2)*0.8
    PDP32 = t_min_EDG+randi(t_max_EDG-t_min_EDG);
    %check ccf
    if rand() < 1.19e-2
        PDP31 = t_min_EDG+randi(t_max_EDG-t_min_EDG);
    end
end
if corrosion_EDG(9) > t_diesel(2)*0.8
    PDP41 = t_min_EDG+randi(t_max_EDG-t_min_EDG);
    %check ccf
    if rand() < 1.19e-2
        PDP42 = t_min_EDG+randi(t_max_EDG-t_min_EDG);
    end
end
if corrosion_EDG(10) > t_diesel(2)*0.8
    PDP42 = t_min_EDG+randi(t_max_EDG-t_min_EDG);
    %check ccf
    if rand() < 1.19e-2
        PDP41 = t_min_EDG+randi(t_max_EDG-t_min_EDG);
    end
end

A = [0.95*9.81;0.95*9.81;ones(8,1).*3.35*9.81];
c_p = min(80,corrosion_EDG'/t_diesel(1)*100);
A = A.*((-0.009611).*c_p+1).*capacity_factor;
bu = 0.26;
br = 0.24;
Q = 0.95;

F = normcdf((log([PGA,PGA,ones(1,8).*PGA2]./A')+bu*normcdf(Q))/br);
%F=0;
rho12 = 0.3158; %obtained with gaussian copula
if rand()<F(1)
    diesel_1 = t_min_EDG+randi(t_max_EDG-t_min_EDG);
    diesel_1_earthquake=1;
    if rand()<rho12
        diesel_2 = t_min_EDG+randi(t_max_EDG-t_min_EDG);
    end
end
if rand()<F(2) 
    diesel_2 = t_min_EDG+randi(t_max_EDG-t_min_EDG);
    if rand()<rho12
        diesel_1 = t_min_EDG+randi(t_max_EDG-t_min_EDG);
    end
end
if rand()<F(3)
    if rand()<F(3)
        PDP11 = t_min_EDG+randi(t_max_EDG-t_min_EDG);
        if rand()<rho12
            PDP12 = t_min_EDG+randi(t_max_EDG-t_min_EDG);
        end
    end
end
if rand()<F(4)
    if rand()<F(4)
        PDP12 = t_min_EDG+randi(t_max_EDG-t_min_EDG);
        if rand()<rho12
            PDP11 = t_min_EDG+randi(t_max_EDG-t_min_EDG);
        end
    end
end

if rand()<F(5)
    if rand()<F(5)
        PDP21 = t_min_EDG+randi(t_max_EDG-t_min_EDG);
        if rand()<rho12
            PDP22 = t_min_EDG+randi(t_max_EDG-t_min_EDG);
        end
    end
end
if rand()<F(6)
    if rand()<F(6)
    PDP22 = t_min_EDG+randi(t_max_EDG-t_min_EDG);
    if rand()<rho12
        PDP21 = t_min_EDG+randi(t_max_EDG-t_min_EDG);
    end
    end
end
if rand()<F(7)
    if rand()<F(7)
    PDP31 = t_min_EDG+randi(t_max_EDG-t_min_EDG);
    if rand()<rho12
        PDP32 = t_min_EDG+randi(t_max_EDG-t_min_EDG);
    end
    end
end
if rand()<F(8)
    if rand()<F(8)
    PDP32 = t_min_EDG+randi(t_max_EDG-t_min_EDG);
    if rand()<rho12
        PDP31 = t_min_EDG+randi(t_max_EDG-t_min_EDG);
    end
    end
end
if rand()<F(9)
    if rand()<F(9)
    PDP41 = t_min_EDG+randi(t_max_EDG-t_min_EDG);
    if rand()<rho12
        PDP42 = t_min_EDG+randi(t_max_EDG-t_min_EDG);
    end
    end
end
if rand()<F(10)
    if rand()<F(10)
    PDP42 = t_min_EDG+randi(t_max_EDG-t_min_EDG);
    if rand()<rho12
        PDP41 = t_min_EDG+randi(t_max_EDG-t_min_EDG);
    end
    end
end

LOCA1 = 0;
LOCA2 = 0;
LOCA3 = 0;
LOCA4 = 0;
t_pipes=t_ACS;
if corrosion_pipes(1)>t_pipes*0.8
    CP_corrosion=1;
    LOCA1 = 1;
end
if corrosion_pipes(2)>t_pipes*0.8
    LOCA2 = 1;
end
if corrosion_pipes(3)>t_pipes*0.8
    LOCA3 = 1;
end
if corrosion_pipes(4)>t_pipes*0.8
    LOCA4 = 1;
end

%PGA_fit = [1,1.5,2,2.5,3].*9.81;
%R = [523.05,655.63,806.04,950.15,1060];
C = 1337;
c_p = min(90,corrosion_pipes/t*100);
C = C*((-0.009611)*c_p+1)*capacity_factor;
R_f = (-0.1537).*PGA.^2+33.9285.*PGA+199.8460;
%R_f = fit(PGA_fit',R','poly2'); fit with data from [Seismic Fragility Evaluation of Main Steam Piping of Isolated APR1400 NPP Considering the Actual Failure Mode]
%%%%%%%%%%%%%%MAYBE TO BE CHANGED TO THE DATA FROM [Seismic fragility
%%%%%%%%%%%%%%analysis of seismically isolated nuclear power plantspiping system]
F = 1-normcdf((log(C)-log(R_f))/sqrt(sigma_C_C^2+sigma_R_C^2));

if rand()<F(1)
    CP_earthquake=1;
    LOCA1 = 1;
end
if rand()<F(2)
    LOCA2 = 1;
end
if rand()<F(3)
    LOCA3 = 1;
end
if rand()<F(4)
    LOCA4 = 1;
end

C = 1337;
c_p = min(90,corrosion_steam_pipes/t*100);
C = C*((-0.009611)*c_p+1)*capacity_factor;
R_f = (-0.1537).*PGA2.^2+33.9285.*PGA2+199.8460;
%R_f = fit(PGA_fit',R','poly2'); fit with data from [Seismic Fragility Evaluation of Main Steam Piping of Isolated APR1400 NPP Considering the Actual Failure Mode]
%%%%%%%%%%%%%%MAYBE TO BE CHANGED TO THE DATA FROM [Seismic fragility
%%%%%%%%%%%%%%analysis of seismically isolated nuclear power plantspiping system]
F = 1-normcdf((log(C)-log(R_f))/sqrt(sigma_C_C^2+sigma_R_C^2));

MSLB1_check = 0;
MSLB2_check = 0;
MSLB3_check = 0;
MSLB4_check = 0;

if rand()<F(1) || corrosion_steam_pipes(1)>t_pipes*0.8
    MSLB1_check = 1;
end
if rand()<F(2) || corrosion_steam_pipes(2)>t_pipes*0.8
    MSLB2_check = 1;
end
if rand()<F(3) || corrosion_steam_pipes(3)>t_pipes*0.8
    MSLB3_check = 1;
end
if rand()<F(4) || corrosion_steam_pipes(4)>t_pipes*0.8
    MSLB4_check = 1;
end

if LOCA1
    LOCA1_time = fix(randsample(x_loca, 1, true, lognorm_pdf));
else
    LOCA1_time = 1200;
end
if LOCA2
    LOCA2_time = fix(randsample(x_loca, 1, true, lognorm_pdf));
else
    LOCA2_time = 1200;
end
if LOCA3
    LOCA3_time = fix(randsample(x_loca, 1, true, lognorm_pdf));
else
    LOCA3_time = 1200;
end
if LOCA4
    LOCA4_time = fix(randsample(x_loca, 1, true, lognorm_pdf));
else
    LOCA4_time = 1200;
end

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
if MSLB1_check || LHS || (v22 && (v23 == 2 || v24 == 2)) || (h1 && (v7 && (v2 || v22))) || (h2 && (v27 && (v2 || v22))) || ((v2 || v22) && (v7 == 2 || v27 == 2)) || (v2 && (v3 == 2 || v4 == 2)) || v2 == 2 || v22 == 2
    MSLB1 = fix(randsample(x_loca, 1, true, lognorm_pdf));
    if MSLB1_check == 1
        MSLB_earthquake = 1;
    else
        if LHS == 0
            MSLB1_hydrogen = 1;
        end
    end
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

if MSLB2_check || LHS || (v22 && (v23 == 2 || v24 == 2)) || (h1 && (v7 && (v2 || v22))) || (h2 && (v27 && (v2 || v22))) || ((v2 || v22) && (v7 == 2 || v27 == 2)) || (v2 && (v3 == 2 || v4 == 2)) || v2 == 2 || v22 == 2
    MSLB2 = fix(randsample(x_loca, 1, true, lognorm_pdf));
    if MSLB2_check == 1
        MSLB_earthquake = 1;
    else
        if LHS == 0
            MSLB2_hydrogen = 1;
        end
    end
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

if MSLB3_check || LHS || (v22 && (v23 == 2 || v24 == 2)) || (h1 && (v7 && (v2 || v22))) || (h2 && (v27 && (v2 || v22))) || ((v2 || v22) && (v7 == 2 || v27 == 2)) || (v2 && (v3 == 2 || v4 == 2)) || v2 == 2 || v22 == 2
    MSLB3 = fix(randsample(x_loca, 1, true, lognorm_pdf));
    if MSLB3_check == 1
        MSLB_earthquake = 1;
    else
        if LHS == 0
            MSLB3_hydrogen = 1;
        end
    end
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

if MSLB4_check || LHS || (v22 && (v23 == 2 || v24 == 2)) || (h1 && (v7 && (v2 || v22))) || (h2 && (v27 && (v2 || v22))) || ((v2 || v22) && (v7 == 2 || v27 == 2)) || (v2 && (v3 == 2 || v4 == 2)) || v2 == 2 || v22 == 2
    MSLB4 = fix(randsample(x_loca, 1, true, lognorm_pdf));
    if MSLB4_check == 1
        MSLB_earthquake = 1;
    else
        if LHS == 0
            MSLB4_hydrogen = 1;
        end
    end
else
    MSLB4=1200;
end

%thermal stress failure in HTE
thermal_failure = 0;
pipe_failure_prob = 0.646;
if (MSLB1<1200 || MSLB2<1200 || MSLB3<1200 || MSLB4<1200) && rand()<heat_f
    for i=1:length(pipe_failure_prob)
        if rand()<pipe_failure_prob(i)
            thermal_failure = 1;
            thermal_failure_time = 60*i;
        end
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