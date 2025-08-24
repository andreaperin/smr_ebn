clear
close all
clc

total_iter = 1;
start = 1;
N = 10^3; %number of montecarlo trials
static = 0;
HTE_dynamic = 1;
design_improvements = 0;
saving=1;

if static
    age_to_sim = 0.1;
else
    age_to_sim = 59.9;
end

for iter=start:total_iter

    tic

    load('Fragility_interp.mat','sw_fun','turbine_fun','tr_fun')
    load('US_hourly_electricity_demand','demand')
    load('Joint_Failure_Functon.mat','jointFailure')
    load('Seismic_fragility_electrical.mat','cb_fun','transformer_fun')

    warning('off','all')

    maint_effectiveness = 0.1;
    save_out = 0;
    model = 'SMDFR';
    All_Outputs = 1;

    %% INPUTS
    tsim = 1200; %simulation time
    distance = 250;

    %% DISTRIBUTIONS PARAMETERS
    %PGA [m/s^2]
    PGA_min = 0.1; %minimum value of PGA
    PGA_max = 2*9.81; %maximum value of PGA
    shape_PGA = 1.11; %shape parameter of the frechet distribution of the PGA
    scale_PGA = 0.000133; %scale parameter of the frechet distribution of the PGA

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

    %Response time of ACS [s]
    t_min_r_ACS = 30; %minimum value of the response time
    t_max_r_ACS = 90; %maximum value of the response time

    %Response time of EDG [s]
    t_min_r_EDG = 10; %minimum value of the response time
    t_max_r_EDG = 180; %maximum value of the response time

    %% INPUTS PRE-PROCESSING
    x_PGA = linspace(PGA_min,PGA_max,1000);
    PGA_PDF = gevpdf(x_PGA, shape_PGA, scale_PGA, 0); %PGA distribution
    PGA_PDF = PGA_PDF/sum(PGA_PDF);
    file = 'SMDFR_Parameters';
    evalc(file);
    age=zeros(N,1);
    PGA=zeros(N,1);
    PGA2=zeros(N,1);
    ACS1_response_time=zeros(N,1);
    ACS2_response_time=zeros(N,1);
    ACS3_response_time=zeros(N,1);
    ACS4_response_time=zeros(N,1);
    EDG1_response_time=zeros(N,1);
    EDG2_response_time=zeros(N,1);
    PDP11_response_time=zeros(N,1);
    PDP12_response_time=zeros(N,1);
    PDP21_response_time=zeros(N,1);
    PDP22_response_time=zeros(N,1);
    PDP31_response_time=zeros(N,1);
    PDP32_response_time=zeros(N,1);
    PDP41_response_time=zeros(N,1);
    PDP42_response_time=zeros(N,1);
    ACS_1=zeros(N,1);
    ACS_2=zeros(N,1);
    ACS_3=zeros(N,1);
    ACS_4=zeros(N,1);
    PDP11=zeros(N,1);
    PDP12=zeros(N,1);
    PDP21=zeros(N,1);
    PDP22=zeros(N,1);
    PDP31=zeros(N,1);
    PDP32=zeros(N,1);
    PDP41=zeros(N,1);
    PDP42=zeros(N,1);
    MSLB1=zeros(N,1);
    MSLB2=zeros(N,1);
    MSLB3=zeros(N,1);
    MSLB4=zeros(N,1);
    LOCA1_time=zeros(N,1);
    LOCA2_time=zeros(N,1);
    LOCA3_time=zeros(N,1);
    LOCA4_time=zeros(N,1);
    EDG_1=zeros(N,1);
    EDG_2=zeros(N,1);
    Power=zeros(N,1);
    alpha_1=zeros(N,1);
    alpha_2=zeros(N,1);
    alpha_3=zeros(N,1);
    alpha_4=zeros(N,1);
    ACS_initial=zeros(N,4);
    pipes_initial=zeros(N,4);
    EDG_initial=zeros(N,10);
    EDG_1_corrosion = zeros(N,1);
    EDG_1_earthquake=zeros(N,1);
    ACS_1_corrosion=zeros(N,1);
    ACS_1_earthquake=zeros(N,1);
    CP_corrosion=zeros(N,1);
    CP_earthquake=zeros(N,1);
    LOP_earthquake=zeros(N,1);
    LOP_hydrogen=zeros(N,1);
    MSLB1_hydrogen=zeros(N,1);
    MSLB2_hydrogen=zeros(N,1);
    MSLB3_hydrogen=zeros(N,1);
    MSLB4_hydrogen=zeros(N,1);
    MSLB_earthquake=zeros(N,1);
    LHS = zeros(N,1);
    LHS_explosion = zeros(N,1);
    LHS_overcurrent = zeros(N,1);
    thermal_failure = zeros(N,1);
    thermal_failure_time = zeros(N,1);

    parfor i=1:N
        %uncertain variables sampling
        age(i) = age_to_sim;
        PGA(i) = randsample(x_PGA, 1, true, PGA_PDF);
        PGA2(i) = exp(log(PGA(i))+normrnd(0,5e-5*distance+0.36,1,1));
        PGA2(i) = PGA(i);
        ACS1_response_time(i) = randi(t_max_r_ACS-t_min_r_ACS)+t_min_r_ACS;
        ACS2_response_time(i) = randi(t_max_r_ACS-t_min_r_ACS)+t_min_r_ACS;
        ACS3_response_time(i) = randi(t_max_r_ACS-t_min_r_ACS)+t_min_r_ACS;
        ACS4_response_time(i) = randi(t_max_r_ACS-t_min_r_ACS)+t_min_r_ACS;
        EDG1_response_time(i) = randi(t_max_r_EDG-t_min_r_EDG)+t_min_r_EDG;
        EDG2_response_time(i) = randi(t_max_r_EDG-t_min_r_EDG)+t_min_r_EDG;
        PDP11_response_time(i) = randi(t_max_r_EDG-t_min_r_EDG)+t_min_r_EDG;
        PDP12_response_time(i) = randi(t_max_r_EDG-t_min_r_EDG)+t_min_r_EDG;
        PDP21_response_time(i) = randi(t_max_r_EDG-t_min_r_EDG)+t_min_r_EDG;
        PDP22_response_time(i) = randi(t_max_r_EDG-t_min_r_EDG)+t_min_r_EDG;
        PDP31_response_time(i) = randi(t_max_r_EDG-t_min_r_EDG)+t_min_r_EDG;
        PDP32_response_time(i) = randi(t_max_r_EDG-t_min_r_EDG)+t_min_r_EDG;
        PDP41_response_time(i) = randi(t_max_r_EDG-t_min_r_EDG)+t_min_r_EDG;
        PDP42_response_time(i) = randi(t_max_r_EDG-t_min_r_EDG)+t_min_r_EDG;
        if static
            [ACS_1(i),ACS_2(i),ACS_3(i),ACS_4(i),LOCA1_time(i),LOCA2_time(i),LOCA3_time(i),...
                LOCA4_time(i),EDG_1(i),EDG_2(i),Power(i),alpha_1(i),alpha_2(i),alpha_3(i),alpha_4(i),...
                ACS_initial(i,:),pipes_initial(i,:),EDG_initial(i,:),EDG_1_corrosion(i),EDG_1_earthquake(i),...
                ACS_1_corrosion(i),ACS_1_earthquake(i),CP_corrosion(i),CP_earthquake(i),PDP11(i),PDP12(i),...
                PDP21(i),PDP22(i),PDP31(i),PDP32(i),PDP41(i),PDP42(i),MSLB1(i),MSLB2(i),MSLB3(i),MSLB4(i),...
                LOP_earthquake(i),LOP_hydrogen(i),MSLB1_hydrogen(i),MSLB2_hydrogen(i),MSLB3_hydrogen(i),...
                MSLB4_hydrogen(i), MSLB_earthquake(i),LHS(i),LHS_explosion(i),LHS_overcurrent(i),...
                thermal_failure(i)] = Failure_sim_static (age(i),PGA(i),PGA2(i),sigma_C_C,sigma_R_C,sigma_C_ACS,sigma_R_ACS, ...
                t_min_EDG,t_max_EDG,t_min_ACS,t_max_ACS, ...
                maint_effectiveness,mu_l_in_ACS,mu_l_in_CP,sigma_l_in_ACS,sigma_l_in_CP,mu_l_in_EDG, ...
                sigma_l_in_EDG,distance,sw_fun,cb_fun,transformer_fun,tr_fun,turbine_fun,demand,jointFailure);
        else
            if HTE_dynamic
                model = 'SMDFR_HTE_model';
                [ACS_1(i),ACS_2(i),ACS_3(i),ACS_4(i),LOCA1_time(i),LOCA2_time(i),LOCA3_time(i),...
                    LOCA4_time(i),EDG_1(i),EDG_2(i),Power(i),alpha_1(i),alpha_2(i),alpha_3(i),alpha_4(i),...
                    ACS_initial(i,:),pipes_initial(i,:),EDG_initial(i,:),EDG_1_corrosion(i),EDG_1_earthquake(i),...
                    ACS_1_corrosion(i),ACS_1_earthquake(i),CP_corrosion(i),CP_earthquake(i),PDP11(i),PDP12(i),...
                    PDP21(i),PDP22(i),PDP31(i),PDP32(i),PDP41(i),PDP42(i),MSLB1(i),MSLB2(i),MSLB3(i),MSLB4(i),...
                    LOP_earthquake(i),LOP_hydrogen(i),MSLB1_hydrogen(i),MSLB2_hydrogen(i),MSLB3_hydrogen(i),...
                    MSLB4_hydrogen(i), MSLB_earthquake(i),LHS(i),LHS_explosion(i),LHS_overcurrent(i),...
                    thermal_failure(i), thermal_failure_time(i)] = Failure_sim_HTE_model (age(i),PGA(i),PGA2(i),sigma_C_C,sigma_R_C,sigma_C_ACS,sigma_R_ACS, ...
                    t_min_EDG,t_max_EDG,t_min_ACS,t_max_ACS, ...
                    maint_effectiveness,mu_l_in_ACS,mu_l_in_CP,sigma_l_in_ACS,sigma_l_in_CP,mu_l_in_EDG, ...
                    sigma_l_in_EDG,distance,sw_fun,cb_fun,transformer_fun,tr_fun,turbine_fun,demand,jointFailure);
            else   
                [ACS_1(i),ACS_2(i),ACS_3(i),ACS_4(i),LOCA1_time(i),LOCA2_time(i),LOCA3_time(i),...
                    LOCA4_time(i),EDG_1(i),EDG_2(i),Power(i),alpha_1(i),alpha_2(i),alpha_3(i),alpha_4(i),...
                    ACS_initial(i,:),pipes_initial(i,:),EDG_initial(i,:),EDG_1_corrosion(i),EDG_1_earthquake(i),...
                    ACS_1_corrosion(i),ACS_1_earthquake(i),CP_corrosion(i),CP_earthquake(i),PDP11(i),PDP12(i),...
                    PDP21(i),PDP22(i),PDP31(i),PDP32(i),PDP41(i),PDP42(i),MSLB1(i),MSLB2(i),MSLB3(i),MSLB4(i),...
                    LOP_earthquake(i),LOP_hydrogen(i),MSLB1_hydrogen(i),MSLB2_hydrogen(i),MSLB3_hydrogen(i),...
                    MSLB4_hydrogen(i), MSLB_earthquake(i),LHS(i),LHS_explosion(i),LHS_overcurrent(i),...
                    thermal_failure(i)] = Failure_sim (age(i),PGA(i),PGA2(i),sigma_C_C,sigma_R_C,sigma_C_ACS,sigma_R_ACS, ...
                    t_min_EDG,t_max_EDG,t_min_ACS,t_max_ACS, ...
                    maint_effectiveness,mu_l_in_ACS,mu_l_in_CP,sigma_l_in_ACS,sigma_l_in_CP,mu_l_in_EDG, ...
                    sigma_l_in_EDG,distance,sw_fun,cb_fun,transformer_fun,tr_fun,turbine_fun,demand,jointFailure,design_improvements);
            end
        end
        %injection of inputs
        ACS1_flow = 25;
        ACS2_flow = 25;
        ACS3_flow = 25;
        ACS4_flow = 25;
        if (LOCA1_time(i) < 1200) || (LOCA2_time(i) < 1200) || (LOCA3_time(i) < 1200) || (LOCA4_time(i) < 1200) || (MSLB1(i) < 1200) || (MSLB2(i) < 1200) || (MSLB3(i) < 1200) || (MSLB4(i) < 1200)
            in(i) = Simulink.SimulationInput(model);
            in(i) = in(i).setModelParameter('StartTime','1','StopTime',num2str(tsim));
            in(i) = in(i).setVariable('All_Outputs',All_Outputs);
            in(i) = in(i).setVariable('LOCA1_time',LOCA1_time(i));
            in(i) = in(i).setVariable('LOCA2_time',LOCA2_time(i));
            in(i) = in(i).setVariable('LOCA3_time',LOCA3_time(i));
            in(i) = in(i).setVariable('LOCA4_time',LOCA4_time(i));
            in(i) = in(i).setVariable('ACS_1',ACS_1(i));
            in(i) = in(i).setVariable('ACS_2',ACS_2(i));
            in(i) = in(i).setVariable('ACS_3',ACS_3(i));
            in(i) = in(i).setVariable('ACS_4',ACS_4(i));
            in(i) = in(i).setVariable('EDG_1',EDG_1(i));
            in(i) = in(i).setVariable('EDG_2',EDG_2(i));
            in(i) = in(i).setVariable('Power',Power(i));
            in(i) = in(i).setVariable('ACS1_response_time',ACS1_response_time(i));
            in(i) = in(i).setVariable('ACS2_response_time',ACS2_response_time(i));
            in(i) = in(i).setVariable('ACS3_response_time',ACS3_response_time(i));
            in(i) = in(i).setVariable('ACS4_response_time',ACS4_response_time(i));
            in(i) = in(i).setVariable('EDG1_response_time',EDG1_response_time(i));
            in(i) = in(i).setVariable('EDG2_response_time',EDG2_response_time(i));
            in(i) = in(i).setVariable('PDP11_response_time',PDP11_response_time(i));
            in(i) = in(i).setVariable('PDP12_response_time',PDP12_response_time(i));
            in(i) = in(i).setVariable('PDP21_response_time',PDP21_response_time(i));
            in(i) = in(i).setVariable('PDP22_response_time',PDP22_response_time(i));
            in(i) = in(i).setVariable('PDP31_response_time',PDP31_response_time(i));
            in(i) = in(i).setVariable('PDP32_response_time',PDP32_response_time(i));
            in(i) = in(i).setVariable('PDP41_response_time',PDP41_response_time(i));
            in(i) = in(i).setVariable('PDP42_response_time',PDP42_response_time(i));
            in(i) = in(i).setVariable('PDP11',PDP11(i));
            in(i) = in(i).setVariable('PDP12',PDP12(i));
            in(i) = in(i).setVariable('PDP21',PDP21(i));
            in(i) = in(i).setVariable('PDP22',PDP22(i));
            in(i) = in(i).setVariable('PDP31',PDP31(i));
            in(i) = in(i).setVariable('PDP32',PDP32(i));
            in(i) = in(i).setVariable('PDP41',PDP41(i));
            in(i) = in(i).setVariable('PDP42',PDP42(i));
            in(i) = in(i).setVariable('MSLB1',MSLB1(i));
            in(i) = in(i).setVariable('MSLB2',MSLB2(i));
            in(i) = in(i).setVariable('MSLB3',MSLB3(i));
            in(i) = in(i).setVariable('MSLB4',MSLB4(i));
            in(i) = in(i).setVariable('ACS1_flow',ACS1_flow*alpha_1(i));
            in(i) = in(i).setVariable('ACS2_flow',ACS2_flow*alpha_2(i));
            in(i) = in(i).setVariable('ACS3_flow',ACS3_flow*alpha_3(i));
            in(i) = in(i).setVariable('ACS4_flow',ACS4_flow*alpha_4(i));

            if HTE_dynamic
                in(i) = in(i).setVariable('thermal_failure',thermal_failure(i));
                in(i) = in(i).setVariable('LHS',LHS(i));
                in(i) = in(i).setVariable('thermal_failure_time',thermal_failure_time(i));
                in(i) = in(i).setVariable('PGA',PGA(i));
            end
        end
        % if mod(i,dN) == 0
        %     send(D,[]);
        % end

    end
    %delete(w);
    if exist("in","var")
        indexes = arrayfun(@(x) ~isempty(x.ModelName), in);
        in_filtered = in(indexes(1:length(in)));
    else
        in_filtered = [];
        indexes = [];
    end
    clear in
    toc
    safe = N-length(in_filtered);
    age = age(indexes);
    PGA = PGA(indexes);
    PGA2 = PGA2(indexes);
    ACS1_response_time = ACS1_response_time(indexes);
    ACS2_response_time = ACS2_response_time(indexes);
    ACS3_response_time = ACS3_response_time(indexes);
    ACS4_response_time = ACS4_response_time(indexes);
    EDG1_response_time = EDG1_response_time(indexes);
    EDG2_response_time = EDG2_response_time(indexes);
    ACS_1 = ACS_1(indexes);
    ACS_2 = ACS_2(indexes);
    ACS_3 = ACS_3(indexes);
    ACS_4 = ACS_4(indexes);
    LOCA1_time = LOCA1_time(indexes);
    LOCA2_time = LOCA2_time(indexes);
    LOCA3_time = LOCA3_time(indexes);
    LOCA4_time = LOCA4_time(indexes);
    MSLB1 = MSLB1(indexes);
    MSLB2 = MSLB2(indexes);
    MSLB3 = MSLB3(indexes);
    MSLB4 = MSLB4(indexes);
    PDP11 = PDP11(indexes);
    PDP12 = PDP12(indexes);
    PDP21 = PDP21(indexes);
    PDP22 = PDP22(indexes);
    PDP31 = PDP31(indexes);
    PDP32 = PDP32(indexes);
    PDP41 = PDP41(indexes);
    PDP42 = PDP42(indexes);
    PDP11_response_time = PDP11_response_time(indexes);
    PDP12_response_time = PDP12_response_time(indexes);
    PDP21_response_time = PDP21_response_time(indexes);
    PDP22_response_time = PDP22_response_time(indexes);
    PDP31_response_time = PDP31_response_time(indexes);
    PDP32_response_time = PDP32_response_time(indexes);
    PDP41_response_time = PDP41_response_time(indexes);
    PDP42_response_time = PDP42_response_time(indexes);
    EDG_1 = EDG_1(indexes);
    EDG_2 = EDG_2(indexes);
    Power = Power(indexes);
    alpha_1 = alpha_1(indexes);
    alpha_2 = alpha_2(indexes);
    alpha_3 = alpha_3(indexes);
    alpha_4 = alpha_4(indexes);
    ACS_initial = ACS_initial(indexes,:);
    pipes_initial = pipes_initial(indexes,:);
    EDG_initial = EDG_initial(indexes,:);
    LOP_earthquake = LOP_earthquake(indexes);
    LOP_hydrogen = LOP_hydrogen(indexes);
    ACS_1_corrosion = ACS_1_corrosion(indexes);
    ACS_1_earthquake = ACS_1_earthquake(indexes);
    CP_corrosion = CP_corrosion(indexes);
    CP_earthquake = CP_earthquake(indexes);
    EDG_1_corrosion = EDG_1_corrosion(indexes);
    EDG_1_earthquake = EDG_1_earthquake(indexes);
    LHS = LHS(indexes);
    MSLB1_hydrogen = MSLB1_hydrogen(indexes);
    MSLB2_hydrogen = MSLB2_hydrogen(indexes);
    MSLB3_hydrogen = MSLB3_hydrogen(indexes);
    MSLB4_hydrogen = MSLB4_hydrogen(indexes);

    %count_MSLB = length(find(MSLB1<1200 | MSLB2<1200 | MSLB3<1200 | MSLB4<1200));
    LOCA_percentage = length(find(LOCA1_time<1200|LOCA2_time<1200|LOCA3_time<1200|LOCA4_time<1200))/N
    MSLB_hydrogen_percentage = (sum(MSLB1_hydrogen)+sum(MSLB2_hydrogen)+sum(MSLB3_hydrogen)+sum(MSLB4_hydrogen))/N
    MSLB_earthquake_percentage = sum(MSLB_earthquake)/N
    LOP_hydrogen_percentage = sum(LOP_hydrogen)/N
    LOP_earthquake_percentage = sum(LOP_earthquake)/N
    LOP_hydrogen_percentage_increase = sum(LOP_hydrogen)/sum(LOP_earthquake)
    LHS_percentage = sum(LHS)/N
    LHS_explosion_percentage = sum(LHS_explosion)/N
    LHS_overcurrent_percentage = sum(LHS_overcurrent)/N
    thermal_failure_percentage = sum(thermal_failure)/N

    date = string(datetime('now', 'Format', 'dd-MMM-yyyy_HH-mm-ss'));
    %save(strcat('All_Inputs_',date,'.mat'))
    
    if saving
        save(fullfile(pwd, strcat('Input_', num2str(iter))),'-v7.3');
    end

    disp('Number of scenarios: ')
    disp(length(in_filtered))

    toc
end