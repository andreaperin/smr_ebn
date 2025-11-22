clear
close all
clc

tic

N=10^6;
age=50;

%distribution parameters
PGA_min = 0.1; %minimum value of PGA
PGA_max = 2*9.81; %maximum value of PGA
shape_PGA = 1.11; %shape parameter of the frechet distribution of the PGA
scale_PGA = 0.000133; %scale parameter of the frechet distribution of the PGA
x_PGA = linspace(PGA_min,PGA_max,1000);
PGA_PDF = gevpdf(x_PGA, shape_PGA, scale_PGA, 0); %PGA distribution
PGA_PDF = PGA_PDF/sum(PGA_PDF);
%Response time of ACS [s]
t_min_r_ACS = 30; %minimum value of the response time
t_max_r_ACS = 90; %maximum value of the response time
%Response time of EDG [s] and PDP
t_min_r_EDG = 10; %minimum value of the response time
t_max_r_EDG = 180; %maximum value of the response time

%vectors initialization
PGA=zeros(N,1); ACS_rtime=zeros(N,1); 
EDG_rtime=zeros(N,1); %EDG response time 
pdp_rtime=zeros(N,1); %PDP response time
LOCA_time=zeros(N,1); %LOCA failure time
LOOP_time=zeros(N,1); %LOOP failure time
LHS_time=zeros(N,1); %LHS failure time
MSLB_time=zeros(N,1); %MSLB failure time
MSLBH2_time=zeros(N,1); %MSLBH2 failure time
LOOPH2_time=zeros(N,1); %LOOPH2 failure time
ACS_time=zeros(N,1); %ACS failure time
EDG_time=zeros(N,1); %EDG failure time
pdp_time=zeros(N,1); %PDP failure time

%sampling
parfor i=1:N
    %sample the response times
    ACS_rtime(i) = randi(t_max_r_ACS-t_min_r_ACS,1,1)+t_min_r_ACS;
    EDG_rtime(i) = randi(t_max_r_EDG-t_min_r_EDG,1,1)+t_min_r_EDG;
    pdp_rtime(i) = randi(t_max_r_EDG-t_min_r_EDG,1,1)+t_min_r_EDG;
    %sample the PGA
    PGA(i) = randsample(x_PGA, 1, true, PGA_PDF);

    %get the failure times from the age and PGA
    [LOCA_time(i), LOOP_time(i), LHS_time(i), MSLB_time(i), MSLBH2_time(i), LOOPH2_time(i),ACS_time(i), EDG_time(i),pdp_time(i)] = times_function (PGA(i),age);
end

%simulation
[T_W1, ~, ~, ~] = Simulation_model_hydrogen_parallel( ...
    LOCA_time, ...
    LOOP_time, ...
    LHS_time, ...
    MSLB_time, ...
    MSLBH2_time, ...
    LOOPH2_time, ...
    ACS_time, ...
    ACS_rtime, ...
    EDG_time, ...
    EDG_rtime, ...
    pdp_time, ...
    pdp_rtime);

%post-processing
fail=0;
for i=1:N
    if max(T_W1(i,:)) > 1243.9
        fail=fail+1;
    end
end

disp("Failure probability: " + string(fail/N));

toc