function params = default_params(design_improvements)
%DEFAULT_PARAMS Builds default parameter structs for closed-form calculators.
%   params = DEFAULT_PARAMS(design_improvements) returns a struct with fields
%   ACS, EDG, LOCA, each containing the constants used by Failure_sim.
%
%   design_improvements = 0 (default) uses capacity_factor = 1.
%   design_improvements = 1 uses capacity_factor = 1.5747 (as in Failure_sim).

    if nargin<1
        design_improvements = 0;
    end
    capacity_factor = design_improvements*1.5747 + (1-design_improvements)*1;
    params.ACS = default_acs_params(capacity_factor);
    params.EDG = default_edg_params(capacity_factor);
    params.LOCA = default_loca_params(capacity_factor);
    params.MSLB = default_mslb_params(capacity_factor);
end

function p = default_acs_params(capacity_factor)
    p.capacity_factor = capacity_factor;
    p.Cr_m = 8.769*0.375; rho = 7750;
    p.Cr = p.Cr_m/rho;
    p.Cr_red = p.Cr/50;
    p.t = 0.02/2.5;
    p.t_ACS = 0.015;
    p.dt = 0.5; % ACS maintenance interval years
    p.maint_effectiveness = 0.01; % same as Only_Inputs default (1-effectiveness)
    p.mu0 = 0.00159/10;             % initial corrosion mean
    p.sigma0 = 0.000619/10;         % initial corrosion std
    p.sigma_C = 0.14;
    p.sigma_R = 0.07;
    p.C0 = 687;
    p.R_poly = [-0.1537, 33.9285, 199.8460];
end

function p = default_edg_params(capacity_factor)
    p.capacity_factor = capacity_factor;
    p.rate = 1.22e-3/50;     % m/y (approx for diesel_1 path)
    p.dt = 1.5;           % maintenance interval years (maint_interval(2))
    p.maint_effectiveness = 0.01;
    p.mu0 = 0.00159;      % initial corrosion mean
    p.sigma0 = 0.000619;  % initial corrosion std
    p.t_diesel = 0.015;
    p.A0 = 0.95*9.81*2;     % diesel_1 seismic capacity base
    p.bu = 0.26;
    p.br = 0.24;
    p.Q = 0.95;
end

function p = default_loca_params(capacity_factor)
    p.capacity_factor = capacity_factor;
    p.Cr_m = 8.769*0.375; rho = 7750;
    p.Cr = p.Cr_m/rho/60;     % m/y
    p.t = 0.02/2.5;
    p.t_pipe = 0.015; % same thickness used in Failure_sim for coolant pipes
    p.dt = 3;             % maint_interval(3) in Failure_sim for pipes
    p.maint_effectiveness = 0.01;
    p.mu0 = 0.00159;      % initial corrosion mean
    p.sigma0 = 0.000619;  % initial corrosion std
    p.sigma_C = 0.14;
    p.sigma_R = 0.07;
    p.C0 = 1337;
    p.R_poly = [-0.1537, 33.9285, 199.8460];
end

function p = default_mslb_params(capacity_factor)
    p.capacity_factor = capacity_factor;
    p.Cr_steam = 30e-6*8.76/60; % m/y (scaled similar to other rates here)
    p.t = 0.02/2.5;
    p.t_pipe = 0.015;     % steam pipe thickness used in Failure_sim logic
    p.dt = 3;             % maintenance interval (maint_interval(3))
    p.maint_effectiveness = 0.01;
    p.mu0 = 0.00159;
    p.sigma0 = 0.000619;
    p.sigma_C = 0.14;
    p.sigma_R = 0.07;
    p.C0 = 1337;
    p.R_poly = [-0.1537, 33.9285, 199.8460];
end
