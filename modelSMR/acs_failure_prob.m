function p_fail = acs_failure_prob(age, PGA, params)
%ACS_FAILURE_PROB Deterministic approximation of ACS_1 failure probability.
%   p_fail = ACS_FAILURE_PROB(age, PGA, params)
%   age in years, PGA in m/s^2, params from default_params().ACS

    p = params.ACS;
    mu_corr = p.mu0 + corrosion_fun(age, p.Cr_red, p.dt, p.maint_effectiveness);
    sigma_corr = p.sigma0; % keep initial dispersion; mean updated by corrosion_fun
    p_corr_exceed = 1 - normcdf(0.8*p.t_ACS, mu_corr, sigma_corr);
    % probability margin (0.8*t - corrosion) is below zero
    p_ccf = normcdf(0, p.t_ACS*0.8 - mu_corr, sigma_corr);

    c_p = min(90, mu_corr/p.t * 100);
    C = p.C0 * ((-0.009611)*c_p + 1) * p.capacity_factor;
    R_f = poly_r(p.R_poly, PGA);
    F_seismic = 1 - normcdf((log(C) - log(R_f)) / sqrt(p.sigma_C^2 + p.sigma_R^2));

    p_branch = p_corr_exceed + (1 - p_corr_exceed) * F_seismic;
    p_fail = p_ccf + (1 - p_ccf) * p_branch;
end

function R = poly_r(poly, PGA)
    R = poly(1).*PGA.^2 + poly(2).*PGA + poly(3);
end
