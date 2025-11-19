function p_fail = mslb_failure_prob(age, PGA, params)
%MSLB_FAILURE_PROB Deterministic approximation of MSLB1 probability.
%   p_fail = MSLB_FAILURE_PROB(age, PGA, params)
%   age in years, PGA in m/s^2, params from default_params().MSLB

    p = params.MSLB;
    mu_corr = p.mu0 + corrosion_fun(age, p.Cr_steam, p.dt, p.maint_effectiveness);
    sigma_corr = p.sigma0;
    p_corr_exceed = 1 - normcdf(0.8*p.t_pipe, mu_corr, sigma_corr);

    c_p = min(90, mu_corr/p.t * 100);
    C = p.C0 * ((-0.009611)*c_p + 1) * p.capacity_factor;
    R_f = poly_r(p.R_poly, PGA);
    F_seismic = 1 - normcdf((log(C) - log(R_f)) / sqrt(p.sigma_C^2 + p.sigma_R^2));

    p_fail = p_corr_exceed + (1 - p_corr_exceed) * F_seismic;
end

function R = poly_r(poly, PGA)
    R = poly(1).*PGA.^2 + poly(2).*PGA + poly(3);
end
