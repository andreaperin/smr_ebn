function p_fail = edg_failure_prob(age, PGA, params)
%EDG_FAILURE_PROB Deterministic approximation of EDG_1 failure probability.
%   p_fail = EDG_FAILURE_PROB(age, PGA, params)
%   age in years, PGA in m/s^2, params from default_params().EDG

    p = params.EDG;
    mu_corr = p.mu0 + corrosion_fun(age, p.rate, p.dt, p.maint_effectiveness);
    sigma_corr = p.sigma0;
    p_corr_exceed = 1 - normcdf(0.8*p.t_diesel, mu_corr, sigma_corr);

    c_p = min(80, mu_corr/p.t_diesel*100);
    A = p.A0 * ((-0.009611)*c_p + 1) * p.capacity_factor;
    F_seismic = normcdf((log(PGA./A) + p.bu*normcdf(p.Q)) / p.br);

    p_fail = p_corr_exceed + (1 - p_corr_exceed) * F_seismic;
end
