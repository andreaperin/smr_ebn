using Distributed
n_workers = 6
addprocs(n_workers)

@everywhere begin
    using UncertaintyQuantification
    using MATLAB
    using DataFrames
    using JLD2

    const mc_sim = 10^3
    const threshold = 1243.9

    age = Parameter(50, :AGE)

    shape_pga = 1.11
    scale_pga = 0.000133
    min_pga = 0.1
    max_pga = 2 * 9.81
    pga = RandomVariable(truncated(Frechet(shape_pga, scale_pga), min_pga, max_pga), :PGA)

    rt_acs = RandomVariable(Uniform(30.0, 90), :rt_acs)
    rt_edg = RandomVariable(Uniform(10, 180), :rt_edg)
    rt_pdp = RandomVariable(Uniform(10, 180), :rt_pdp)


    include("networks/model_T.jl")
    model_times = Model(df -> failure_times.(df.PGA, df.AGE), :failure_times)

    f_loca = df -> df.LOCA_time[1]
    model_loca = Model(df -> f_loca.(df.failure_times), :t_loca)

    f_loop = df -> df.LOOP_time[1]
    model_loop = Model(df -> f_loop.(df.failure_times), :t_loop)

    f_lhs = df -> df.LHS_time[1]
    model_LHS = Model(df -> f_lhs.(df.failure_times), :t_LHS)

    f_mslb = df -> df.MSLB_time[1]
    model_MSLB = Model(df -> f_mslb.(df.failure_times), :t_mslb)

    f_mslbh2 = df -> df.MSLBH2_time[1]
    model_MSLBH2 = Model(df -> f_mslbh2.(df.failure_times), :t_MSLBH2)

    f_looph2 = df -> df.LOOPH2_time[1]
    model_LOOPH2 = Model(df -> f_looph2.(df.failure_times), :t_LOOPH2)

    f_acs = df -> df.ACS_time[1]
    model_ACS = Model(df -> f_acs.(df.failure_times), :t_acs)

    f_edg = df -> df.EDG_time[1]
    model_EDG = Model(df -> f_edg.(df.failure_times), :t_edg)

    f_pdp = df -> df.pdp_time[1]
    model_pdp = Model(df -> f_pdp.(df.failure_times), :t_pdp)

    model_temp = Model(df -> model_temperatures_parallel(df.t_loca, df.t_loop, df.t_LHS, df.t_mslb, df.t_MSLBH2, df.t_LOOPH2, df.t_acs, df.rt_acs, df.t_edg, df.rt_edg, df.t_pdp, df.rt_pdp), :maxT_W1)

    models = [model_times, model_loca, model_loop, model_LHS, model_MSLB, model_MSLBH2, model_LOOPH2, model_ACS, model_EDG, model_pdp, model_temp]
    inputs = [age, pga, rt_acs, rt_edg, rt_pdp]
    performance = df -> threshold .- df.maxT_W1

    sim = MonteCarlo(mc_sim)
end

using Dates

@time p_f, var, samples = probability_of_failure(models, performance, inputs, sim)
res = [p_f, var, samples]

path_to_simulation = joinpath(pwd(), "simulations")
mkpath(path_to_simulation)
sim_name = Dates.format(now(), "yyyy_mm_dd_HH_MM") * "_" * string(sim) * ".jld2"
@save joinpath(path_to_simulation, sim_name) res