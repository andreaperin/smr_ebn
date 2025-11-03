# if Sys.islinux()
#     ENV["MATLAB_HOME"] = expanduser("~/Matlab/R2025b")
# elseif Sys.isapple()
#     ENV["MATLAB_HOME"] = "/Applications/MATLAB_R2024b.app"
# else
#     error("OS not set up")
# end

using MATLAB

using DataFrames

function model_temperatures(
    LOCA_time::Float64,
    LOOP_time::Float64,
    LHS_time::Float64,
    MSLB_time::Float64,
    MSLBH2_time::Float64,
    LOOPH2_time::Float64,
    ACS_time::Float64,
    ACS_rtime::Float64,
    EDG_time::Float64,
    EDG_rtime::Float64,
    pdp_time::Float64,
    pdp_rtime::Float64
)

    mat"close all; clear all;"

    mat"""
    addpath('./modelSMR')
    """
    # Suppress all MATLAB warnings
    mat"warning('off', 'all');"

    # Call MATLAB function with multiple outputs
    T_W1, T_W2, T_W3, T_W4 = mxcall(:Simulation_model_hydrogen, 4,
        LOCA_time,
        LOOP_time,
        LHS_time,
        MSLB_time,
        MSLBH2_time,
        LOOPH2_time,
        ACS_time,
        ACS_rtime,
        EDG_time,
        EDG_rtime,
        pdp_time,
        pdp_rtime
    )
    df = DataFrame(
        max_T_W1=maximum(T_W1),
        max_T_W2=maximum(T_W2),
        max_T_W3=maximum(T_W3),
        max_T_W4=maximum(T_W4)
    )
    return df
end


function model_temperatures_parallel(
    LOCA_time::Vector{Float64},
    LOOP_time::Vector{Float64},
    MSLB_time::Vector{Float64},
    ACS_time::Vector{Float64},
    ACS_rtime::Vector{Float64},
    EDG_time::Vector{Float64},
    EDG_rtime::Vector{Float64},
)

    mat"close all; clear all;"

    mat"""
    addpath('./modelSMR')
    """
    # Suppress all MATLAB warnings
    mat"warning('off', 'all');"

    # Call MATLAB function with multiple outputs
    T_W1, T_W2, T_W3, T_W4 = mxcall(:Simulation_model_hydrogen_parallel, 4,
        LOCA_time,
        LOOP_time,
        1200,
        MSLB_time,
        1200,
        1200,
        ACS_time,
        ACS_rtime,
        EDG_time,
        EDG_rtime,
        1200,
        0
    )

    # df = DataFrame(
    #     max_T_W1=vec(maximum(T_W1, dims=2)),
    #     max_T_W2=vec(maximum(T_W2, dims=2)),
    #     max_T_W3=vec(maximum(T_W3, dims=2)),
    #     max_T_W4=vec(maximum(T_W4, dims=2)),
    # )

    max_T_W1 = vec(maximum(T_W1, dims=2))
    # max_T_W2 = vec(maximum(T_W2, dims=2))
    # max_T_W3 = vec(maximum(T_W3, dims=2))
    # max_T_W4 = vec(maximum(T_W4, dims=2))

    return max_T_W1
end

# a = model_temperatures(20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0)

# val = [20.0, 30.0, 40.0, 50.0, 60.0, 20.0, 25.0, 30.0, 10.0, 155.0, 20.0, 30.0, 40.0, 50.0, 60.0, 20.0, 25.0, 30.0, 10.0, 155.0, 20.0, 30.0, 40.0, 50.0, 60.0, 20.0, 25.0, 30.0, 10.0, 155.0, 20.0, 30.0, 40.0, 50.0, 60.0, 20.0, 25.0, 30.0, 10.0, 155.0, 20.0, 30.0, 40.0, 50.0, 60.0, 20.0, 25.0, 30.0, 10.0, 155.0, 20.0, 30.0, 40.0, 50.0, 60.0, 20.0, 25.0, 30.0, 10.0, 155.0, 20.0, 30.0, 40.0, 50.0, 60.0, 20.0, 25.0, 30.0, 10.0, 155.0, 20.0, 30.0, 40.0, 50.0, 60.0, 20.0, 25.0, 30.0, 10.0, 155.0, 20.0, 30.0, 40.0, 50.0, 60.0, 20.0, 25.0, 30.0, 10.0, 155.0, 20.0, 30.0, 40.0, 50.0, 60.0, 20.0, 25.0, 30.0, 10.0, 155.0, 20.0, 30.0, 40.0, 50.0, 60.0, 20.0, 25.0, 30.0, 10.0, 155.0, 20.0, 30.0, 40.0, 50.0, 60.0, 20.0, 25.0, 30.0, 10.0, 155.0, 20.0, 30.0, 40.0, 50.0, 60.0, 20.0, 25.0, 30.0, 10.0, 155.0, 20.0, 30.0, 40.0, 50.0, 60.0, 20.0, 25.0, 30.0, 10.0, 155.0, 20.0, 30.0, 40.0, 50.0, 60.0, 20.0, 25.0, 30.0, 10.0, 155.0, 20.0, 30.0, 40.0, 50.0, 60.0, 20.0, 25.0, 30.0, 10.0, 155.0, 20.0, 30.0, 40.0, 50.0, 60.0, 20.0, 25.0, 30.0, 10.0, 155.0, 20.0, 30.0, 40.0, 50.0, 60.0, 20.0, 25.0, 30.0, 10.0, 155.0, 20.0, 30.0, 40.0, 50.0, 60.0, 20.0, 25.0, 30.0, 10.0, 155.0, 20.0, 30.0, 40.0, 50.0, 60.0, 20.0, 25.0, 30.0, 10.0, 155.0, 20.0, 30.0, 40.0, 50.0, 60.0, 20.0, 25.0, 30.0, 10.0, 155.0, 20.0, 30.0, 40.0, 50.0, 60.0, 20.0, 25.0, 30.0, 10.0, 155.0, 20.0, 30.0, 40.0, 50.0, 60.0, 20.0, 25.0, 30.0, 10.0, 155.0, 20.0, 30.0, 40.0, 50.0, 60.0, 20.0, 25.0, 30.0, 10.0, 155.0]
# val = [20.0, 30.0, 40.0, 50.0]

# LOCA_time = val
# LOOP_time = val
# LHS_time = val
# MSLB_time = val
# MSLBH2_time = val
# LOOPH2_time = val
# ACS_time = val
# ACS_rtime = val
# EDG_time = val
# EDG_rtime = val
# pdp_time = val
# pdp_rtime = val
# b = model_temperatures_parallel(val, val, val, val, val, val, val, val, val, val, val, val)

# @show(a)


# using UncertaintyQuantification

# t_loca = RandomVariable(Uniform(20, 100), :t_loca)
# t_loop = RandomVariable(Uniform(20, 100), :t_loop)
# t_LHS = RandomVariable(Uniform(20, 100), :t_LHS)
# t_mslb = RandomVariable(Uniform(20, 100), :t_mslb)
# t_MSLBH2 = RandomVariable(Uniform(20, 100), :t_MSLBH2)
# t_LOOPH2 = RandomVariable(Uniform(20, 100), :t_LOOPH2)
# t_acs = RandomVariable(Uniform(20, 100), :t_acs)
# rt_acs = RandomVariable(Uniform(20, 100), :rt_acs)
# t_edg = RandomVariable(Uniform(20, 100), :t_edg)
# rt_edg = RandomVariable(Uniform(20, 100), :rt_edg)
# t_pdp = RandomVariable(Uniform(20, 100), :t_pdp)
# rt_pdp = RandomVariable(Uniform(20, 100), :rt_pdp)

# inputs = [t_loca, t_loop, t_LHS, t_mslb, t_MSLBH2, t_LOOPH2, t_acs, rt_acs, t_edg, rt_edg, t_pdp, rt_pdp]

# model_temp = Model(df -> model_temperatures_parallel(df.t_loca, df.t_loop, df.t_LHS, df.t_mslb, df.t_MSLBH2, df.t_LOOPH2, df.t_acs, df.rt_acs, df.t_edg, df.rt_edg, df.t_pdp, df.rt_pdp), :max_Ts)

# samples = sample(inputs, 10)
# UncertaintyQuantification.evaluate!(model_temp, samples)