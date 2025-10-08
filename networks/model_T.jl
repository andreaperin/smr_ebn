if Sys.islinux()
    ENV["MATLAB_HOME"] = "/usr/local/MATLAB/R2023a"
elseif Sys.isapple()
    ENV["MATLAB_HOME"] = "/Applications/MATLAB_R2024b.app"
else
    error("OS not set up")
end

using MATLAB

using DataFrames

function model_temperatures(LOCA1_time::Float64, ACS1_time::Float64)
    # Add MATLAB function folder to path
    LOCA1_time = Float64(LOCA1_time)
    ACS1_time = Float64(ACS1_time)

    mat"close all; clear all;"

    mat"""
    addpath('./modelSMR')
    """
    # Suppress all MATLAB warnings
    mat"warning('off', 'all');"
    # Disable Simulink autosave recovery (prevents the message)
    mat"set_param(0, 'RecoverAutosave', 'off');"

    # Call MATLAB function with multiple outputs
    T_W1, T_W2, T_W3, T_W4 = mxcall(:Simulation_model, 4, LOCA1_time, ACS1_time)
    df = DataFrame(
        max_T_W1=maximum(T_W1),
        max_T_W2=maximum(T_W2),
        max_T_W3=maximum(T_W3),
        max_T_W4=maximum(T_W4)
    )
    return df
end

# a = model_temperatures(20.0, 20.0)

# @show(a)