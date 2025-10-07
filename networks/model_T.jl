#ENV["MATLAB_HOME"] = expanduser("~/Matlab/R2025b")
ENV["MATLAB_HOME"] = "/Applications/MATLAB_R2024b.app"

using MATLAB
using DataFrames

function model_temperatures(LOCA1_time::Float64, ACS1_time::Float64)
    # Add MATLAB function folder to path
    LOCA1_time = Float64(LOCA1_time)
    ACS1_time = Float64(ACS1_time)
    
    mat"close all; clear all;"
   # mat"""
    #addpath('/home/perin/Documents/projects/work/code/smr_ebn/modelSMR')
    #"""

    mat"""
    addpath('/Users/stefanomarchetti/Library/CloudStorage/OneDrive-PolitecnicodiMilano/Python/Cursor/smr_ebn/modelSMR')
    """

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


# a = model_temperatures(20.0)

# TW1, TW2, TW3, TW4 = model_temperatures(20.0)

# ## TEST

# using UncertaintyQuantification
# loca = RandomVariable(Normal(), :LOCA)
# model = Model(df -> model_temperatures.(df.LOCA), :T_W1)

# samples = sample([loca], 20)
# evaluate!(model, samples)