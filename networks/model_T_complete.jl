#ENV["MATLAB_HOME"] = expanduser("~/Matlab/R2025b")
ENV["MATLAB_HOME"] = "/Applications/MATLAB_R2024b.app"

using MATLAB
using DataFrames

function model_temperatures(
    t_LOCA1::Float64, t_LOCA2::Float64, t_LOCA3::Float64, t_LOCA4::Float64,
    t_ACS1::Float64,  t_ACS2::Float64,  t_ACS3::Float64,  t_ACS4::Float64,
    tr_ACS1::Float64, tr_ACS2::Float64, tr_ACS3::Float64, tr_ACS4::Float64,
    t_EDG1::Float64,  t_EDG2::Float64,
    tr_EDG1::Float64, tr_EDG2::Float64,
    t_LOOP::Float64
)
    # Ensure numeric types are Float64 (redundant with annotations, but explicit)
    t_LOCA1 = Float64(t_LOCA1); t_LOCA2 = Float64(t_LOCA2); t_LOCA3 = Float64(t_LOCA3); t_LOCA4 = Float64(t_LOCA4)
    t_ACS1  = Float64(t_ACS1);  t_ACS2  = Float64(t_ACS2);  t_ACS3  = Float64(t_ACS3);  t_ACS4  = Float64(t_ACS4)
    tr_ACS1 = Float64(tr_ACS1); tr_ACS2 = Float64(tr_ACS2); tr_ACS3 = Float64(tr_ACS3); tr_ACS4 = Float64(tr_ACS4)
    t_EDG1  = Float64(t_EDG1);  t_EDG2  = Float64(t_EDG2)
    tr_EDG1 = Float64(tr_EDG1); tr_EDG2 = Float64(tr_EDG2)
    t_LOOP  = Float64(t_LOOP)

    # Reset MATLAB state (optional) and ensure model path is available
    mat"close all; clear all;"
    mat"addpath('/Users/stefanomarchetti/Library/CloudStorage/OneDrive-PolitecnicodiMilano/Python/Cursor/smr_ebn/modelSMR')"

    # Call MATLAB Simulation_model with the expanded signature
    T_W1, T_W2, T_W3, T_W4 = mxcall(
        :Simulation_model_complete, 4,
        t_LOCA1, t_LOCA2, t_LOCA3, t_LOCA4,
        t_ACS1,  t_ACS2,  t_ACS3,  t_ACS4,
        tr_ACS1, tr_ACS2, tr_ACS3, tr_ACS4,
        t_EDG1,  t_EDG2,
        tr_EDG1, tr_EDG2,
        t_LOOP
    )

    df = DataFrame(
        max_T_W1 = maximum(T_W1),
        max_T_W2 = maximum(T_W2),
        max_T_W3 = maximum(T_W3),
        max_T_W4 = maximum(T_W4)
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