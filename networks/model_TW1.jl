
ENV["MATLAB_HOME"] = expanduser("~/Matlab/R2025b")

using MATLAB
using DataFrames

function model_T_W1(LOCA1_time::Float64)
    # Add MATLAB function folder to path
    mat"""
    addpath('/home/perin/Documents/projects/work/code/smr_ebn/modelSMR')
    """
    # Call MATLAB function
    T_W1 = mat"""
    Simulation_model($LOCA1_time)
    """
    return T_W1
end


## TEST

using UncertaintyQuantification
loca = RandomVariable(Normal(), :LOCA)
model = Model(df -> model_T_W1.(df.LOCA), :T_W1)

samples = sample([loca], 20)
evaluate!(model, samples)