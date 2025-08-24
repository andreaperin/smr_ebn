using UncertaintyQuantification
using CSV
using DataFrames

t_loca = RandomVariable(LogNormal(3.3, 1), :t_loca)
path = "/Applications/MATLAB_R2025a.app/bin/matlab"
source = "Simulation_model"
args = "-batch"
solver1 = Solver(path, source, args)
sourcedir = "/Users/andreaperin_macos/Documents/Code/Stefano_SMR/modelSMR" #{absolute path of /m file }
sources = ["Failure_model_outputs.csv"]
solver = solver1
workdir = "/Users/andreaperin_macos/Documents/Code/Stefano_SMR/runs/temp"
extras = ["Simulation_model.m", "SMDFR_Parameters.m", "SMDFR_HTE_model.slx", "msfcn_indirectps_v1.m", "msfcn_limintm_v3.m", "msfcn_schedule.m", "SMDFR_lib.slx"]
cleanup = false
function extract_function(base_path::String)
    data = CSV.read(joinpath(base_path, "Simulation_model_outputs.csv"), DataFrame)
    return data.T_W1
end
extrator = Extractor(extract_function, :T_W1)
model = ExternalModel(sourcedir, sources, extrator, solver; extras=extras, workdir=workdir, cleanup=false)

df = sample([t_loca], 2)

## SOLVER Definition

path = "/Applications/MATLAB_R2025a.app/bin/matlab"
source = "Simulation_model"
args = "-batch"
solver1 = Solver(path, source, args)

sourcedir = "/Users/andreaperin_macos/Documents/Code/Stefano_SMR/modelSMR" #{absolute path of /m file }
sources = ["Failure_model_outputs.csv"]
solver = solver1
workdir = "/Users/andreaperin_macos/Documents/Code/Stefano_SMR/modelSMR/temp"
extras = ["Simulation_model.m", "SMDFR_Parameters.m", "SMDFR_HTE_model.slx", "msfcn_indirectps_v1.m", "msfcn_limintm_v3.m", "msfcn_schedule.m", "SMDFR_lib.slx"]
cleanup = false

function extract_function(base_path::String)
    data = CSV.read(joinpath(base_path, "Simulation_model_outputs.csv"), DataFrame)
    return data.T_W1
end

extrator = Extractor(extract_function, :T_W1)


model = ExternalModel(sourcedir, sources, extrator, solver; extras=extras, workdir=workdir, cleanup=true)
performance = df -> 1000 .- maximum(df.T_W1)

inputs = [t_loca]
@time probability_of_failure(model, performance, inputs, MonteCarlo(10))

