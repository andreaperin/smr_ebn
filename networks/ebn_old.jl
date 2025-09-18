using EnhancedBayesianNetworks
using EnhancedBayesianNetworks: evaluate!
using CSV
using JLD2
using Dates

# --- initial time ---
t0 = time_ns()

cpt_pga = DiscreteConditionalProbabilityTable{PreciseDiscreteProbability}(:PGA)
cpt_pga[:PGA=>:PGA_00] = 0.05
cpt_pga[:PGA=>:PGA_01] = 0.05
cpt_pga[:PGA=>:PGA_02] = 0.05
cpt_pga[:PGA=>:PGA_03] = 0.05
cpt_pga[:PGA=>:PGA_04] = 0.05
cpt_pga[:PGA=>:PGA_05] = 0.05
cpt_pga[:PGA=>:PGA_06] = 0.05
cpt_pga[:PGA=>:PGA_07] = 0.05
cpt_pga[:PGA=>:PGA_08] = 0.05
cpt_pga[:PGA=>:PGA_09] = 0.05
cpt_pga[:PGA=>:PGA_10] = 0.05
cpt_pga[:PGA=>:PGA_11] = 0.05
cpt_pga[:PGA=>:PGA_12] = 0.05
cpt_pga[:PGA=>:PGA_13] = 0.05
cpt_pga[:PGA=>:PGA_14] = 0.05
cpt_pga[:PGA=>:PGA_15] = 0.05
cpt_pga[:PGA=>:PGA_16] = 0.05
cpt_pga[:PGA=>:PGA_17] = 0.05
cpt_pga[:PGA=>:PGA_18] = 0.05
cpt_pga[:PGA=>:PGA_19] = 0.05
pga_node = DiscreteNode(:PGA, cpt_pga)

age_cpt = DiscreteConditionalProbabilityTable{PreciseDiscreteProbability}(:AGE)
age_cpt[:AGE=>:AGE_0] = 1 / 6
age_cpt[:AGE=>:AGE_10] = 1 / 6
age_cpt[:AGE=>:AGE_20] = 1 / 6
age_cpt[:AGE=>:AGE_30] = 1 / 6
age_cpt[:AGE=>:AGE_40] = 1 / 6
age_cpt[:AGE=>:AGE_50] = 1 / 6

age_node = DiscreteNode(:AGE, age_cpt)

data = CSV.read("networks/LOCA_probability.csv", DataFrame)
rename!(data, :Column1 => :PGA)
df_long1 = stack(data, Not(:PGA), variable_name=:AGE, value_name=:Π)
df_long2 = deepcopy(df_long1)
df_long2[!, :Π] = 1 .- df_long2[!, :Π]
df_long1[!, [:PGA, :AGE]] = Symbol.(df_long1[!, [:PGA, :AGE]])
df_long2[!, [:PGA, :AGE]] = Symbol.(df_long1[!, [:PGA, :AGE]])
newcol1 = fill(:YES_LOCA, nrow(df_long1))
insertcols!(df_long1, 3, :LOCA => newcol1)
newcol2 = fill(:NO_LOCA, nrow(df_long1))
insertcols!(df_long2, 3, :LOCA => newcol2)
df_loca = sort(vcat(df_long1, df_long2))

loca_cpt = DiscreteConditionalProbabilityTable{PreciseDiscreteProbability}(df_loca)
loca_node = DiscreteNode(:LOCA, loca_cpt)

t_loca_cpt = ContinuousConditionalProbabilityTable{PreciseContinuousInput}(:LOCA)
t_loca_cpt[:LOCA=>:YES_LOCA] = LogNormal(3.3, 1)
t_loca_cpt[:LOCA=>:NO_LOCA] = Normal(1200, 0)
t_loca_node = ContinuousNode(:t_loca, t_loca_cpt)

### MODEL node
current_dir = pwd()
path = "/Applications/MATLAB_R2024b.app/bin/matlab"
source = "Simulation_model"
args = "-batch"
solver1 = Solver(path, source, args)
sourcedir = joinpath(current_dir, "modelSMR")
#sourcedir = "/Users/andreaperin_macos/Documents/Code/Stefano_SMR/modelSMR" #{absolute path of /m file }
sources = ["Failure_model_outputs.csv"]
solver = solver1
workdir = joinpath(current_dir, "runs/temp")


#workdir = "/Users/andreaperin_macos/Documents/Code/Stefano_SMR/runs/temp"
extras = ["Simulation_model.m", "SMDFR_Parameters.m", "SMDFR_HTE_model.slx", "msfcn_indirectps_v1.m", "msfcn_limintm_v3.m", "msfcn_schedule.m", "SMDFR_lib.slx"]
cleanup = true
function extract_function(base_path::String)
    data = CSV.read(joinpath(base_path, "Simulation_model_outputs.csv"), DataFrame)
    return data.T_W1
end
extrator = Extractor(extract_function, :T_W1)
model = ExternalModel(sourcedir, sources, extrator, solver; extras=extras, workdir=workdir, cleanup=false)
performance = df -> 1244 .- maximum(df.T_W1)  # [K]
sim = MonteCarlo(5)
model_node = DiscreteFunctionalNode(:Reactor, [model], performance, sim)

nodes = [loca_node, age_node, pga_node, t_loca_node, model_node]
ebn = EnhancedBayesianNetwork(nodes)
add_child!(ebn, :AGE, :LOCA)
add_child!(ebn, :PGA, :LOCA)
add_child!(ebn, :LOCA, :t_loca)
add_child!(ebn, :t_loca, :Reactor)
order!(ebn)

evaluate!(ebn)

path_to_ebn = joinpath(current_dir, "networks/ebn_jld2")
ebn_name = Dates.format(now(), "yyyy_mm_dd_HH_MM") * "_" * string(model_node.simulation) * ".jld2"

@save joinpath(path_to_ebn, ebn_name) ebn

seconds = (time_ns() - t0) / 1e9
println("Elapsed time: $(round(seconds, digits=3)) s")