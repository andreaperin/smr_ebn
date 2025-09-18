using EnhancedBayesianNetworks
using EnhancedBayesianNetworks: evaluate!
using CSV
using JLD2
using Dates

const MATLAB_BIN = "/Applications/MATLAB_R2024b.app/bin/matlab"

const SIMULATIONS = 500

# --- initial time ---
t0 = time_ns()

# Prior on peak ground acceleration
cpt_pga = DiscreteConditionalProbabilityTable{PreciseDiscreteProbability}(:PGA)
for i in 0:19
    cpt_pga[:PGA => Symbol(:PGA_, lpad(i, 2, '0'))] = 0.05
end
pga_node = DiscreteNode(:PGA, cpt_pga)

# Prior on reactor age
age_cpt = DiscreteConditionalProbabilityTable{PreciseDiscreteProbability}(:AGE)
for a in (0, 10, 20, 30, 40, 50)
    age_cpt[:AGE => Symbol(:AGE_, a)] = 1/6
end
age_node = DiscreteNode(:AGE, age_cpt)

# Read LOCA probabilities from CSV, stack and combine into CPT
data = CSV.read("networks/LOCA_probability.csv", DataFrame)
rename!(data, :Column1 => :PGA)
df_long1 = stack(data, Not(:PGA), variable_name=:AGE, value_name=:Π)
df_long2 = deepcopy(df_long1)
df_long2[!, :Π] = 1 .- df_long2[!, :Π]
df_long1[!, [:PGA, :AGE]] .= Symbol.(df_long1[!, [:PGA, :AGE]])
df_long2[!, [:PGA, :AGE]] .= Symbol.(df_long2[!, [:PGA, :AGE]])
insertcols!(df_long1, 3, :LOCA => fill(:YES_LOCA, nrow(df_long1)))
insertcols!(df_long2, 3, :LOCA => fill(:NO_LOCA, nrow(df_long2)))
df_loca = sort(vcat(df_long1, df_long2))

loca_cpt = DiscreteConditionalProbabilityTable{PreciseDiscreteProbability}(df_loca)
loca_node = DiscreteNode(:LOCA, loca_cpt)

# Time-to-LOCA conditional distribution
t_loca_cpt = ContinuousConditionalProbabilityTable{PreciseContinuousInput}(:LOCA)
t_loca_cpt[:LOCA => :YES_LOCA] = LogNormal(3.3, 1)
t_loca_cpt[:LOCA => :NO_LOCA] = Normal(1200, 0)
t_loca_node = ContinuousNode(:t_loca, t_loca_cpt)

### MODEL node
current_dir = pwd()

# Directory containing Simulation_model.m and other extras
sourcedir = joinpath(current_dir, "modelSMR")
# If you need an absolute path instead, uncomment and modify the next line
# sourcedir = "/absolute/path/to/modelSMR"

# Files produced by the model that must be copied back by ExternalModel
sources = ["Failure_model_outputs.csv"]

# Use a local, non-synced temp directory to avoid OneDrive/ICloud delays
workdir = joinpath(tempdir(), "smr_ebn_runs")
mkpath(workdir)

println("workdir = ", workdir)

# Solver: create a job file and then wait for Simulation_model_outputs.csv
path = "/bin/bash"
args = "-lc"
source = "rm -f Simulation_model_outputs.csv && touch run.job && until [ -f Simulation_model_outputs.csv ]; do sleep 0.05; done"
solver1 = Solver(path, source, args)
solver  = solver1

# Files needed to run the Simulink model
extras = [
    "Simulation_model.m",
    "SMDFR_Parameters.m",
    "SMDFR_HTE_model.slx",
    "msfcn_indirectps_v1.m",
    "msfcn_limintm_v3.m",
    "msfcn_schedule.m",
    # "SMDFR_lib.slx"
]
cleanup = true

# Persistent MATLAB server (one-time launcher)
function is_matlab_running()
    try
        # pgrep -x MATLAB matches the MATLAB process on macOS and Linux
        run(pipeline(`pgrep -x MATLAB`, stdout = devnull))
        return true
    catch
        return false
    end
end

function start_matlab_server!(workdir::String, matlab_path::String, sourcedir::String)
    ready_flag = joinpath(workdir, ".server_ready")
    if isfile(ready_flag) && is_matlab_running()
        return  # MATLAB is already running
    end
    # remove stale flag and launch server
    rm(ready_flag; force=true)
    mkpath(workdir)
    matlab_arg = "cd('" * sourcedir * "'); Simulation_model('server','" * workdir * "')"
    cmd = `$(matlab_path) -batch $(matlab_arg)`
    @async run(cmd)
    open(ready_flag, "w") do io
        write(io, string(Dates.now()))
    end
end

# Start the persistent server
ready_flag = joinpath(workdir, ".server_ready")
if isfile(ready_flag)
    rm(ready_flag)
end

start_matlab_server!(workdir, MATLAB_BIN, sourcedir)

# Extractor: read the T_W1 column from the simulation outputs
function extract_function(base_path::String)
    file_path = joinpath(base_path, "Simulation_model_outputs.csv")
    data = DataFrame(CSV.File(file_path; select = ["T_W1"]))
    return data.T_W1
end
extrator = Extractor(extract_function, :T_W1)

# Define the external model
model = ExternalModel(
    sourcedir,
    sources,
    extrator,
    solver;
    extras = extras,
    workdir = workdir,
    cleanup = false
)
# Performance function: compute margin from 1244 K
performance = df -> 1244 .- maximum(df.T_W1)  # [K]

# Monte Carlo sampling
sim = MonteCarlo(SIMULATIONS)
model_node = DiscreteFunctionalNode(:Reactor, [model], performance, sim)

# Build the Bayesian network
nodes = [loca_node, age_node, pga_node, t_loca_node, model_node]
ebn = EnhancedBayesianNetwork(nodes)
add_child!(ebn, :AGE, :LOCA)
add_child!(ebn, :PGA, :LOCA)
add_child!(ebn, :LOCA, :t_loca)
add_child!(ebn, :t_loca, :Reactor)
order!(ebn)

# Evaluate the network
evaluate!(ebn)

# Save the network to disk
path_to_ebn = joinpath(current_dir, "networks", "ebn_jld2")
mkpath(path_to_ebn)
ebn_name = Dates.format(now(), "yyyy_mm_dd_HH_MM") * "_" *
           string(model_node.simulation) * ".jld2"
@save joinpath(path_to_ebn, ebn_name) ebn

# Print elapsed time
seconds = (time_ns() - t0) / 1e9
println("Elapsed time: $(round(seconds, digits=3)) s")