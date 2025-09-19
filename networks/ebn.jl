using EnhancedBayesianNetworks
using EnhancedBayesianNetworks: evaluate!
using CSV
using JLD2
using Dates

const MATLAB_BIN   = "/Applications/MATLAB_R2024b.app/bin/matlab"
const SIMULATIONS  = 1000
const nan_counter = Ref(0)  # counter for failed simulations
const simulation_index = Ref(0)

#need to add ACS nodes

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

# Read LOCA probabilities from CSV and create the LOCA CPT
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

loca_cpt  = DiscreteConditionalProbabilityTable{PreciseDiscreteProbability}(df_loca)
loca_node = DiscreteNode(:LOCA, loca_cpt)

# Time-to-LOCA conditional distribution
t_loca_cpt = ContinuousConditionalProbabilityTable{PreciseContinuousInput}(:LOCA)
t_loca_cpt[:LOCA => :YES_LOCA] = LogNormal(3.3, 1)
t_loca_cpt[:LOCA => :NO_LOCA] = Normal(1200, 0)
t_loca_node = ContinuousNode(:t_loca, t_loca_cpt)

# Read ACS1 probabilities from CSV and create the ACS1 CPT
data = CSV.read("networks/ACS1_probability.csv", DataFrame)
rename!(data, :Column1 => :PGA)
df_long1 = stack(data, Not(:PGA), variable_name=:AGE, value_name=:Π)
df_long2 = deepcopy(df_long1)
df_long2[!, :Π] = 1 .- df_long2[!, :Π]
df_long1[!, [:PGA, :AGE]] .= Symbol.(df_long1[!, [:PGA, :AGE]])
df_long2[!, [:PGA, :AGE]] .= Symbol.(df_long2[!, [:PGA, :AGE]])
insertcols!(df_long1, 3, :ACS1 => fill(:YES_ACS1, nrow(df_long1)))
insertcols!(df_long2, 3, :ACS1 => fill(:NO_ACS1, nrow(df_long2)))
df_acs1 = sort(vcat(df_long1, df_long2))

acs1_cpt  = DiscreteConditionalProbabilityTable{PreciseDiscreteProbability}(df_acs1)
acs1_node = DiscreteNode(:ACS1, acs1_cpt)

# Time-to-ACS1 failure conditional distribution
t_acs1_cpt = ContinuousConditionalProbabilityTable{PreciseContinuousInput}(:ACS1)
t_acs1_cpt[:ACS1 => :YES_ACS1] = Uniform(1, 1200)
t_acs1_cpt[:ACS1 => :NO_ACS1] = Normal(1200, 0)
t_acs1_node = ContinuousNode(:t_acs1, t_acs1_cpt)

### MODEL node
current_dir = pwd()

# Directory containing Simulation_model.m and other extras
sourcedir = joinpath(current_dir, "modelSMR")

# Files produced by the model that must be copied back
sources = ["Failure_model_outputs.csv"]

# Use a local, non-synced temp directory to avoid OneDrive/ICloud delays
workdir = joinpath(tempdir(), "smr_ebn_runs")
mkpath(workdir)
println("workdir = ", workdir)

# Solver: create a .job file then wait for Simulation_model_outputs.csv
path   = "/bin/bash"
args   = "-lc"
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
    "msfcn_schedule.m"
]
cleanup = true

# Helper to detect if MATLAB is running
function is_matlab_running()
    try
        run(pipeline(`pgrep -x MATLAB`, stdout = devnull))
        return true
    catch
        return false
    end
end

# Persistent MATLAB server (one-time launcher)
function start_matlab_server!(workdir::String, matlab_path::String, sourcedir::String)
    ready_flag = joinpath(workdir, ".server_ready")
    if isfile(ready_flag) && is_matlab_running()
        return
    end
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
rm(ready_flag; force=true)
start_matlab_server!(workdir, MATLAB_BIN, sourcedir)

# Extractor: return the maximum T_W1 value as a scalar
function extract_function(base_path::String)
    file_path = joinpath(base_path, "Simulation_model_outputs.csv")
    simulation_index[] += 1
    try
        data = DataFrame(CSV.File(file_path; select=["T_W1"]))
        if maximum(data.T_W1) > 1243
            println("Maximum temperature exceeds 1243 K")
        end
        return maximum(data.T_W1)   # scalar maximum
    catch
        sleep(1)
        println("Failed loading, waiting 1 second and retrying")
        try
            data = DataFrame(CSV.File(file_path; select=["T_W1"]))
            return maximum(data.T_W1)   # scalar maximum
        catch
            nan_counter[] += 1      # increment counter on failure
            println("Index of failed simulation: ", simulation_index[])
            println("Failed to read or process file: ", file_path)
            return NaN
        end
    end
end
extrator = Extractor(extract_function, :T_W1)

# Define the external model
model = ExternalModel(
    sourcedir,
    sources,
    extrator,
    solver;
    extras  = extras,
    workdir = workdir,
    cleanup = true # true to not save stuff
)

# Performance function: handle empty output by returning NaN
performance = df -> begin
    val_vec = df.T_W1
    if isempty(val_vec) || isnan(val_vec[1])
        return NaN
    else
        return 1243 - val_vec[1]   # subtract scalar
    end
end

# Monte Carlo sampling
sim = MonteCarlo(SIMULATIONS)
model_node = DiscreteFunctionalNode(:Reactor, [model], performance, sim)

# Build the Bayesian network
nodes = [loca_node, age_node, pga_node, t_loca_node, model_node, t_acs1_node, acs1_node]
ebn = EnhancedBayesianNetwork(nodes)
add_child!(ebn, :AGE, :LOCA)
add_child!(ebn, :PGA, :LOCA)
add_child!(ebn, :LOCA, :t_loca)
add_child!(ebn, :t_loca, :Reactor)
add_child!(ebn, :AGE, :ACS1)
add_child!(ebn, :PGA, :ACS1)
add_child!(ebn, :ACS1, :t_acs1)
add_child!(ebn, :t_acs1, :Reactor)
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
println("Number of failed simulations: ", nan_counter[])