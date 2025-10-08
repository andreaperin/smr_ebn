using Distributed
numProcs = 10
addprocs(numProcs)

@everywhere begin
    using EnhancedBayesianNetworks
    using EnhancedBayesianNetworks: evaluate!
    using CSV
    using DataFrames
    using JLD2
    using Dates

    const SIMULATIONS = 10 # number of Monte Carlo simulations
    const threshold = 1243.9

    const current_dir = pwd()


    ``` PGA node peak ground acceleration
    ```
    cpt_pga = DiscreteConditionalProbabilityTable{PreciseDiscreteProbability}(:PGA)
    pga_probs = CSV.read("networks/csv/PGA_probability.csv", DataFrame)[:, 1]
    pga_states = ["PGA_$(lpad(i, 2, '0'))" for i in 0:length(pga_probs)]
    map((p, st) -> cpt_pga[:PGA=>Symbol(st)] = p, pga_probs, pga_states)
    pga_node = DiscreteNode(:PGA, cpt_pga)


    ``` AGE node
    ```
    cpt_age = DiscreteConditionalProbabilityTable{PreciseDiscreteProbability}(:AGE)
    age_states = ["AGE_$(i)0" for i in 0:5]
    map(st -> cpt_age[:AGE=>Symbol(st)] = 1 / length(age_states), age_states)
    age_node = DiscreteNode(:AGE, cpt_age)


    ``` LOCA node
    ```
    data = CSV.read("networks/csv/LOCA_probability.csv", DataFrame)
    rename!(data, :Column1 => :PGA)
    rename!(data, :AGE_0 => :AGE_00)
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


    ``` LOCA-TIME node
    ```
    t_loca_cpt = ContinuousConditionalProbabilityTable{PreciseContinuousInput}(:LOCA)
    t_loca_cpt[:LOCA=>:YES_LOCA] = LogNormal(3.3, 1)
    t_loca_cpt[:LOCA=>:NO_LOCA] = Normal(1200, 0)
    t_loca_node = ContinuousNode(:t_loca, t_loca_cpt)


    ``` ACS1 node
    ```
    data = CSV.read("networks/csv/ACS1_probability.csv", DataFrame)
    rename!(data, :Column1 => :PGA)
    rename!(data, :AGE_0 => :AGE_00)
    df_long1 = stack(data, Not(:PGA), variable_name=:AGE, value_name=:Π)
    df_long2 = deepcopy(df_long1)
    df_long2[!, :Π] = 1 .- df_long2[!, :Π]
    df_long1[!, [:PGA, :AGE]] .= Symbol.(df_long1[!, [:PGA, :AGE]])
    df_long2[!, [:PGA, :AGE]] .= Symbol.(df_long2[!, [:PGA, :AGE]])
    insertcols!(df_long1, 3, :ACS1 => fill(:YES_ACS1, nrow(df_long1)))
    insertcols!(df_long2, 3, :ACS1 => fill(:NO_ACS1, nrow(df_long2)))
    df_acs1 = sort(vcat(df_long1, df_long2))
    acs1_cpt = DiscreteConditionalProbabilityTable{PreciseDiscreteProbability}(df_acs1)
    acs1_node = DiscreteNode(:ACS1, acs1_cpt)


    ``` ACS1-TIME node
    ```
    t_acs1_cpt = ContinuousConditionalProbabilityTable{PreciseContinuousInput}(:ACS1)
    t_acs1_cpt[:ACS1=>:YES_ACS1] = Uniform(1, 1200)
    t_acs1_cpt[:ACS1=>:NO_ACS1] = Normal(1200, 0)
    t_acs1_node = ContinuousNode(:t_acs1, t_acs1_cpt)



    ``` MODEL node
    ```
    include("model_T.jl")

    # --- Parallel wrapper for the plant model (t_loca and t_acs1) ---
    function _run_model_temperatures(t_loca::Real, t_acs1::Real)
        return model_temperatures(t_loca, t_acs1)
    end

    # Map over aligned vectors with pmap (each pair (t_loca, t_acs1) is one simulation)
    parallel_model_temperatures(ts_loca::AbstractVector, ts_acs1::AbstractVector) =
        pmap((tl, ta) -> _run_model_temperatures(tl, ta), ts_loca, ts_acs1)

    model_temp = Model(df -> parallel_model_temperatures(df.t_loca, df.t_acs1), :max_Ts)


    function performance_function(threshold::Real, df::DataFrame)
        maxval = maximum(Matrix(df))
        return threshold - maxval
    end

    performance = df -> performance_function.(threshold, df.max_Ts)
    sim = MonteCarlo(SIMULATIONS)
    model_node = DiscreteFunctionalNode(:Reactor, [model_temp], performance, sim)


    ``` Enhanced Bayesian Networks
    ```
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

    # --- initial time ---
    t0 = time_ns()

end

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
