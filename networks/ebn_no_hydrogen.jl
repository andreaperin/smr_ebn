# using Distributed
# numProcs = 6
# addprocs(numProcs)

# @everywhere begin
    using EnhancedBayesianNetworks
    using EnhancedBayesianNetworks: evaluate!
    using CSV
    using DataFrames
    using JLD2
    using Dates
    
    const SIMULATIONS = 50 # number of Monte Carlo simulations
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
    
    
    """ LOCA node
    """
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
    
    
    ``` LOOP node
    ```
    data = CSV.read("networks/csv/LOOP_probability.csv", DataFrame)
    rename!(data, :Column1 => :PGA)
    df_loop = sort(vcat(
        DataFrame(PGA=Symbol.(data.PGA), LOOP=:YES_LOOP, Π=Real.(data.LOOP)),
        DataFrame(PGA=Symbol.(data.PGA), LOOP=:NO_LOOP, Π=Real.(1 .- data.LOOP))
    ))
    loop_cpt = DiscreteConditionalProbabilityTable{PreciseDiscreteProbability}(df_loop)
    loop_node = DiscreteNode(:LOOP, loop_cpt)
    
    
    ``` LOOP-TIME node
    ```
    t_loop_cpt = ContinuousConditionalProbabilityTable{PreciseContinuousInput}(:LOOP)
    t_loop_cpt[:LOOP=>:YES_LOOP] = LogNormal(3.3, 1)
    t_loop_cpt[:LOOP=>:NO_LOOP] = Normal(1200, 0)
    t_loop_node = ContinuousNode(:t_loop, t_loop_cpt)
    
    
    ``` MSLB node
    ```
    data = CSV.read("networks/csv/MSLB_probability.csv", DataFrame)
    rename!(data, :Column1 => :PGA)
    rename!(data, :AGE_0 => :AGE_00)
    df_long1 = stack(data, Not(:PGA), variable_name=:AGE, value_name=:Π)
    df_long2 = deepcopy(df_long1)
    df_long2[!, :Π] = 1 .- df_long2[!, :Π]
    df_long1[!, [:PGA, :AGE]] .= Symbol.(df_long1[!, [:PGA, :AGE]])
    df_long2[!, [:PGA, :AGE]] .= Symbol.(df_long2[!, [:PGA, :AGE]])
    insertcols!(df_long1, 3, :MSLB => fill(:YES_MSLB, nrow(df_long1)))
    insertcols!(df_long2, 3, :MSLB => fill(:NO_MSLB, nrow(df_long2)))
    df_mslb = sort(vcat(df_long1, df_long2))
    mslb_cpt = DiscreteConditionalProbabilityTable{PreciseDiscreteProbability}(df_mslb)
    mslb_node = DiscreteNode(:MSLB, mslb_cpt)
    
    
    ``` MSLB-TIME node
    ```
    t_mslb_cpt = ContinuousConditionalProbabilityTable{PreciseContinuousInput}(:MSLB)
    t_mslb_cpt[:MSLB=>:YES_MSLB] = LogNormal(3.3, 1)
    t_mslb_cpt[:MSLB=>:NO_MSLB] = Normal(1200, 0)
    t_mslb_node = ContinuousNode(:t_mslb, t_mslb_cpt)
    
    
    ``` ACS node
    ```
    data = CSV.read("networks/csv/ACS_probability.csv", DataFrame)
    rename!(data, :Column1 => :PGA)
    rename!(data, :AGE_0 => :AGE_00)
    df_long1 = stack(data, Not(:PGA), variable_name=:AGE, value_name=:Π)
    df_long2 = deepcopy(df_long1)
    df_long2[!, :Π] = 1 .- df_long2[!, :Π]
    df_long1[!, [:PGA, :AGE]] .= Symbol.(df_long1[!, [:PGA, :AGE]])
    df_long2[!, [:PGA, :AGE]] .= Symbol.(df_long2[!, [:PGA, :AGE]])
    insertcols!(df_long1, 3, :ACS => fill(:YES_ACS, nrow(df_long1)))
    insertcols!(df_long2, 3, :ACS => fill(:NO_ACS, nrow(df_long2)))
    df_acs = sort(vcat(df_long1, df_long2))
    acs_cpt = DiscreteConditionalProbabilityTable{PreciseDiscreteProbability}(df_acs)
    acs_node = DiscreteNode(:ACS, acs_cpt)
    
    
    ``` ACS-TIME node
    ```
    t_acs_cpt = ContinuousConditionalProbabilityTable{PreciseContinuousInput}(:ACS)
    t_acs_cpt[:ACS=>:YES_ACS] = Uniform(1, 1200)
    t_acs_cpt[:ACS=>:NO_ACS] = Normal(1200, 0)
    t_acs_node = ContinuousNode(:t_acs, t_acs_cpt)
    
    ``` ACS-rTIME node
    ```
    rt_acs_cpt = ContinuousConditionalProbabilityTable{PreciseContinuousInput}()
    rt_acs_cpt[] = Uniform(30.0, 90.0)
    rt_acs_node = ContinuousNode(:rt_acs, rt_acs_cpt)
    
    
    ``` EDG node
    ```
    data = CSV.read("networks/csv/EDG_probability.csv", DataFrame)
    rename!(data, :Column1 => :PGA)
    rename!(data, :AGE_0 => :AGE_00)
    df_long1 = stack(data, Not(:PGA), variable_name=:AGE, value_name=:Π)
    df_long2 = deepcopy(df_long1)
    df_long2[!, :Π] = 1 .- df_long2[!, :Π]
    df_long1[!, [:PGA, :AGE]] .= Symbol.(df_long1[!, [:PGA, :AGE]])
    df_long2[!, [:PGA, :AGE]] .= Symbol.(df_long2[!, [:PGA, :AGE]])
    insertcols!(df_long1, 3, :EDG => fill(:YES_EDG, nrow(df_long1)))
    insertcols!(df_long2, 3, :EDG => fill(:NO_EDG, nrow(df_long2)))
    df_edg = sort(vcat(df_long1, df_long2))
    edg_cpt = DiscreteConditionalProbabilityTable{PreciseDiscreteProbability}(df_edg)
    edg_node = DiscreteNode(:EDG, edg_cpt)
    
    
    ``` EDG-TIME node
    ```
    t_edg_cpt = ContinuousConditionalProbabilityTable{PreciseContinuousInput}(:EDG)
    t_edg_cpt[:EDG=>:YES_EDG] = Uniform(1, 1200)
    t_edg_cpt[:EDG=>:NO_EDG] = Normal(1200, 0)
    t_edg_node = ContinuousNode(:t_edg, t_edg_cpt)
    
    ``` EDG-rTIME node
    ```
    rt_edg_cpt = ContinuousConditionalProbabilityTable{PreciseContinuousInput}()
    rt_edg_cpt[] = Uniform(10.0, 180.0)
    rt_edg_node = ContinuousNode(:rt_edg, rt_edg_cpt)
    
    
    
    
    ``` MODEL node
    ```
    include("model_T_no_hydrogen.jl")
    
    # model_temp = ParallelModel(df -> model_temperatures(df.t_loca, df.t_loop, df.t_LHS, df.t_mslb, df.t_MSLBH2, df.t_LOOPH2, df.t_acs, df.rt_acs, df.t_edg, df.rt_edg, df.t_pdp, df.rt_pdp), :max_Ts)
    
    model_temp = Model(df -> model_temperatures_parallel(df.t_loca, df.t_loop, df.t_LHS, df.t_mslb, df.t_MSLBH2, df.t_LOOPH2, df.t_acs, df.rt_acs, df.t_edg, df.rt_edg, df.t_pdp, df.rt_pdp), :maxT_W1)
    
    # function performance_function(threshold::Real, df::DataFrame)
    #     maxval = maximum(Matrix(df))
    #     return threshold - maxval
    # end
    
    # performance = df -> performance_function.(threshold, df.max_Ts)
    
    performance = df -> threshold .- df.maxT_W1
    sim = MonteCarlo(SIMULATIONS)
    model_node = DiscreteFunctionalNode(:Reactor, [model_temp], performance, sim)
    
    
    ``` Enhanced Bayesian Networks
    ```
    nodes = [
        age_node, pga_node,
        loca_node, t_loca_node,
        loop_node, t_loop_node,
        mslb_node, t_mslb_node,
        acs_node, t_acs_node, rt_acs_node,
        edg_node, t_edg_node, rt_edg_node,
        model_node
    ]
    
    ebn = EnhancedBayesianNetwork(nodes)
    add_child!(ebn, :AGE, :LOCA)
    add_child!(ebn, :PGA, :LOCA)
    add_child!(ebn, :LOCA, :t_loca)
    
    add_child!(ebn, :PGA, :LOOP)
    add_child!(ebn, :LOOP, :t_loop)
    
    add_child!(ebn, :AGE, :MSLB)
    add_child!(ebn, :PGA, :MSLB)
    add_child!(ebn, :MSLB, :t_mslb)
    
    add_child!(ebn, :AGE, :ACS)
    add_child!(ebn, :PGA, :ACS)
    add_child!(ebn, :ACS, :t_acs)
    
    add_child!(ebn, :AGE, :EDG)
    add_child!(ebn, :PGA, :EDG)
    add_child!(ebn, :EDG, :t_edg)
    
    add_child!(ebn, :AGE, :PDP)
    add_child!(ebn, :PGA, :PDP)
    
    add_child!(ebn, :t_loca, :Reactor)
    add_child!(ebn, :t_loop, :Reactor)
    add_child!(ebn, :t_mslb, :Reactor)
    add_child!(ebn, :t_acs, :Reactor)
    add_child!(ebn, :rt_acs, :Reactor)
    add_child!(ebn, :t_edg, :Reactor)
    add_child!(ebn, :rt_edg, :Reactor)
    
    
    order!(ebn)
    gplot(ebn; NODESIZEFACTOR=0.1, ARROWLENGTH=0.05, NODELABELSIZE=2.5)
    # --- initial time ---
    t0 = time_ns()
    
    # end
    
    evaluate!(ebn, false, false)
    
    # Save the network to disk
    path_to_ebn = joinpath(current_dir, "networks", "ebn_jld2")
    mkpath(path_to_ebn)
    ebn_name = Dates.format(now(), "yyyy_mm_dd_HH_MM") * "_" *
               string(model_node.simulation) * ".jld2"
    @save joinpath(path_to_ebn, ebn_name) ebn
    
    # Print elapsed time
    seconds = (time_ns() - t0) / 1e9
    println("Elapsed time: $(round(seconds, digits=3)) s")
    
    # rmprocs(workers())