using Distributed
numProcs = 5
addprocs(numProcs)

@everywhere begin
    using EnhancedBayesianNetworks
    using EnhancedBayesianNetworks: evaluate!
    using CSV
    using DataFrames
    using JLD2
    using Dates

    const SIMULATIONS = 5 # number of Monte Carlo simulations
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


    ``` LOCA nodes
    ```

    function build_binary_child_from_csv(child_sym::Symbol, csv_path::AbstractString; yes_label::Symbol, no_label::Symbol)
        raw = CSV.read(csv_path, DataFrame)
        rename!(raw, :Column1 => :PGA)
        rename!(raw, :AGE_0 => :AGE_00)

        # long tables for YES and NO
        df_yes = stack(raw, Not(:PGA), variable_name=:AGE, value_name=:Π)  # Π = P(YES)
        df_no = deepcopy(df_yes)
        df_no.Π .= 1 .- df_no.Π                                           # Π = P(NO)

        # convert parent labels to Symbols to match node state ids
        df_yes[!, [:PGA, :AGE]] .= Symbol.(df_yes[!, [:PGA, :AGE]])
        df_no[!, [:PGA, :AGE]] .= Symbol.(df_no[!, [:PGA, :AGE]])

        # insert child outcome column
        insertcols!(df_yes, 3, child_sym => fill(yes_label, nrow(df_yes)))
        insertcols!(df_no, 3, child_sym => fill(no_label, nrow(df_no)))

        df_all = sort(vcat(df_yes, df_no))  # sorted for determinism

        cpt = DiscreteConditionalProbabilityTable{PreciseDiscreteProbability}(df_all)
        node = DiscreteNode(child_sym, cpt)
        return node
    end

    # -------------------------------
    # LOCA1..4 (binary)
    # -------------------------------
    loca1_node = build_binary_child_from_csv(:LOCA1, "networks/csv/LOCA_probability.csv";
        yes_label=:YES_LOCA1, no_label=:NO_LOCA1)

    # If LOCA2..4 share the same CSV structure (possibly same file or different files), call the helper:
    loca2_node = build_binary_child_from_csv(:LOCA2, "networks/csv/LOCA_probability.csv";
        yes_label=:YES_LOCA2, no_label=:NO_LOCA2)

    loca3_node = build_binary_child_from_csv(:LOCA3, "networks/csv/LOCA_probability.csv";
        yes_label=:YES_LOCA3, no_label=:NO_LOCA3)

    loca4_node = build_binary_child_from_csv(:LOCA4, "networks/csv/LOCA_probability.csv";
        yes_label=:YES_LOCA4, no_label=:NO_LOCA4)

    # -------------------------------
    # LOCA-TIME nodes (continuous, child of LOCAi)
    # -------------------------------

    t_loca1_cpt = ContinuousConditionalProbabilityTable{PreciseContinuousInput}(:LOCA1)
    t_loca1_cpt[:LOCA1=>:YES_LOCA1] = LogNormal(3.3, 1.0)
    t_loca1_cpt[:LOCA1=>:NO_LOCA1] = Normal(1200.0, 0)
    t_loca1_node = ContinuousNode(:t_LOCA1, t_loca1_cpt)

    t_loca2_cpt = ContinuousConditionalProbabilityTable{PreciseContinuousInput}(:LOCA2)
    t_loca2_cpt[:LOCA2=>:YES_LOCA2] = LogNormal(3.3, 1.0)
    t_loca2_cpt[:LOCA2=>:NO_LOCA2] = Normal(1200.0, 0)
    t_loca2_node = ContinuousNode(:t_LOCA2, t_loca2_cpt)

    t_loca3_cpt = ContinuousConditionalProbabilityTable{PreciseContinuousInput}(:LOCA3)
    t_loca3_cpt[:LOCA3=>:YES_LOCA3] = LogNormal(3.3, 1.0)
    t_loca3_cpt[:LOCA3=>:NO_LOCA3] = Normal(1200.0, 0)
    t_loca3_node = ContinuousNode(:t_LOCA3, t_loca3_cpt)

    t_loca4_cpt = ContinuousConditionalProbabilityTable{PreciseContinuousInput}(:LOCA4)
    t_loca4_cpt[:LOCA4=>:YES_LOCA4] = LogNormal(3.3, 1.0)
    t_loca4_cpt[:LOCA4=>:NO_LOCA4] = Normal(1200.0, 0)
    t_loca4_node = ContinuousNode(:t_LOCA4, t_loca4_cpt)

    # -------------------------------
    # ACS1..4 (binary), plus ACS1 time
    # -------------------------------
    acs1_node = build_binary_child_from_csv(:ACS1, "networks/csv/ACS1_probability.csv";
        yes_label=:YES_ACS1, no_label=:NO_ACS1)

    acs2_node = build_binary_child_from_csv(:ACS2, "networks/csv/ACS1_probability.csv";
        yes_label=:YES_ACS2, no_label=:NO_ACS2)

    acs3_node = build_binary_child_from_csv(:ACS3, "networks/csv/ACS1_probability.csv";
        yes_label=:YES_ACS3, no_label=:NO_ACS3)

    acs4_node = build_binary_child_from_csv(:ACS4, "networks/csv/ACS1_probability.csv";
        yes_label=:YES_ACS4, no_label=:NO_ACS4)

    t_acs1_cpt = ContinuousConditionalProbabilityTable{PreciseContinuousInput}(:ACS1)
    t_acs1_cpt[:ACS1=>:YES_ACS1] = Uniform(1.0, 1200.0)
    t_acs1_cpt[:ACS1=>:NO_ACS1] = Normal(1200.0, 0)
    t_acs1_node = ContinuousNode(:t_ACS1, t_acs1_cpt)

    t_acs2_cpt = ContinuousConditionalProbabilityTable{PreciseContinuousInput}(:ACS2)
    t_acs2_cpt[:ACS2=>:YES_ACS2] = Uniform(1.0, 1200.0)
    t_acs2_cpt[:ACS2=>:NO_ACS2] = Normal(1200.0, 0)
    t_acs2_node = ContinuousNode(:t_ACS2, t_acs2_cpt)

    t_acs3_cpt = ContinuousConditionalProbabilityTable{PreciseContinuousInput}(:ACS3)
    t_acs3_cpt[:ACS3=>:YES_ACS3] = Uniform(1.0, 1200.0)
    t_acs3_cpt[:ACS3=>:NO_ACS3] = Normal(1200.0, 0)
    t_acs3_node = ContinuousNode(:t_ACS3, t_acs3_cpt)

    t_acs4_cpt = ContinuousConditionalProbabilityTable{PreciseContinuousInput}(:ACS4)
    t_acs4_cpt[:ACS4=>:YES_ACS4] = Uniform(1.0, 1200.0)
    t_acs4_cpt[:ACS4=>:NO_ACS4] = Normal(1200.0, 0)
    t_acs4_node = ContinuousNode(:t_ACS4, t_acs4_cpt)

    tr_acs1_cpt = ContinuousConditionalProbabilityTable{PreciseContinuousInput}(:ACS1)
    tr_acs1_cpt[:ACS1=>:YES_ACS1] = Normal(1200.0, 0)
    tr_acs1_cpt[:ACS1=>:NO_ACS1] = Uniform(30.0, 90.0)
    tr_acs1_node = ContinuousNode(:tr_ACS1, tr_acs1_cpt)

    tr_acs2_cpt = ContinuousConditionalProbabilityTable{PreciseContinuousInput}(:ACS2)
    tr_acs2_cpt[:ACS2=>:YES_ACS2] = Normal(1200.0, 0)
    tr_acs2_cpt[:ACS2=>:NO_ACS2] = Uniform(30.0, 90.0)
    tr_acs2_node = ContinuousNode(:tr_ACS2, tr_acs2_cpt)

    tr_acs3_cpt = ContinuousConditionalProbabilityTable{PreciseContinuousInput}(:ACS3)
    tr_acs3_cpt[:ACS3=>:YES_ACS3] = Normal(1200.0, 0)
    tr_acs3_cpt[:ACS3=>:NO_ACS3] = Uniform(30.0, 90.0)
    tr_acs3_node = ContinuousNode(:tr_ACS3, tr_acs3_cpt)

    tr_acs4_cpt = ContinuousConditionalProbabilityTable{PreciseContinuousInput}(:ACS4)
    tr_acs4_cpt[:ACS4=>:YES_ACS4] = Normal(1200.0, 0)
    tr_acs4_cpt[:ACS4=>:NO_ACS4] = Uniform(30.0, 90.0)
    tr_acs4_node = ContinuousNode(:tr_ACS4, tr_acs4_cpt)


    # -------------------------------
    # EDG1..2 (binary), plus EDG time
    # -------------------------------
    edg1_node = build_binary_child_from_csv(:EDG1, "networks/csv/EDG_probability.csv";
        yes_label=:YES_EDG1, no_label=:NO_EDG1)

    edg2_node = build_binary_child_from_csv(:EDG2, "networks/csv/EDG_probability.csv";
        yes_label=:YES_EDG2, no_label=:NO_EDG2)

    t_edg1_cpt = ContinuousConditionalProbabilityTable{PreciseContinuousInput}(:EDG1)
    t_edg1_cpt[:EDG1=>:YES_EDG1] = Uniform(1.0, 1200.0)
    t_edg1_cpt[:EDG1=>:NO_EDG1] = Normal(1200.0, 0)
    t_edg1_node = ContinuousNode(:t_EDG1, t_edg1_cpt)

    t_edg2_cpt = ContinuousConditionalProbabilityTable{PreciseContinuousInput}(:EDG2)
    t_edg2_cpt[:EDG2=>:YES_EDG2] = Uniform(1.0, 1200.0)
    t_edg2_cpt[:EDG2=>:NO_EDG2] = Normal(1200.0, 0)
    t_edg2_node = ContinuousNode(:t_EDG2, t_edg2_cpt)

    tr_edg1_cpt = ContinuousConditionalProbabilityTable{PreciseContinuousInput}(:EDG1)
    tr_edg1_cpt[:EDG1=>:YES_EDG1] = Normal(1200.0, 0)
    tr_edg1_cpt[:EDG1=>:NO_EDG1] = Uniform(10.0, 180.0)
    tr_edg1_node = ContinuousNode(:tr_EDG1, tr_edg1_cpt)

    tr_edg2_cpt = ContinuousConditionalProbabilityTable{PreciseContinuousInput}(:EDG2)
    tr_edg2_cpt[:EDG2=>:YES_EDG2] = Normal(1200.0, 0)
    tr_edg2_cpt[:EDG2=>:NO_EDG2] = Uniform(10.0, 180.0)
    tr_edg2_node = ContinuousNode(:tr_EDG2, tr_edg2_cpt)


    # ----- LOOP node (depends only on PGA) -----
    raw = CSV.read("networks/csv/LOOP_probability.csv", DataFrame)
    rename!(raw, names(raw)[1] => :PGA)
    # If your second column isn't already :Π, rename it:

    rename!(raw, names(raw)[2] => :Π)
    # Ensure parent state labels are Symbols (as required by CPT)
    raw.PGA = Symbol.(strip.(String.(raw.PGA)))

    # Build YES / NO rows
    # (Safety) Ensure probability is numeric and PGA remains Symbol after selection
    raw.Π = Float64.(raw.Π)

    df_yes = select(raw, :PGA, :Π)
    insertcols!(df_yes, 2, :LOOP => fill(:YES_LOOP, nrow(df_yes)))

    df_no = select(raw, :PGA)
    df_no.Π = 1 .- raw.Π
    insertcols!(df_no, 2, :LOOP => fill(:NO_LOOP, nrow(df_no)))

    df_loop = sort(vcat(df_yes, df_no))

    loop_cpt = DiscreteConditionalProbabilityTable{PreciseDiscreteProbability}(df_loop)
    loop_node = DiscreteNode(:LOOP, loop_cpt)

    # ----- LOOP time node -----
    t_loop_cpt = ContinuousConditionalProbabilityTable{PreciseContinuousInput}(:LOOP)
    t_loop_cpt[:LOOP=>:YES_LOOP] = LogNormal(3.3, 1.0)
    t_loop_cpt[:LOOP=>:NO_LOOP] = Normal(1200.0, 0)
    t_loop_node = ContinuousNode(:t_LOOP, t_loop_cpt)



    ``` MODEL node
    ```
    include("model_T_complete.jl")

    # --- Parallel wrapper for the plant model (all time inputs) ---
    parallel_model_temperatures(df::DataFrame) =
        pmap(eachrow(df)) do r
            model_temperatures(
                r.t_LOCA1, r.t_LOCA2, r.t_LOCA3, r.t_LOCA4,
                r.t_ACS1, r.t_ACS2, r.t_ACS3, r.t_ACS4,
                r.tr_ACS1, r.tr_ACS2, r.tr_ACS3, r.tr_ACS4,
                r.t_EDG1, r.t_EDG2,
                r.tr_EDG1, r.tr_EDG2,
                r.t_LOOP
            )
        end

    # The Model will call our function with a DataFrame containing all parent values;
    # we return a Vector{DataFrame} (one per simulation). EBN stores it under :max_Ts
    model_temp = Model(parallel_model_temperatures, :max_Ts)


    function performance_function(threshold::Real, df::DataFrame)
        maxval = maximum(Matrix(df))
        return threshold - maxval
    end
    performance = df -> performance_function.(threshold, df.max_Ts)
    sim = MonteCarlo(SIMULATIONS)
    model_node = DiscreteFunctionalNode(:Reactor, [model_temp], performance, sim)


    ``` Enhanced Bayesian Networks
    ```
    nodes = [
        # roots
        age_node, pga_node,

        # LOCA and times
        loca1_node, loca2_node, loca3_node, loca4_node,
        t_loca1_node, t_loca2_node, t_loca3_node, t_loca4_node,

        # ACS and times
        acs1_node, acs2_node, acs3_node, acs4_node,
        t_acs1_node, t_acs2_node, t_acs3_node, t_acs4_node,
        tr_acs1_node, tr_acs2_node, tr_acs3_node, tr_acs4_node,

        # EDG and times
        edg1_node, edg2_node,
        t_edg1_node, t_edg2_node,
        tr_edg1_node, tr_edg2_node,

        # LOOP and times
        loop_node,
        t_loop_node,

        # external model
        model_node,
    ]
    ebn = EnhancedBayesianNetwork(nodes)

    # LOCA parents
    add_child!(ebn, :AGE, :LOCA1)
    add_child!(ebn, :PGA, :LOCA1)
    add_child!(ebn, :AGE, :LOCA2)
    add_child!(ebn, :PGA, :LOCA2)
    add_child!(ebn, :AGE, :LOCA3)
    add_child!(ebn, :PGA, :LOCA3)
    add_child!(ebn, :AGE, :LOCA4)
    add_child!(ebn, :PGA, :LOCA4)

    # Link each LOCA to its time node
    add_child!(ebn, :LOCA1, :t_LOCA1)
    add_child!(ebn, :LOCA2, :t_LOCA2)
    add_child!(ebn, :LOCA3, :t_LOCA3)
    add_child!(ebn, :LOCA4, :t_LOCA4)

    # ACS parents and time
    add_child!(ebn, :AGE, :ACS1)
    add_child!(ebn, :PGA, :ACS1)
    add_child!(ebn, :AGE, :ACS2)
    add_child!(ebn, :PGA, :ACS2)
    add_child!(ebn, :AGE, :ACS3)
    add_child!(ebn, :PGA, :ACS3)
    add_child!(ebn, :AGE, :ACS4)
    add_child!(ebn, :PGA, :ACS4)

    # Link ACSs to their time nodes
    add_child!(ebn, :ACS1, :t_ACS1)
    add_child!(ebn, :ACS2, :t_ACS2)
    add_child!(ebn, :ACS3, :t_ACS3)
    add_child!(ebn, :ACS4, :t_ACS4)
    add_child!(ebn, :ACS1, :tr_ACS1)
    add_child!(ebn, :ACS2, :tr_ACS2)
    add_child!(ebn, :ACS3, :tr_ACS3)
    add_child!(ebn, :ACS4, :tr_ACS4)

    # EDG parents (same parents pattern as ACS/LOCA since you used the same CSV helper)
    add_child!(ebn, :AGE, :EDG1)
    add_child!(ebn, :PGA, :EDG1)
    add_child!(ebn, :AGE, :EDG2)
    add_child!(ebn, :PGA, :EDG2)

    # Link EDG to their time nodes
    add_child!(ebn, :EDG1, :t_EDG1)
    add_child!(ebn, :EDG2, :t_EDG2)
    add_child!(ebn, :EDG1, :tr_EDG1)
    add_child!(ebn, :EDG2, :tr_EDG2)

    add_child!(ebn, :PGA, :LOOP)     # parent
    add_child!(ebn, :LOOP, :t_LOOP)   # time node

    # Reactor depends on all time nodes
    add_child!(ebn, :t_LOCA1, :Reactor)
    add_child!(ebn, :t_LOCA2, :Reactor)
    add_child!(ebn, :t_LOCA3, :Reactor)
    add_child!(ebn, :t_LOCA4, :Reactor)
    add_child!(ebn, :t_ACS1, :Reactor)
    add_child!(ebn, :t_ACS2, :Reactor)
    add_child!(ebn, :t_ACS3, :Reactor)
    add_child!(ebn, :t_ACS4, :Reactor)
    add_child!(ebn, :t_EDG1, :Reactor)
    add_child!(ebn, :t_EDG2, :Reactor)
    add_child!(ebn, :t_LOOP, :Reactor)
    add_child!(ebn, :tr_ACS1, :Reactor)
    add_child!(ebn, :tr_ACS2, :Reactor)
    add_child!(ebn, :tr_ACS3, :Reactor)
    add_child!(ebn, :tr_ACS4, :Reactor)
    add_child!(ebn, :tr_EDG1, :Reactor)
    add_child!(ebn, :tr_EDG2, :Reactor)
    order!(ebn)

    # --- initial time ---
    t0 = time_ns()

end

evaluate!(ebn, false, true)

# Save the network to disk
path_to_ebn = joinpath(current_dir, "networks", "ebn_jld2")
mkpath(path_to_ebn)
ebn_name = Dates.format(now(), "yyyy_mm_dd_HH_MM") * "_" *
           string(model_node.simulation) * ".jld2"
@save joinpath(path_to_ebn, ebn_name) ebn

# Print elapsed time
seconds = (time_ns() - t0) / 1e9
println("Elapsed time: $(round(seconds, digits=3)) s")

rmprocs(workers())