using EnhancedBayesianNetworks
using JLD2
using DataFrames
using Plots
using CSV

const NETWORK_DIR = joinpath(@__DIR__, "ebn_jld2")
const PLOT_DIR = joinpath(@__DIR__, "imgs")
const DEFAULT_AGE_STATE = :AGE_00
const DEFAULT_DISTANCE_STATE = :Distance_250

function load_last_node_cpt(path::AbstractString)
    data = load(path)
    ebn = haskey(data, "ebn") ? data["ebn"] : error("File $(path) does not contain an `ebn` object")
    node = ebn.nodes[end]
    return (
        path = path,
        node_name = node.name,
        node_type = typeof(node),
        table = node.cpt.data,
    )
end

values_equal(a, b; atol=1e-10, rtol=1e-10) =
    a isa Number && b isa Number ? isapprox(a, b; atol=atol, rtol=rtol) : a == b

function tables_equal(df1::DataFrame, df2::DataFrame; atol=1e-10, rtol=1e-10)
    names(df1) == names(df2) || return false
    nrow(df1) == nrow(df2) || return false
    for row in 1:nrow(df1), col in names(df1)
        values_equal(df1[row, col], df2[row, col]; atol=atol, rtol=rtol) || return false
    end
    return true
end

function mismatching_rows(df1::DataFrame, df2::DataFrame; atol=1e-10, rtol=1e-10)
    names(df1) == names(df2) && nrow(df1) == nrow(df2) || return collect(1:max(nrow(df1), nrow(df2)))
    bad = Int[]
    for row in 1:nrow(df1)
        for col in names(df1)
            values_equal(df1[row, col], df2[row, col]; atol=atol, rtol=rtol) || (push!(bad, row); break)
        end
    end
    return bad
end

function label_slice(str::AbstractString; maxlen=60)
    length(str) <= maxlen && return str
    return string(first(str, maxlen - 3), "…")
end

function row_labels(df::DataFrame)
    keys_mask = names(df) .!= :Π
    map(1:nrow(df)) do i
        pairs = join(["$(names(df)[j])=$(df[i, j])" for j in eachindex(names(df)) if keys_mask[j]], " | ")
        isempty(pairs) ? "row $(i)" : pairs
    end
end

function probability_of_state(ϕ::Factor, node::Symbol, state::Symbol)
    mapping = get(ϕ.states_mapping, node) do
        error("Node $(node) not found in factor dimensions $(ϕ.dimensions)")
    end
    idx = get(mapping, state) do
        available = join(string.(keys(mapping)), ", ")
        error("State $(state) not available for $(node). Available: $(available)")
    end
    return float(ϕ.potential[idx])
end

function reactor_failure_probability(
    path::AbstractString;
    failure_state::Symbol=:Reactor_fail,
    age_state::Symbol=DEFAULT_AGE_STATE,
    distance_state::Symbol=DEFAULT_DISTANCE_STATE,
)
    data = load(path)
    ebn = haskey(data, "ebn") ? data["ebn"] : error("File $(path) does not contain an `ebn` object")
    bn = dispatch(ebn)
    evidence = Evidence(:AGE => age_state, :DISTANCE => distance_state)
    ϕ = infer(bn, :Reactor, evidence)
    return probability_of_state(ϕ, :Reactor, failure_state)
end

function plot_reactor_failure(
    entries;
    failure_state::Symbol=:Reactor_fail,
    age_state::Symbol=DEFAULT_AGE_STATE,
    distance_state::Symbol=DEFAULT_DISTANCE_STATE,
)
    mkpath(PLOT_DIR)
    labels = String[]
    probs = Float64[]
    for entry in entries
        try
            p = reactor_failure_probability(
                entry.path;
                failure_state=failure_state,
                age_state=age_state,
                distance_state=distance_state,
            )
            push!(labels, basename(entry.path))
            push!(probs, p)
            println("Reactor failure probability for $(basename(entry.path)) with AGE=$(age_state), DISTANCE=$(distance_state): $(p)")
        catch err
            @warn "Skipping reactor failure probability for $(entry.path)" error=err
        end
    end
    isempty(probs) && return
    plt = bar(
        labels,
        probs;
        xlabel = "Network",
        ylabel = "P(Reactor failure)",
        title = "Reactor failure probability comparison",
        legend = false,
        size = (900, 500),
    )
    fname = joinpath(PLOT_DIR, "reactor_failure_comparison.png")
    savefig(plt, fname)
    println("Saved reactor failure comparison plot: $(fname)")
end

function plot_cpt_differences(entries; atol=1e-10, rtol=1e-10, top_n::Int=15, threshold::Float64=0.0)
    mkpath(PLOT_DIR)
    ref = entries[1]
    ref_labels = row_labels(ref.table)
    for entry in entries[2:end]
        if entry.node_name != ref.node_name
            @warn "Skip plotting: node name differs" entry.path entry.node_name ref_node=ref.node_name
            continue
        end
        if names(entry.table) != names(ref.table)
            @warn "Skip plotting: CPT columns differ" entry.path
            continue
        end
        ref_probs = ref.table[:, :Π]
        cur_probs = entry.table[:, :Π]
        if length(ref_probs) != length(cur_probs)
            @warn "Skip plotting: CPT row counts differ" entry.path
            continue
        end
        diffs = cur_probs .- ref_probs
        keep = findall(abs.(diffs) .> threshold)
        isempty(keep) && begin
            println("No rows exceed threshold $(threshold) for $(basename(entry.path)); skipping plot.")
            continue
        end
        order = sortperm(abs.(diffs[keep]); rev = true)
        sel = keep[order[1:min(top_n, length(order))]]
        ylabels = label_slice.(ref_labels[sel])
        ys = 1:length(sel)
        plt = bar(
            diffs[sel],
            ys;
            orientation = :horizontal,
            yticks = (ys, ylabels),
            xlabel = "ΔΠ (current - reference)",
            ylabel = "CPT row (parents | state)",
            title = "Top $(length(sel)) CPT diffs: $(basename(entry.path)) vs $(basename(ref.path))",
            legend = false,
            size = (1000, 600),
        )
        vline!(plt, [0.0]; color=:black, lw=1, ls=:dash)
        fname = joinpath(PLOT_DIR, "cpt_diff_top_$(basename(entry.path)).png")
        savefig(plt, fname)
        println("Saved CPT difference plot: $(fname)")
        summary_df = DataFrame(
            row = sel,
            label = ref_labels[sel],
            delta = diffs[sel],
            ref_Pi = ref.table[sel, :Π],
            cur_Pi = entry.table[sel, :Π],
        )
        csv_name = joinpath(PLOT_DIR, "cpt_diff_top_$(basename(entry.path)).csv")
        CSV.write(csv_name, summary_df)
        println("Saved CPT difference table: $(csv_name)")
    end
end

function compare_last_node_cpts(
    dir::AbstractString=NETWORK_DIR;
    atol=1e-10,
    rtol=1e-10,
    max_report::Int=5,
    make_plots::Bool=true,
    top_n::Int=15,
    threshold::Float64=0.0,
    failure_state::Symbol=:Reactor_fail,
    age_state::Symbol=DEFAULT_AGE_STATE,
    distance_state::Symbol=DEFAULT_DISTANCE_STATE,
)
    paths = sort(filter(p -> endswith(p, ".jld2"), readdir(dir; join=true)))
    isempty(paths) && error("No .jld2 files found in $(dir)")
    entries = load_last_node_cpt.(paths)
    ref = entries[1]
    println("Reference file: $(basename(ref.path)) | last node: $(ref.node_name) :: $(ref.node_type)")
    for entry in entries[2:end]
        same_node = entry.node_name == ref.node_name
        equal_tables = same_node && tables_equal(entry.table, ref.table; atol=atol, rtol=rtol)
        println("Comparing $(basename(entry.path))")
        if !same_node
            println("  node differs: $(entry.node_name) vs $(ref.node_name)")
            continue
        end
        println("  CPT $(equal_tables ? "matches" : "differs") reference")
        equal_tables && continue
        diff_rows = mismatching_rows(ref.table, entry.table; atol=atol, rtol=rtol)
        for row in diff_rows[1:min(length(diff_rows), max_report)]
            println("  row $(row): reference=$(ref.table[row, :]) | current=$(entry.table[row, :])")
        end
        length(diff_rows) > max_report && println("  … $(length(diff_rows) - max_report) more differing rows omitted")
    end
    if make_plots
        plot_cpt_differences(entries; atol=atol, rtol=rtol, top_n=top_n, threshold=threshold)
        plot_reactor_failure(
            entries;
            failure_state=failure_state,
            age_state=age_state,
            distance_state=distance_state,
        )
    end
    return entries
end

compare_last_node_cpts(top_n=10, threshold=0.001, make_plots=true)
