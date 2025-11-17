using EnhancedBayesianNetworks
using JLD2
using Printf
using Plots

const MODEL_PATH = "/Users/stefanomarchetti/Library/CloudStorage/OneDrive-PolitecnicodiMilano/Python/Cursor/smr_ebn/networks/ebn_jld2/2025_11_05_21_16_MonteCarlo(50).jld2"
const FAILURE_STATE = :Reactor_fail
const PLOT_PATH = "networks/imgs/hydrogen_accident_ranking.png"
const HYDROGEN_ACCIDENTS = [
    (node = :OC, state = :YES_OC, label = "OC"),
    (node = :EXP, state = :YES_EXP, label = "EXP"),
    (node = :MSLBH2, state = :YES_MSLBH2, label = "MSLB_H2"),
]

@load MODEL_PATH ebn
bn = dispatch(ebn)

function probability_of_state(ϕ::Factor, node::Symbol, state::Symbol)
    mapping = get(ϕ.states_mapping, node) do
        error("Node $(node) not found in factor dimensions $(ϕ.dimensions)")
    end
    idx = get(mapping, state) do
        available = join(string.(keys(mapping)), ", ")
        error("State $(state) not available for $(node). Available: $(available)")
    end
    return ϕ.potential[idx]
end

function rank_hydrogen_accidents(bn::BayesianNetwork, accidents, failure_state::Symbol)
    failure_factor = infer(bn, :Reactor)
    p_failure = probability_of_state(failure_factor, :Reactor, failure_state)
    failure_evidence = Evidence(:Reactor => failure_state)

    results = Vector{NamedTuple{(
        :label,
        :node,
        :prior,
        :posterior,
        :failure_given_accident,
        :uplift
    ),Tuple{String,Symbol,Float64,Float64,Float64,Float64}}}()
    for spec in accidents
        node = spec.node
        accident_state = spec.state
        label = spec.label

        prior_factor = infer(bn, node)
        posterior_factor = infer(bn, node, failure_evidence)

        p_accident = float(probability_of_state(prior_factor, node, accident_state))
        p_accident_given_failure = float(probability_of_state(posterior_factor, node, accident_state))

        impact = isnothing(p_accident) || isapprox(p_accident, 0.0; atol=eps(1.0)) ?
            NaN :
            p_accident_given_failure * p_failure / p_accident

        push!(results, (
            label = label,
            node = node,
            prior = p_accident,
            posterior = p_accident_given_failure,
            failure_given_accident = impact,
            uplift = impact - p_failure,
        ))
    end

    order = sortperm(results; by = r -> isfinite(r.failure_given_accident) ? r.failure_given_accident : -Inf, rev = true)
    return results[order], p_failure
end

function plot_accident_impacts(ranked_impacts, failure_state::Symbol; output_path::String=PLOT_PATH)
    isempty(ranked_impacts) && return nothing
    labels = [r.label for r in ranked_impacts]
    raw_values = [r.failure_given_accident for r in ranked_impacts]
    values = [isfinite(v) ? v : 0.0 for v in raw_values]
    xs = collect(1:length(labels))

    max_val = maximum(values)
    ymax = iszero(max_val) ? 1.0 : max_val * 1.1

    bar_kwargs = (
        xlabel = "Hydrogen accident",
        ylabel = "P(" * string(failure_state) * "|accident)",
        label = "",
        color = :royalblue4,
        grid = :y,
        legend = false,
        bar_width = 0.55,
        dpi = 300,
        size = (900, 600),
        ylims = (0, ymax),
        xticks = (xs, labels),
        xtickfont = Plots.font(10, rotation = 20),
    )

    plt = bar(xs, values; bar_kwargs...)

    mkpath(dirname(output_path))
    savefig(plt, output_path)
    return output_path
end

ranked_impacts, p_failure = rank_hydrogen_accidents(bn, HYDROGEN_ACCIDENTS, FAILURE_STATE)

println("Hydrogen accident ranking via backward inference (evidence: Reactor = $(FAILURE_STATE))")
@printf "Baseline reactor failure probability P(F): %.6e\n\n" p_failure
@printf "%-10s %12s %16s %16s %16s\n" "Accident" "P(A)" "P(A|F)" "P(F|A)" "ΔP(F)"
separator = repeat("-", 74)
println(separator)
for r in ranked_impacts
    failure_given_accident = isfinite(r.failure_given_accident) ? r.failure_given_accident : NaN
    uplift = isfinite(r.uplift) ? r.uplift : NaN
    @printf "%-10s %12.6e %16.6e %16.6e %16.6e\n" r.label r.prior r.posterior failure_given_accident uplift
end

plot_path = plot_accident_impacts(ranked_impacts, FAILURE_STATE; output_path=PLOT_PATH)
if plot_path === nothing
    println("No ranking results available for plotting.")
else
    println("Saved ranking plot to $(plot_path)")
end
