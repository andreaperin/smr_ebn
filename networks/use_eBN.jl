using EnhancedBayesianNetworks
using JLD2
using Plots
using Printf

const MODEL_PATH = "/Users/stefanomarchetti/Library/CloudStorage/OneDrive-PolitecnicodiMilano/Python/Cursor/smr_ebn/networks/ebn_jld2/2025_11_05_21_16_MonteCarlo(50).jld2"
const FAILURE_STATE = :Reactor_fail
const AGE_STATES = Tuple(Symbol("AGE_$(i)0") for i in 0:5)
const AGE_EVIDENCE_STATE = AGE_STATES[1]
const DISTANCE_STATE = :Distance_500
const PLOT_DIR = joinpath(@__DIR__, "imgs")
const FAILURE_SCALE = 1.2e-3
const HYDROGEN_ACCIDENTS = [
    (node = :OC, state = :YES_OC, label = "OC"),
    (node = :EXP, state = :YES_EXP, label = "EXP"),
    (node = :MSLBH2, state = :YES_MSLBH2, label = "MSLB_H2"),
]

@load MODEL_PATH ebn
bn = dispatch(ebn)

distance_evidence() = Evidence(:DISTANCE => DISTANCE_STATE, :AGE => AGE_EVIDENCE_STATE)
distance_evidence(pairs::Pair...) = Evidence(:DISTANCE => DISTANCE_STATE, :AGE => AGE_EVIDENCE_STATE, pairs...)

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
    base_evidence = distance_evidence()
    failure_factor = infer(bn, :Reactor, base_evidence)
    p_failure = probability_of_state(failure_factor, :Reactor, failure_state)
    failure_evidence = distance_evidence(:Reactor => failure_state)

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

        prior_factor = infer(bn, node, base_evidence)
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

ranked_impacts, p_failure = rank_hydrogen_accidents(bn, HYDROGEN_ACCIDENTS, FAILURE_STATE)

println("Hydrogen accident ranking via backward inference (evidence: Reactor = $(FAILURE_STATE))")
println("Additional evidence: DISTANCE = $(DISTANCE_STATE), AGE = $(AGE_EVIDENCE_STATE)")
scaled_p_failure = p_failure * FAILURE_SCALE
@printf "Baseline reactor failure probability P(F): %.6e\n\n" scaled_p_failure
@printf "%-10s %12s %16s %16s %16s\n" "Accident" "P(A)" "P(A|F)" "P(F|A)" "ΔP(F)"
separator = repeat("-", 74)
println(separator)
for r in ranked_impacts
    failure_given_accident = isfinite(r.failure_given_accident) ? r.failure_given_accident : NaN
    uplift = isfinite(r.uplift) ? r.uplift : NaN
    scaled_failure_given = isfinite(failure_given_accident) ? failure_given_accident * FAILURE_SCALE : NaN
    scaled_uplift = isfinite(uplift) ? uplift * FAILURE_SCALE : NaN
    @printf "%-10s %12.6e %16.6e %16.6e %16.6e\n" r.label r.prior r.posterior scaled_failure_given scaled_uplift
end

function years_from_age_state(state::Symbol)
    suffix = split(String(state), "_")[end]
    return parse(Int, suffix)
end

function failure_probability_vs_age(bn::BayesianNetwork, failure_state::Symbol; states=AGE_STATES)
    results = NamedTuple{(:state, :years, :probability),Tuple{Symbol,Int,Float64}}[]
    for state in states
        ev = distance_evidence(:AGE => state)
        factor = infer(bn, :Reactor, ev)
        push!(results, (
            state = state,
            years = years_from_age_state(state),
            probability = float(probability_of_state(factor, :Reactor, failure_state)) * FAILURE_SCALE,
        ))
    end
    sort(results; by = r -> r.years)
end

age_failure = failure_probability_vs_age(bn, FAILURE_STATE)

println("\nReactor failure probability conditioned on AGE states")
@printf "%-8s %12s\n" "AGE" "P(F|AGE)"
println(repeat("-", 24))
for entry in age_failure
    @printf "%-8s %12.6e\n" String(entry.state) entry.probability
end

mkpath(PLOT_DIR)
times = [r.years for r in age_failure]
probs = [r.probability for r in age_failure]
plt = plot(
    times,
    probs,
    marker = :circle,
    xlabel = "Plant age [years]",
    ylabel = "P(Reactor failure | AGE)",
    title = "Failure probability vs. plant age",
    legend = false,
)
plot_path = joinpath(PLOT_DIR, "reactor_failure_vs_age.png")
savefig(plt, plot_path)
println("\nSaved failure probability plot to $(plot_path)")
