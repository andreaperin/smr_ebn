using EnhancedBayesianNetworks
using Cairo
using PGFPlotsX

@load "networks/ebn_jld2/2025_09_15_09_19_MonteCarlo(50).jld2"

fig1 = gplot(ebn, NODESIZEFACTOR=0.2, NODELABELSIZE=6)
draw(PNG("networks/imgs/1_ebn.png", 16cm, 16cm), fig1)
