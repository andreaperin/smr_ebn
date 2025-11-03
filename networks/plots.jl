using EnhancedBayesianNetworks
using JLD2
using Cairo
using PGFPlotsX

mkpath("networks/imgs")  # creates the folder if it doesn’t exist

@load "/Users/stefanomarchetti/Library/CloudStorage/OneDrive-PolitecnicodiMilano/Python/Cursor/smr_ebn/networks/ebn_jld2/2025_10_29_17_04_MonteCarlo(50).jld2"

fig1 = gplot(ebn, NODESIZEFACTOR=0.1, NODELABELSIZE=6)
draw(PNG("networks/imgs/1_ebn.png", 16cm, 16cm), fig1)
