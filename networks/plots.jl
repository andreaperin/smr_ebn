using EnhancedBayesianNetworks
using Cairo
using PGFPlotsX

mkpath("networks/imgs")  # creates the folder if it doesn’t exist

@load "/Users/stefanomarchetti/Library/CloudStorage/OneDrive-PolitecnicodiMilano/Python/Cursor/smr_ebn/networks/ebn_jld2/2025_09_23_16_33_MonteCarlo(150).jld2"

fig1 = gplot(ebn, NODESIZEFACTOR=0.2, NODELABELSIZE=6)
draw(PNG("networks/imgs/1_ebn.png", 16cm, 16cm), fig1)
