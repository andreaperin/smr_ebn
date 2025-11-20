#!/bin/bash
#SBATCH --account=andrea.perin
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --output=slurm/mc_output.txt
#SBATCH --error=slurm/mc_error.txt

julia -p $SLURM_CPUS_PER_TASK --project=/home/andrea.perin/smr_ebn -e 'using Pkg; Pkg.instantiate(); include("./standard_mc.jl")'
