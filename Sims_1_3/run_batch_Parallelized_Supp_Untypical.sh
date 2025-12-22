#!/bin/bash

#SBATCH --account=SSCM013902 #EDIT
#SBATCH --job-name NCC
#SBATCH -o /data/Epic/subprojects/BiocratesGWAS/work/WeightedNCC_Final/Simulations/shlog/%A_%a.log #EDIT
#SBATCH -e /data/Epic/subprojects/BiocratesGWAS/work/WeightedNCC_Final/Simulations/shlog/%A_%a.err #EDIT
#SBATCH --ntasks=1                     # EDIT One task per job - pointless/wasteful to increase (1. mclapply runs the parallelization on the same node 2. if using a slurm array, 1 ntasks will be automatically allocated to each job in the array vector)
#SBATCH --cpus-per-task=16              # EDIT CPUs per task (set the min (num_simulations, max(available cpus_per_node)))
#SBATCH --mem-per-cpu=1500
#SBATCH --array=0-3 # EDIT set len(array) as len(parameter combo grid)
#SBATCH -p low_p

# Load environment variables from config.env EDIT loads workdir used in Rscript
if [ -f config.env ]; then
    source config.env
else
    echo "Error: config.env file not found!"
    exit 1
fi

module load languages/R/4.4.1
Rscript Sim_Supp_Untypical.R "$SLURM_ARRAY_TASK_ID" "$WORKDIR"
