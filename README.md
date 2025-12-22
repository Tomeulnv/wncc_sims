# Replication Package for Simulation Studies  

This repository contains the code required to **replicate all simulation results and figures** presented in the paper:
**“On the Estimation of Inclusion Probabilities for Weighted Analyses of Nested Case-Control Studies”**

#### DISCLAIMER
This replication package is undergoing improvements to facilitate smoother replication.

## Directory Structure

```
.
├── README.md
├── Sims_1_3
│   ├── analyse_results.R
│   ├── analyse_results_SuppMat_AssoExpo.R
│   ├── analyse_results_SuppMat_Untypical.R
│   ├── config.env
│   ├── figure1.R
│   ├── Figures                                             # Folder containing generated figures
│   ├── Res                                                 # Folder containing simulation outputs and result datasets
│   ├── RES                                                 # Folder containing datasets used to produce Figure 1
│   ├── RFunctions                                          # Folder containing auxiliary R functions
│   ├── run_batch_Parallelized_Sim_1_3.sh                   # Bash script to run simulations for Sections 3.2 and 3.4
│   ├── run_batch_Parallelized_Supp_AssoExpo.sh             # Bash script to run simulations for the Supplementary Material (association/exposure)
│   ├── run_batch_Parallelized_Supp_Untypical.sh            # Bash script to run simulations for the Supplementary Material (untypical NCC)
│   ├── shlog                                               # Folder containing SLURM log files
│   ├── Sim_1_3.R
│   ├── Sim_Supp_AssoExpo_Inter.R
│   ├── Sim_Supp_AssoExpo.R
│   └── Sim_Supp_Untypical.R
└── Sims_2
    ├── display_results.R
    ├── InterMatching_Par.R
    ├── Res                                                 # Folder containing simulation outputs and result datasets
    ├── RFunctions                                          # Folder containing auxiliary R functions
    └── run_batch_par.sh                                    # Bash script to run simulations for Section 3.3

```

## Replication of the Results from the Simulations

### Figure 1: Motivating Example
Run the script `/Sims_1_3/figure1.R` to produce figure 1.

### Section 3.2: NN vs. Caliper Matching and Section 3.4: When Matching Factors May Be Excluded from Weight Computation
The replication of the results for sections 3.2 and 3.4 can be found in the folder Sims_1_3 by following these steps:

1. Run the script `run_batch_Parallelized_Sim_1_3.sh`. 

    (*Note: this script is designed to run on a SLURM-based system.*)

2. Run the script: `analyse_results.R` to obtain the figures included in the paper

###  Section 3.3: Interaction Between Matching Factors and Disease Risk
For replication of the results of section 3.3, run the script `run_batch_par.sh`.
Results are outputted into the folder `Sims_2/Res`

### Supplementary Material
- Scripts `/Sims_1_3/run_batch_Parallelized_Supp_AssoExpo.sh`and `/Sims_1_3/run_batch_Parallelized_Supp_Untypical.sh` produce the simulations and
results presented in the supplementary material. 

- To produce the figures run `analyse_results_SuppMat_AssoExpo.R` and `analyse_results_SuppMat_Untypical.R`