# pathvv <- "/data/Epic/subprojects/BiocratesGWAS/work/WeightedNCC_Final/Simulations/"
# setwd(pathvv)


# Emtpy R working space

rm(list=ls())

#Libraries
lop <- c("parallel","tidyverse","survival","survminer","ggplot2", "multipleNCC", "mgcv", "pbapply")
newp <- lop[!(lop %in% installed.packages()[,"Package"])]
if(length(newp)) install.packages(newp,  repos="https://cloud.r-project.org/")
lapply(lop, require, character.only = TRUE)


#setwd("/user/home/nu23752/work/wNCC_R")

################################################################################
# Initialisation
################################################################################

# Set number of simulations and cores to use
NSIM <- 100
RNGkind("L'Ecuyer-CMRG") #reproducibility

# Fetch relevant SLURM environment variables
numNodes <- as.numeric(Sys.getenv("SLURM_NNODES", unset=1))  # Total nodes
cat("Total number of nodes:", numNodes, "\n")
numTasks <- as.numeric(Sys.getenv("SLURM_NTASKS", unset=NA))  # Total tasks
cat("Total number of tasks:", numTasks, "\n")
tasksPerNode <- as.numeric(Sys.getenv("SLURM_NTASKS_PER_NODE", unset=NA))  # Tasks per node
cat("Total number of tasks per node:", tasksPerNode, "\n")
cpusPerTask <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset=1))  # CPUs per task
cat("Total number of cpus per task:", cpusPerTask, "\n")

# Compute the total number of workers
if (!is.na(numTasks)) {
  numWorkers <- numTasks * cpusPerTask
} else if (!is.na(tasksPerNode)) {
  numWorkers <- numNodes * tasksPerNode * cpusPerTask
} else {
  numWorkers <- parallel::detectCores(logical = FALSE) - 1  # Fallback: use all physical cores
}

numWorkers <- ifelse(!is.na(numTasks), numTasks * cpusPerTask, parallel::detectCores(logical = FALSE) - 1)

cat("Total number of workers:", numWorkers, "\n")

args <- commandArgs(trailingOnly = TRUE)
cat("args", args, "\n")

#ind <- as.integer(args[1]) + 1  # SLURM_ARRAY_TASK_ID is 0-based, so we add 1 to match the row index
ind<-as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))+1
cat("ind", ind, "\n")


workdir<-args[2]
cat("wd", workdir, "\n")

setwd(workdir)


################################################################################
# Initialisation
################################################################################
# Set WD and Import the required functions
source("RFunctions/Functions_Simulations_T.R")
source("RFunctions/Functions_Computation_Weights_T.R")
source("RFunctions/Functions_Survival_CondSurvival_Estimates_T.R")

# Set simulation parameters
ALLR2SNPBMI           <- 0.1 # c(0, 10)/100
ALLlogHRBMI           <- log(c(1, 2))
ALLlogHRSNP           <- log(c(1, 2))
ALLlogHRMatching      <- log(c(2))
# ALLlogHRMatchingSNP   <- log(c(1, 2))
ALLMAF                <- c(0.25)
PSELCASES             <- c(1)
AgeMinCens            <- c(0, 20)
AgeMaxCens            <- c(30, 50)
NNMatching            <- c(0)

# Define all combinations of parameters
ALLPRMTR0            <- expand.grid(ALLR2SNPBMI, ALLlogHRBMI, ALLlogHRSNP, log(2), ALLMAF, PSELCASES, AgeMinCens, AgeMaxCens, NNMatching, 1) 
                             
colnames(ALLPRMTR0)  <- c("R2SNPBMI", "logHRBMI", "logHRSNP", "logHRMatching", "MAF", "PSELCASES", 
                                                           "AgeMinCens", "AgeMaxCens", "NNmatching", "Asso_M_SNP")
ALLPRMTR             <- ALLPRMTR0[!(((ALLPRMTR0[, "AgeMinCens"]==0) & (ALLPRMTR0[, "AgeMaxCens"]==50))|
                                         ((ALLPRMTR0[, "logHRBMI"]==0) & (ALLPRMTR0[, "logHRSNP"]==0))), ]


CohortSize    <- 1e5
logHRSNP      <- ALLPRMTR$logHRSNP[ind]
logHRBMI      <- ALLPRMTR$logHRBMI[ind]
logHRMatching <- ALLPRMTR$logHRMatching[ind]
R2SNPBMI      <- ALLPRMTR$R2SNPBMI[ind]
MAF           <- ALLPRMTR$MAF[ind]
pselCases     <- ALLPRMTR$PSELCASES[ind]
agemincens    <- ALLPRMTR$AgeMinCens[ind]
agemaxcens    <- ALLPRMTR$AgeMaxCens[ind]
nnmatching    <- ALLPRMTR$NNmatching[ind]
assoMsnp_num  <- ALLPRMTR$Asso_M_SNP[ind]
assoMsnp      <- ifelse(assoMsnp_num==1, TRUE, FALSE)

DirError      <- "Res/"
suffix <- paste0("Sim_Supp_Asso_CohortSize_", CohortSize, "_MAF_", MAF, "_R2SNPBMI_", R2SNPBMI, "_logHRSNP_", round(logHRSNP,2),
                 "_logHRBMI_", round(logHRBMI,2), "_logHRMatching_", round(logHRMatching,2), 
                 "_pselCases_", pselCases, "_ageCens_", agemincens, "_", agemaxcens, 
                 "_NNmatching_", nnmatching, "_AssoMsnp_", assoMsnp_num)

RECAP_RES_ASSO <- NULL
Counts <- NULL

set.seed(1234)


tryCatch({
  tempcohort_sim1 <- simul_cohort_univM(n=CohortSize, 
                                        R2SNPBMI, logHRSNP, logHRBMI, logHRMatching, 
                                        MeanBMI=NA, SDBMI=NA, MeanSNP=NA, SDSNP=NA, 
                                        MAF=MAF, agemincens, agemaxcens, assoMsnp)
  
  
  cat("tempcohort cols:", colnames(tempcohort_sim1$MeanSNP), "\n")
  
  
  MeanBMI          <- tempcohort_sim1$MeanBMI
  SDBMI            <- tempcohort_sim1$SDBMI
  MeanSNP          <- tempcohort_sim1$MeanSNP
  SDSNP            <- tempcohort_sim1$SDSNP
  
  
  cat("tempcohort_sim1:", DirError, "\n")
  
}, error = function(e) {
  writeLines(as.character(e), paste0(DirError, "ErrorSimCohort_", ind, "_", 1, ".txt"))
})


internal_task <- function(nsim) {
  
  set.seed(nsim)
  cat("Starting simulation", nsim, "\n")
  tryCatch(
    {
      if (nsim == 1)
      {
        cat("tempcohort_sim1:", DirError, "\n")
        tempcohort <- tempcohort_sim1
      }else
      {
        tempcohort <- simul_cohort_univM(n=CohortSize, 
                                         R2SNPBMI, logHRSNP, logHRBMI, logHRMatching, 
                                         MeanBMI=NA, SDBMI=NA, MeanSNP=NA, SDSNP=NA, 
                                         MAF=MAF, agemincens, agemaxcens, assoMsnp)
      }
    },
    error  = function(e){writeLines(as.character(e), paste0(DirError, "ErrorSimCohort_", ind, "_", nsim, ".txt"))}
  )
  
  # Sample a matched case-control study using density sampling
  tryCatch(
    {
      if (nnmatching==1) 
      {
        forNCC  <- sample_case_control_univM_sample_NNMatching(cohort = tempcohort$CohortData, pselCases = pselCases)
        NCC     <- forNCC$NCC
      }else{
        forNCC  <- sample_case_control_univM_sample_CaliperMatching(cohort = tempcohort$CohortData, pselCases = pselCases, returnRiskSet=F)
        NCC     <- forNCC$NCC  
      }
    }, 
    error              = function(e){writeLines(as.character(e), paste0(DirError, "ErrorNCC_", ind, "_", nsim, ".txt"))}
  )
  # Compute the Weights
  
  ### can we have something based on your function (compute_weights) rather than mine?
  ### and also run it "in parallel" for NCC_NN and NCC_Caliper?
  tryCatch({
    cohort0             <- prepare_for_weights(tempcohort$CohortData, NCC, univar=T)
    cohort_withweights  <- compute_weights(cohort = cohort0, NCC = NCC,
                                           pselCases = pselCases,
                                           # ListRiskSets_Matching = forNCC$ListRiskSets_Matching, 
                                           univariate = T, 
                                           nnmatch = (nnmatching==1), 
                                           withTDontop = T)
  }, 
  error = function(e){writeLines(as.character(e), paste0(DirError, "ErrorWeights_", ind, "_", nsim, ".txt"))
  })
  
  rm(forNCC)
  if (nsim == 1) {
    listres <- list(cohort_withweights = cohort_withweights, NCC = NCC)
    cat("listres", DirError, "\n")
    saveRDS(listres, file = paste0("Res/Data_test_", suffix, ".rds"))
  }
  tryCatch({
    
    
    #### LINEAR ASSOCIATION BETWEEN THE SNP AND BMI
    mod_allcohort                             <- lm(BMI ~ SNP, data= cohort_withweights)
    mod_casectrlWeightedKM_matching           <- lm(BMI ~ SNP, weights= WeightsKM_Matching, data= cohort_withweights)
    mod_casectrlWeightedGAM_matching          <- lm(BMI ~ SNP, weights= WeightsGAM_Matching, data= cohort_withweights)
    mod_casectrlWeightedGAM_matching_thresh   <- lm(BMI ~ SNP, weights= WeightsGAM_Matching_thresh, data= cohort_withweights)
    mod_casectrlWeightedGAM_noTD              <- lm(BMI ~ SNP, weights= WeightsGAM_noTD, data= cohort_withweights)
    mod_NCC_ctrl                              <- lm(BMI ~ SNP, data= subset(NCC, Ind==0))
    mod_NCC                                   <- lm(BMI ~ SNP, data= NCC %>% distinct(Id, .keep_all = T))
    mod_subsamplecohort                       <- lm(BMI ~ SNP, data= cohort_withweights[sample(1:nrow(cohort_withweights), nrow(NCC), replace=F), ])
  },
  error = function(e){writeLines(as.character(e), paste0(DirError, "ErrorEstimAssos_", ind, "_", nsim, ".txt"))
  })
  cat("Return:", nsim, "\n")
  
  tryCatch(
    {
      monster <-  c(nsim, pselCases, MAF, R2SNPBMI, logHRSNP, logHRBMI, logHRMatching, 
                    agemincens, agemaxcens, assoMsnp_num, 
                    coef(mod_allcohort)[2],
                    coef(mod_subsamplecohort)[2],
                    coef(mod_NCC_ctrl)[2],
                    coef(mod_NCC)[2],
                    coef(mod_casectrlWeightedKM_matching)[2],
                    coef(mod_casectrlWeightedGAM_matching)[2],
                    coef(mod_casectrlWeightedGAM_matching_thresh)[2],
                    coef(mod_casectrlWeightedGAM_noTD)[2],
                    sum(cohort_withweights$Ind),
                    sum(cohort_withweights$WeightsGAM_Matching > 100),
                    sum(cohort_withweights$WeightsGAM_Matching > 1000))
      
      cat("Finished simulation", monster, "\n")
      return(                           monster)
      
    }, 
    error              = function(e){writeLines(as.character(e), paste0(DirError, "ErrorReturn_", ind, "_", nsim, ".txt"))}
  )
  
}

RECAP_RES_ASSO <- do.call(rbind, mclapply(1:NSIM, internal_task, mc.cores = max(1, numWorkers - 1), mc.preschedule = FALSE))
gc()
cat("length cols", ncol(RECAP_RES_ASSO), "\n")
cat("length rows", nrow(RECAP_RES_ASSO), "\n")
allnames <- c("nsim", "pselCases","MAF","R2SNPBMI", "logHRSNP", "logHRBMI", "logHRM", "Agemincens", "Agemaxcens", "AssoMsnp",
              "beta_allcohort", "beta_allcohort_subsample",
              "beta_NCC_ctrl",
              "beta_NCC",
              "beta_NCC_WeightedKM_matching",
              "beta_NCC_WeightedGAM_matching",
              "beta_NCC_WeightedGAM_matching_thresh",
              "beta_NCC_WeightedGAM_noTD",
              "nbCases",
              "nbGAMWeightsMatchingabove100",
              "nbGAMWeightsMatchingabove1000")
cat("length", length(allnames), "\n")

colnames(RECAP_RES_ASSO) <- allnames
cat("length", length(colnames(RECAP_RES_ASSO)), "\n")

AllRes <- list(RECAP_RES_ASSO = RECAP_RES_ASSO, Counts = Counts)
saveRDS(AllRes, file = paste0("Res/Res_", suffix, ".rds"))
