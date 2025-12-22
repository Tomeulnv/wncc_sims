rm(list = ls())
library(parallel)
library(tidyverse)
library(survival)
library(ggplot2)
library(survminer)
library(multipleNCC)
library(mgcv)
library(pbapply)

################################################################################
# Initialisation
################################################################################
# Set WD and Import the required functions
setwd(".")
source("RFunctions/Functions_Simulations.R")
source("RFunctions/Functions_Computation_Weights.R")
source("RFunctions/Functions_Survival_CondSurvival_Estimates.R")

# Set number of simulations and cores to use
nsims       <- 100
numWorkers  <- 25

# Set simulation parameters
R2SNPBMI              <- 0.1
ALLlogHRBMI           <- log(c(2))
ALLlogHRSNP           <- log(c(2))
ALLlogHRM             <- log(c(1, 2))
ALLlogHRM1M2          <- log(c(1, 2))
ALLlogHRSNP_M1M2      <- log(c(1, 2))
ALLMAF                <- c(0.25)
Asso_M_SNP            <- 0.2
PSELCASES             <- c(1)
AgeMinCens            <- c(20)
AgeMaxCens            <- c(50)
NNmatching            <- FALSE
matchMethodLabel      <- if (NNmatching) "NN" else "Caliper"

# Define all combinations of parameters
ALLPRMTR              <- expand.grid(R2SNPBMI, ALLlogHRBMI, ALLlogHRSNP, ALLlogHRM, ALLlogHRM1M2, ALLlogHRSNP_M1M2, ALLMAF, PSELCASES, AgeMinCens, AgeMaxCens)
colnames(ALLPRMTR)    <- c("R2SNPBMI", "logHRBMI", "logHRSNP", "logHRM", "logHRM1M2", "logHRSNP_M1M2", "MAF", "PSELCASES", "AgeMinCens", "AgeMaxCens")

################################################################################
# Function for computation in parallel
################################################################################
fun_all_para <- function(ind) {
  logHRSNP      <- ALLPRMTR$logHRSNP[ind]
  logHRBMI      <- ALLPRMTR$logHRBMI[ind]
  logHRM        <- ALLPRMTR$logHRM[ind]
  logHRM1M2     <- ALLPRMTR$logHRM1M2[ind]
  logHRSNP_M1M2 <- ALLPRMTR$logHRSNP_M1M2[ind]
  R2SNPBMI      <- ALLPRMTR$R2SNPBMI[ind]
  MAF           <- ALLPRMTR$MAF[ind]
  pselCases     <- ALLPRMTR$PSELCASES[ind]
  agemincens    <- ALLPRMTR$AgeMinCens[ind]
  agemaxcens    <- ALLPRMTR$AgeMaxCens[ind]
  DirError      <- "Res/"
  suffix <- paste0(
    "InterMatching_onSurvival_MAF_", MAF, "_Asso_M_SNP_", Asso_M_SNP,
    "_R2SNPBMI_", R2SNPBMI, "_logHRSNP_", round(logHRSNP, 2),
    "_logHRBMI_", round(logHRBMI, 2), "_logHRM_",
    round(logHRM, 2), "_logHRM1M2_", round(logHRM1M2, 2),
    "_logHRSNP_M1M2_", round(logHRSNP_M1M2, 2),
    "_pselCases_", pselCases,
    "_ageCens_", agemincens, "_", agemaxcens,
    "_", matchMethodLabel
  )

  RECAP_RES_ASSO <- NULL

  for (nsim in 1:nsims) {
    set.seed(nsim)
    cat("Starting simulation", nsim, "\n")
    tryCatch(
      {
        if (nsim == 1) {
          # First simulation
          tempcohort <- simul_cohort(n = 1e5, R2SNPBMI, logHRSNP, logHRBMI,
                                     logHRM, logHRM1M2, logHRSNP_M1M2,
                                     MeanBMI = NA, SDBMI = NA,
                                     MeanSNP = NA, SDSNP = NA,
                                     MeanM1 = NA, SDM1 = NA,
                                     MeanM2 = NA, SDM2 = NA,
                                     MeanM1M2 = NA, SDM1M2 = NA,
                                     MeanSNP_M1M2 = NA, SDSNP_M1M2 = NA,
                                     Asso_M_SNP = Asso_M_SNP,
                                     MAF = MAF, agemincens, agemaxcens)

          MeanBMI          <- tempcohort$MeanBMI
          SDBMI            <- tempcohort$SDBMI
          MeanSNP          <- tempcohort$MeanSNP
          SDSNP            <- tempcohort$SDSNP
          MeanM1           <- tempcohort$MeanM1
          SDM1             <- tempcohort$SDM1
          MeanM2           <- tempcohort$MeanM2
          SDM2             <- tempcohort$SDM2
          MeanM1M2         <- tempcohort$MeanM1M2
          SDM1M2           <- tempcohort$SDM1M2
          MeanSNP_M1M2     <- tempcohort$MeanSNP_M1M2
          SDSNP_M1M2       <- tempcohort$SDSNP_M1M2
          # Second simulation onwards:
        } else {

          tempcohort <- simul_cohort(n = 1e5, R2SNPBMI, logHRSNP, logHRBMI,
                                     logHRM, logHRM1M2, logHRSNP_M1M2,
                                     MeanBMI, SDBMI, MeanSNP, SDSNP,
                                     MeanM1, SDM1, MeanM2, SDM2,
                                     MeanM1M2, SDM1M2,
                                     MeanSNP_M1M2, SDSNP_M1M2,
                                     Asso_M_SNP = Asso_M_SNP,
                                     MAF, agemincens, agemaxcens)
        }
      Counts <- tempcohort$Counts
      },
      error  = function(e){writeLines(as.character(e), paste0(DirError, "ErrorSimCohort_", ind, "_", nsim, ".txt"))}
    )

    # Sample a matched case-control study using density sampling
    tryCatch(
      {
        forNCC <- if (NNmatching) {
          sample_case_control_NN(cohort = tempcohort$CohortData, pselCases = pselCases)
        } else {
          sample_case_control(cohort = tempcohort$CohortData, pselCases = pselCases, returnRiskSet = TRUE)
        }

        NCC <- forNCC$NCC
      },
      error = function(e) {
        writeLines(as.character(e), paste0(DirError, "ErrorNCC_", ind, "_", nsim, ".txt"))
      }
    )
    # Compute the Weights
    tryCatch({
      cohort0             <- prepare_for_weights(tempcohort$CohortData, NCC)
      cohort_withweights  <- compute_weights(cohort = cohort0, NCC = NCC,
                                             pselCases = pselCases,
                                             ListRiskSets_Matching = forNCC$ListRiskSets_Matching)
    },
      error = function(e){writeLines(as.character(e), paste0(DirError, "ErrorWeights_", ind, "_", nsim, ".txt"))
      })

    rm(forNCC)

    if (nsim == 1) {
      listres <- list(cohort_withweights = cohort_withweights, NCC = NCC)
      saveRDS(listres, file = paste0("Res/Data_test_", suffix, ".rds"))
    }
    tryCatch({
      cohort_withweights  <- cohort_withweights %>%
        select(-any_of(c(
                         "WeightsGAM_noMatching_thresh",
                         "WeightsGAM_Matching_thresh"))
        )

      #### LINEAR ASSOCIATION BETWEEN THE SNP AND BMI ####
      mod_allcohort                           <- lm(BMI ~ SNP, data= cohort_withweights)
      mod_casectrlWeightedKM_matching         <- lm(BMI ~ SNP, weights= WeightsKM_Matching, data= cohort_withweights)
      mod_casectrlWeightedGAM_nomatching      <- lm(BMI ~ SNP, weights= WeightsGAM_noMatching, data= cohort_withweights)
      mod_casectrlWeightedGAM_matching        <- lm(BMI ~ SNP, weights= WeightsGAM_Matching, data= cohort_withweights)
      mod_casectrlWeightedGAM_matchingInter   <- lm(BMI ~ SNP, weights= WeightsGAM_MatchingInter, data= cohort_withweights)
      mod_NCC_ctrl                            <- lm(BMI ~ SNP, data= subset(NCC, Ind==0))
      mod_NCC                                 <- lm(BMI ~ SNP, data= NCC %>% distinct(Id, .keep_all = T))
      mod_subsamplecohort                     <- lm(BMI ~ SNP, data= cohort_withweights[sample(1:nrow(cohort_withweights), nrow(NCC), replace=F), ])

      ##### COX PROPORTIONAL HAZARD RATE ESTIMATES ####
      modCox_allcohort                         <- coxph(Surv(AgeOut, Ind) ~ SNP + BMI + M1*M2 + M1:M2:SNP, data = cohort_withweights)
      modCox_casectrlWeightedKM_matching       <- coxph(Surv(AgeOut, Ind) ~ SNP + BMI + M1*M2 + M1:M2:SNP, weights= WeightsKM_Matching, data = subset(cohort_withweights, statforweight_02>0))
      modCox_casectrlWeightedGAM_nomatching    <- coxph(Surv(AgeOut, Ind) ~ SNP + BMI + M1*M2 + M1:M2:SNP, weights= WeightsGAM_noMatching, data = subset(cohort_withweights, statforweight_02>0))
      modCox_casectrlWeightedGAM_matching      <- coxph(Surv(AgeOut, Ind) ~ SNP + BMI + M1*M2 + M1:M2:SNP, weights= WeightsGAM_Matching, data = subset(cohort_withweights, statforweight_02>0))
      modCox_casectrlWeightedGAM_matchingInter <- coxph(Surv(AgeOut, Ind) ~ SNP + BMI + M1*M2 + M1:M2:SNP, weights= WeightsGAM_MatchingInter, data = subset(cohort_withweights, statforweight_02>0))
      mod_clogit                               <- clogit(IndCaseCtrlNCC   ~ SNP + BMI + M1*M2 + M1:M2:SNP + strata(Pair), data = NCC)

      ##### KM ESTIMATES ####
      km_est          <- km_estim(cohort_withweights, suffix_file_plotkm = NULL)
      KM_FullCohort                             <- 1 - summary(km_est$fitKM[1], times = agemaxcens, extend=T)$surv
      KM_casectrlWeightedKM_matching            <- 1 - summary(km_est$fitKM[5], times = agemaxcens, extend=T)$surv
      KM_casectrlWeightedGAM_nomatching         <- 1 - summary(km_est$fitKM[4], times = agemaxcens, extend=T)$surv
      KM_casectrlWeightedGAM_matching           <- 1 - summary(km_est$fitKM[2], times = agemaxcens, extend=T)$surv
      KM_casectrlWeightedGAM_matchingInter      <- 1 - summary(km_est$fitKM[3], times = agemaxcens, extend=T)$surv
      km_cond_est      <- km_cond_estim(cohort_withweights, suffix_file_plotkm = NULL)
      KM_cond_FullCohort                        <- 1 - summary(km_cond_est$fitKM[1], times = agemaxcens, extend=T)$surv
      KM_cond_casectrlWeightedKM_matching       <- 1 - summary(km_cond_est$fitKM[5], times = agemaxcens, extend=T)$surv
      KM_cond_casectrlWeightedGAM_nomatching    <- 1 - summary(km_cond_est$fitKM[4], times = agemaxcens, extend=T)$surv
      KM_cond_casectrlWeightedGAM_matching      <- 1 - summary(km_cond_est$fitKM[2], times = agemaxcens, extend=T)$surv
      KM_cond_casectrlWeightedGAM_matchingInter <- 1 - summary(km_cond_est$fitKM[3], times = agemaxcens, extend=T)$surv
    },
    error = function(e){writeLines(
      as.character(e),
      paste0(DirError, "ErrorEstimAssos_", ind, "_", nsim, ".txt"))
    })

    ##### STORE AND APPEND RESULTS
    RECAP_RES_ASSO <- rbind(RECAP_RES_ASSO,
      c(
        nsim, pselCases, MAF, R2SNPBMI,
        logHRSNP, logHRBMI, logHRM,
        logHRM1M2, logHRSNP_M1M2,
        coef(mod_allcohort)[2],
        coef(mod_subsamplecohort)[2],
        coef(mod_NCC_ctrl)[2],
        coef(mod_NCC)[2],
        coef(mod_casectrlWeightedKM_matching)[2],
        coef(mod_casectrlWeightedGAM_nomatching)[2],
        coef(mod_casectrlWeightedGAM_matching)[2],
        coef(mod_casectrlWeightedGAM_matchingInter)[2],
        coef(modCox_allcohort)[1:2],
        coef(mod_clogit)[1:2],
        coef(modCox_casectrlWeightedKM_matching)[1:2],
        coef(modCox_casectrlWeightedGAM_nomatching)[1:2],
        coef(modCox_casectrlWeightedGAM_matching)[1:2],
        coef(modCox_casectrlWeightedGAM_matchingInter)[1:2],
        KM_FullCohort,
        KM_casectrlWeightedKM_matching,
        KM_casectrlWeightedGAM_nomatching,
        KM_casectrlWeightedGAM_matching,
        KM_casectrlWeightedGAM_matchingInter,
        KM_cond_FullCohort,
        KM_cond_casectrlWeightedKM_matching,
        KM_cond_casectrlWeightedGAM_nomatching,
        KM_cond_casectrlWeightedGAM_matching,
        KM_cond_casectrlWeightedGAM_matchingInter,
        sum(cohort_withweights$Ind),
        sum(cohort_withweights$WeightsGAM_noMatching > 100),
        sum(cohort_withweights$WeightsGAM_noMatching > 1000),
        sum(cohort_withweights$WeightsGAM_Matching > 100),
        sum(cohort_withweights$WeightsGAM_Matching > 1000))
    )
  }
  allnames <- c(
    "nsim", "pselCases", "MAF", "R2SNPBMI",
    "logHRSNP", "logHRBMI",
    "logHRM", "logHRM1M2", "logHRSNP_M1M2",
    "beta_allcohort", "beta_allcohort_subsample",
    "beta_NCC_ctrl",
    "beta_NCC",
    "beta_NCC_WeightedKM_nomatching",
    "beta_NCC_WeightedKM_matching",
    "beta_NCC_WeightedGAM_nomatching",
    "beta_NCC_WeightedGAM_matching",
    "beta_NCC_WeightedGAM_matchingInter",
    "alphaSNP_allcohort", "alphaBMI_allcohort",
    "alphaSNP_NCC", "alphaBMI_NCC",
    "alphaSNP_NCC_WeightedKM_nomatching",
    "alphaBMI_NCC_WeightedKM_nomatching",
    "alphaSNP_NCC_WeightedKM_matching",
    "alphaBMI_NCC_WeightedKM_matching",
    "alphaSNP_NCC_WeightedGAM_nomatching",
    "alphaBMI_NCC_WeightedGAM_nomatching",
    "alphaSNP_NCC_WeightedGAM_matching",
    "alphaBMI_NCC_WeightedGAM_matching",
    "alphaSNP_NCC_WeightedGAM_matchingInter",
    "alphaBMI_NCC_WeightedGAM_matchingInter",
    "KM_allcohort",
    "KM_NCC_WeightedKM_nomatching",
    "KM_NCC_WeightedKM_matching",
    "KM_NCC_WeightedGAM_nomatching",
    "KM_NCC_WeightedGAM_matching",
    "KM_NCC_WeightedGAM_matchingInter",
    "KM_cond_allcohort",
    "KM_cond_NCC_WeightedKM_nomatching",
    "KM_cond_NCC_WeightedKM_matching",
    "KM_cond_NCC_WeightedGAM_nomatching",
    "KM_cond_NCC_WeightedGAM_matching",
    "KM_cond_NCC_WeightedGAM_matchingInter",
    "nbCases",
    "nbGAMWeightsnoMatchingabove100",
    "nbGAMWeightsnoMatchingabove1000",
    "nbGAMWeightsMatchingabove100",
    "nbGAMWeightsMatchingabove1000"
  )

  colnames(RECAP_RES_ASSO) <- allnames[-which(grepl("KM_nomatching", allnames))]
  AllRes <- list(RECAP_RES_ASSO = RECAP_RES_ASSO, Counts = Counts)
  saveRDS(AllRes, file = paste0("Res/Res_", suffix, ".rds"))
}

################################################################################
# Run the function
################################################################################
res <- mclapply(1:nrow(ALLPRMTR), fun_all_para, mc.cores = numWorkers)