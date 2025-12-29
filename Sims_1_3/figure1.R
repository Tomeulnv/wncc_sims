library(tidyverse)
library(survival)
library(survminer)

ROOT <- Sys.getenv("WNCC_ROOT")
stopifnot(nzchar(ROOT))

setwd(file.path(ROOT, "wncc_sims", "Sims_1_3"))

source("/RFunctions/Functions_Analysis_EPIC_Biocrates.R")

StudyInterest  <- "MB_CORU_01"
comp_res  <- readRDS(
  paste0("Res/forweights_all_", StudyInterest, "_June2025.rds")
)
forweights_all <- comp_res$forweights_all
forweights_all <- forweights_all %>%
  mutate(
    Cncr_Mal_Interest = haven::zap_labels(Cncr_Mal_Interest),
    Weights_Selectable = case_when(
      Cncr_Mal_Interest == 1 | Idepic %in% comp_res$all_eligible_KM ~ 1,
      TRUE ~ 0
    )
  )

forweights_all <- forweights_all %>%
  mutate(
    WeightsNCC_all_Total = case_when(Sel_all > 0 ~ 1, TRUE ~ 0),
    WeightsNCC_all_Ctls = case_when(Sel_all == 1 ~ 1, TRUE ~ 0),
    WeightsNCC_int_Total = case_when(Sel_Interest > 0 ~ 1, TRUE ~ 0),
    WeightsNCC_int_Ctls = case_when(Sel_Interest == 1 ~ 1, TRUE ~ 0),
    Weights_Selectable = case_when(
      Cncr_Mal_Interest == 1 | Idepic %in% comp_res$all_eligible_KM ~ 1, TRUE ~ 0
    )
)

fun_KM_estimate(
  forweights_all, StudyInterest,
  WeightsToKeep = c("WeightsGAM_NoTek", "WeightsGAM_WithTek", "WeightsKM", "WeightsKM_woTek", "Weights_FullCohort", "Weights_Selectable"), 
  levelsWeights = c("Weights_FullCohort", "Weights_Selectable", "WeightsGAM_WithTek", "WeightsGAM_NoTek", "WeightsKM", "WeightsKM_woTek"), 
  legendWeights = c("Originating Cohort", "Eligible", paste0("GAM-Weights", c("", " wo Tech.")), paste0("KM-Weights", c("", " wo Tech.")))
)
