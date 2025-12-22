simul_cohort <- function(n = 1e5, R2SNPBMI, logHRSNP, logHRBMI,
                         logHRM, logHRM1M2, logHRSNP_M1M2,
                         MeanBMI, SDBMI, MeanSNP, SDSNP, 
                         MeanM1, SDM1, MeanM2, SDM2,
                         MeanM1M2, SDM1M2,
                         MeanSNP_M1M2, SDSNP_M1M2,
                         Asso_M_SNP, MAF,
                         agemincens, agemaxcens) {
  ###### Exposures
  if (prod(is.na(MeanBMI),is.na(SDBMI),is.na(MeanSNP),is.na(SDSNP), is.na(MeanM1), is.na(SDM1), is.na(MeanM2), is.na(SDM2))==1) {
  # Generate Means and SDs for first simulation
    # Matching variables
    ntest         <- 1e6
    M1test        <- scale(runif(ntest))
    M2test        <- runif(ntest)
    M2test        <- ifelse(M2test > 0.5, 1, 0)
    M1M2test      <- M1test * M2test
    # Means and SDs
    MeanM1        <- mean(M1test)
    SDM1          <- sd(M1test)
    MeanM2        <- mean(M2test)
    SDM2          <- sd(M2test)
    MeanM1M2      <- mean(M1M2test)
    SDM1M2        <- sd(M1M2test)
    MAFtest       <- MAF + Asso_M_SNP * (M1test <= -1 & M2test == 0)
    # Simulate the SNPS for the first simulation
    forSNPstest   <- t(sapply(1:ntest, function(i){rmultinom(1, 1, c(1-2*MAFtest[i]+MAFtest[i]^2, 2*MAFtest[i]*(1-MAFtest[i]), MAFtest[i]^2))}))
    SNPstest      <- apply(forSNPstest, 1, function(v){which(v ==1) - 1})
    MeanSNP       <- mean(SNPstest)
    SDSNP         <- sd(SNPstest)
    SNPstest      <- scale(SNPstest)
    # Interaction SNPs with Ms
    SNP_M1M2test  <- SNPstest*M1M2test
    MeanSNP_M1M2  <- mean(SNP_M1M2test)
    SDSNP_M1M2    <- sd(SNP_M1M2test)
    # Interaction SNPs with BMI
    snpbmitest    <- SNPstest * 0.1
    alpha2test    <- var(snpbmitest)
    sigma2bmitest <- alpha2test * (1 - R2SNPBMI) / R2SNPBMI
    if (R2SNPBMI != 0) {
      BMIstest    <- rnorm(ntest, mean = snpbmitest, sd = sqrt(sigma2bmitest))
      MeanBMI     <- mean(BMIstest)
      SDBMI       <- sd(BMIstest)
    }
  }
  # Second Simulation onwards
    # Simulate Matching Variables
    M1              <- scale(runif(n))
    M2              <- runif(n)
    M2              <- ifelse(M2 > 0.5, 1, 0)
    M1M2            <- M1 * M2
    M1M2            <- (M1M2 - MeanM1M2) / SDM1M2
    # Simulate SNPs
    MAFoo           <- MAF + Asso_M_SNP * (M1 <= -1 & M2 == 0)
    forSNPs         <- t(sapply(1:n, function(i){rmultinom(1, 1, c(1-2*MAFoo[i]+MAFoo[i]^2, 2*MAFoo[i]*(1-MAFoo[i]), MAFoo[i]^2))}))
    SNPs            <- (apply(forSNPs, 1, function(v){which(v==1) - 1}) - MeanSNP) / SDSNP
    # Interaction SNPs with Ms
    SNP_M1M2        <- SNPs * M1M2
    SNP_M1M2        <- (SNP_M1M2 - MeanSNP_M1M2) / SDSNP_M1M2
    # Interaction SNPs with BMI
    snpbmi          <- SNPs * 0.1
    alpha2          <- var(snpbmi)
    sigma2bmi       <- alpha2 * (1 - R2SNPBMI) / R2SNPBMI
    if (R2SNPBMI == 0) {BMIs      <- rnorm(n)} else {
      BMIs <- (rnorm(n, mean = snpbmi, sd = sqrt(sigma2bmi)) - MeanBMI) / SDBMI
    }
    BMIs        <- pmin(pmax(BMIs, -3), 3)

  # Compute Survival (https://www.econstor.eu/bitstream/10419/31002/1/483869465.PDF)
  COVS              <- cbind(BMIs, SNPs, M1, M2, M1M2, SNP_M1M2)
  linear_pred       <- COVS %*% matrix(c(logHRBMI, logHRSNP, logHRM, logHRM, logHRM1M2, logHRSNP_M1M2), nrow=6) 
  linear_pred       <- linear_pred - mean(linear_pred)
  lam0              <- - log(1 - 0.01 / median(exp(linear_pred))) / 40 #0.000005
  survival_times    <- rexp(n, rate = 1) / (lam0 * exp(linear_pred))
  censoring_times   <- runif(n, min = agemincens, max = agemaxcens)
  observed_times    <- pmin(survival_times, censoring_times)
  delta             <- 1 * (survival_times <= censoring_times)
  if (mean(delta) > 0.07 | mean(delta) < 0.06) {
    n_iter <- 0
    if (mean(delta) > 0.07) {
      while (mean(delta) > 0.07 & n_iter < 20) {
        n_iter          <- n_iter + 1
        lam0            <- lam0 / 1.1
        survival_times  <- rexp(n, rate = 1) / (lam0 * exp(linear_pred))
        censoring_times <- runif(n, min = agemincens, max = agemaxcens)
        observed_times  <- pmin(survival_times, censoring_times)
        delta           <- 1 * (survival_times <= censoring_times)
      }
    }
    if (mean(delta) < 0.06) {
      while (mean(delta) < 0.06 & n_iter < 20) {
        n_iter         <- n_iter + 1
        lam0           <- lam0 * 1.1
        survival_times <- rexp(n, rate = 1) / (lam0 * exp(linear_pred))
        censoring_times <- runif(n, min = agemincens, max = agemaxcens)
        observed_times  <- pmin(survival_times, censoring_times)
        delta           <- 1 * (survival_times <= censoring_times)
      }
    }
  }
  # Store cohort data
  CohortData  <- tibble(
    Id = 1:n, SNP = SNPs, BMI = BMIs,
    M1 = as.numeric(M1), M2 = as.numeric(M2),
    AgeOut = as.numeric(observed_times),
    Ind = as.numeric(delta)
  )
  Counts      <- CohortData %>% count(Ind, M1, M2)

  return(list(CohortData       = CohortData,
              MeanBMI          = MeanBMI,
              SDBMI            = SDBMI,
              MeanSNP          = MeanSNP,
              SDSNP            = SDSNP,
              MeanM1           = MeanM1,
              SDM1             = SDM1,
              MeanM2           = MeanM2,
              MeanM1M2         = MeanM1M2,
              SDM1M2           = SDM1M2,
              SDM2             = SDM2,
              MeanSNP_M1M2     = MeanSNP_M1M2,
              SDSNP_M1M2       = SDSNP_M1M2,
              Counts           = Counts))
}


sample_case_control <- function(cohort, pselCases, returnRiskSet = FALSE) {
  # Select cases to be chosen
  NTotCases    <- sum(cohort$Ind)
  indCases     <- which(cohort$Ind == 1)
  NSelCases    <- pselCases * NTotCases
  # Extract cases based on Typical or Untypical NCC
  if (pselCases == 1) {
    indCasessel  <- indCases
  } else {
    indCasessel  <- sample(indCases, round(NSelCases), replace = FALSE)
  }
  # Obtain NCC
  NCC          <- cohort %>% mutate(Pair = "aa") %>% slice(0)
  if (pselCases < 1 | returnRiskSet) {
    ListRiskSets_noMatching <- ListRiskSets_Matching <- vector("list", NSelCases)
  } else {
    ListRiskSets_noMatching <- ListRiskSets_Matching <- NULL
  }
  for (ii in 1:NSelCases) {
    icase              <- indCasessel[ii] 
    idcase             <- cohort$Id[icase]
    M1_case            <- cohort$M1[icase]
    M2_case            <- cohort$M2[icase]
    ageout_case        <- cohort$AgeOut[icase]
    # Select Risk Set by Caliper Matching (M2 binary so actually exact)
    RiskSet_noMatching <- cohort %>% filter(AgeOut >= ageout_case)
    RiskSet_Matching   <- cohort %>% filter(AgeOut >= ageout_case,
                                            abs(M1 - M1_case) < sqrt(3) / 5,
                                            abs(M2 - M2_case) < sqrt(3) / 5)
    # Make sure at least one control is in the Risk Set
    # a) For Typical NCC
    if (nrow(RiskSet_Matching) == 1) {
      eps1 <- eps2 <- sqrt(3) / 5
      keps <- 0
      while (nrow(RiskSet_Matching) == 1) {
        whichonetorelax <- sample(c(1, 2), 1)
        if (whichonetorelax == 1) {
          eps1 <- eps1 * 2
        } else {
          eps2 <- eps2 * 2
        }
        # Select with broader tolerance
        RiskSet_Matching   <- cohort %>% filter(AgeOut >= ageout_case,
                                                abs(M1 - M1_case) < eps1,
                                                abs(M2 - M2_case) < eps2)
        # Check again whether at least 1 control
        if (nrow(RiskSet_Matching) == 1) {
          if(whichonetorelax == 2) {
            eps1 <- eps1 * 2
          } else {
            eps2 <- eps2 * 2
          }
          RiskSet_Matching   <- cohort %>% filter(AgeOut >= ageout_case,
                                                  abs(M1 - M1_case) < eps1,
                                                  abs(M2 - M2_case) < eps2)
        }
        keps <- keps + 1
      }
    }
    # b) For Untypical NCC
    if (pselCases < 1 | returnRiskSet) {
      ListRiskSets_Matching[[ii]] <- RiskSet_Matching$Id
    }

    # Extract the sample
    tempsample  <- RiskSet_Matching %>% filter(Id != idcase) %>% select(Id) %>% pull

    if (length(tempsample) > 1) {
      icontrol  <- sample(tempsample, 1)
    } else {
      icontrol  <- tempsample
    }
    # Extract the NCC
    NCC <- bind_rows(NCC, cohort %>% slice(c(icase, icontrol)) %>%
                       mutate(IndCaseCtrlNCC = c(1, 0),
                              Pair = paste0("PairCase_", ii)))
  }
  return(list(NCC = NCC,
              ListRiskSets_noMatching = ListRiskSets_noMatching,
              ListRiskSets_Matching = ListRiskSets_Matching))
}

sample_case_control_NN2 <- function(cohort, pselCases, returnRiskSet = FALSE) {
  # Total cases and indices
  NTotCases <- sum(cohort$Ind)
  indCases  <- which(cohort$Ind == 1)
  NSelCases <- if (pselCases == 1) NTotCases else round(pselCases * NTotCases)

  # Sample cases if needed
  if (pselCases == 1) {
    indCasessel <- indCases
  } else {
    indCasessel <- sample(indCases, NSelCases, replace = FALSE)
  }

  NCC <- cohort %>% mutate(Pair = "aa") %>% slice(0)
  if (pselCases < 1 | returnRiskSet) {
    ListRiskSets_Matching <- vector("list", NSelCases)
  } else {
    ListRiskSets_Matching <- NULL
  }

  for (ii in seq_len(NSelCases)) {
    icase       <- indCasessel[ii]
    idcase      <- cohort$Id[icase]
    M1_case     <- cohort$M1[icase]
    M2_case     <- cohort$M2[icase]
    ageout_case <- cohort$AgeOut[icase]

    # Risk set: subjects at risk at case's time and with exact matching on M2
    RiskSet <- cohort %>% filter(AgeOut >= ageout_case, M2 == M2_case)

    if (pselCases < 1 | returnRiskSet) {
      ListRiskSets_Matching[[ii]] <- RiskSet$Id
    }

    # Remove the case itself
    potential_controls <- RiskSet %>% filter(Id != idcase)

    if (nrow(potential_controls) == 0) {
      warning(paste("No matching control found for case ID", idcase))
      next
    }

    # Nearest neighbor on M1: find control with minimal absolute difference in M1
    potential_controls <- potential_controls %>%
      mutate(diff_M1 = abs(M1 - M1_case)) %>%
      arrange(diff_M1)

    # Select closest match (first one)
    icontrol_id <- potential_controls$Id[1]

    # Bind case and control as pair
    NCC <- bind_rows(
      NCC,
      cohort %>% filter(Id %in% c(idcase, icontrol_id)) %>%
        mutate(IndCaseCtrlNCC = ifelse(Id == idcase, 1, 0),
               Pair = paste0("PairCase_", ii))
    )
  }

  return(list(NCC = NCC,
              ListRiskSets_Matching = ListRiskSets_Matching))
}

sample_case_control_NN <- function(cohort, pselCases) {
  NTotCases   <- sum(cohort$Ind)
  indCases    <- which(cohort$Ind == 1)
  NSelCases   <- round(pselCases * NTotCases)
  indCasessel <- sample(indCases, NSelCases, replace = FALSE)

  # Initialize NCC output
  NCC <- cohort %>% mutate(Pair = "aa") %>% slice(0)
  ListRiskSets <- if (pselCases < 1) vector("list", NSelCases) else NULL

  for (ii in seq_len(NSelCases)) {
    icase       <- indCasessel[ii]
    idcase      <- cohort$Id[icase]

    m1_case     <- cohort$M1[icase]
    m2_case     <- cohort$M2[icase]
    ageout_case <- cohort$AgeOut[icase]
    RiskSet     <- which(cohort$AgeOut >= ageout_case)

    if (pselCases < 1) {
      ListRiskSets[[ii]] <- RiskSet
    }

    # Find eligible controls in the risk set
    potential_controls <- cohort %>%
      slice(RiskSet) %>%
      filter(Id != idcase, M2 == m2_case) %>%
      mutate(diff_match = abs(M1 - m1_case)) %>%
      arrange(diff_match)

    if (nrow(potential_controls) == 0) {
      warning(paste("No matching control found for case ID", idcase))
      next
    }

    # Choose closest match (first one)
    icontrol_id <- potential_controls %>%
      slice(1) %>%
      pull(Id)

    # Bind matched pair to output
    NCC <- bind_rows(
      NCC,
      cohort %>%
        filter(Id %in% c(idcase, icontrol_id)) %>%
        mutate(
          IndCaseCtrlNCC = ifelse(Id == idcase, 1, 0),
          Pair = paste0("PairCase_", ii)
        )
    )
  }

    # Post-matching summary of problematic pairs
  pairwise_variance_summary <- NCC %>%
    group_by(Pair) %>%
    summarise(
      var_M1 = var(M1),
      var_SNP = var(SNP),
      var_M1M2 = var(M1 * M2),
      var_M1M2SNP = var(M1 * M2 * SNP),
      all_M2_zero = all(M2 == 0)
    ) %>%
    mutate(
      zero_M1 = var_M1 == 0,
      zero_SNP = var_SNP == 0,
      zero_M1M2 = var_M1M2 == 0,
      zero_M1M2SNP = var_M1M2SNP == 0
    )

  n_pairs <- nrow(pairwise_variance_summary)
  n_M2_zero <- sum(pairwise_variance_summary$all_M2_zero)
  pct_M2_zero <- round(100 * n_M2_zero / n_pairs, 2)

  # Print summary counts
  message("Summary of zero variance across pairs:")
  message("Pairs with zero var(M1): ", sum(pairwise_variance_summary$zero_M1))
  message("Pairs with zero var(SNP): ", sum(pairwise_variance_summary$zero_SNP))
  message("Pairs with zero var(M1*M2): ", sum(pairwise_variance_summary$zero_M1M2))
  message("Pairs with zero var(M1*M2*SNP): ", sum(pairwise_variance_summary$zero_M1M2SNP))
  message(sprintf("Pairs with all M2 = 0: %d (%.2f%%)", n_M2_zero, pct_M2_zero))

  n_missing_controls <- sum(!(paste0("PairCase_", seq_len(NSelCases)) %in% NCC$Pair[NCC$IndCaseCtrlNCC == 1]))
  message("Cases without matched control: ", n_missing_controls)

  return(list(NCC = NCC, ListRiskSets_Matching = ListRiskSets))
}


prepare_for_weights <- function(cohort, NCC) {
  ncc         <- NCC %>% arrange(Id, - IndCaseCtrlNCC) %>% distinct(Id, .keep_all = T)
  ind_ncccase <- ncc %>% filter(IndCaseCtrlNCC==1) %>% select(Id) %>% pull
  ind_nccctrl <- ncc %>% filter(IndCaseCtrlNCC==0) %>% select(Id) %>% pull
  cohort      <- cohort %>%
    mutate(
           statforweight_03 = case_when(
      (Id %in% ind_ncccase & Id %in% ind_nccctrl) ~ 3,
       Id %in% ind_ncccase ~ 2,
       Id %in% ind_nccctrl ~ 1,
                                        TRUE ~ 0),
      statforweight_02 = case_when(Id %in% ind_ncccase ~ 2,
                                   Id %in% ind_nccctrl ~ 1,
                                   TRUE ~ 0),
      indic_included = 1 * (statforweight_02 > 0), M1M2 = M1 * M2)
}


compute_weights <- function(cohort, NCC, pselCases = 1,
                            ListRiskSets_Matching = ListRiskSets_Matching) {

  TypicalNCC_value <- (pselCases == 1)  # TRUE if pselCases==1, otherwise FALSE

  cohort <- cohort %>% mutate(M2fac = as.factor(M2))
  gam_matched   <- compute_GAM_weights(cohort,
                                       id.var       = "Id",
                                       status.var   = "statforweight_03",
                                       case.ind     = "Ind",
                                       samplestatus = "indic_included",
                                       survtime.var = "s(AgeOut)",
                                       match.vars   = "s(M1) + M2fac",
                                       TD           = TRUE,
                                       cases.vars   = NULL,
                                       typical.ncc  = TypicalNCC_value)

  gam_unmatched <- compute_GAM_weights(cohort,
                                       id.var       = "Id",
                                       status.var   = "statforweight_03",
                                       case.ind     = "Ind",
                                       samplestatus = "indic_included",
                                       survtime.var = "s(AgeOut)",
                                       match.vars   = "1",
                                       TD           = TRUE,
                                       cases.vars   = NULL,
                                       typical.ncc  = TypicalNCC_value)

  gam_matchedinter <- compute_GAM_weights(cohort,
                                          id.var       = "Id",
                                          status.var   = "statforweight_03",
                                          case.ind     = "Ind",
                                          samplestatus = "indic_included",
                                          survtime.var = "s(AgeOut)",
                                          match.vars   = "s(M1) + M2fac + s(M1, by=M2fac)",
                                          TD           = TRUE,
                                          cases.vars   = NULL,
                                          typical.ncc  = TypicalNCC_value)

  if (!is.null(ListRiskSets_Matching)) {
    km_matched <- compute_weights_KM_Untypical(cohort, NCC, pselCases,
                                               ListRiskSets_Matching,
                                               VarIdCohort = "Id",
                                               VarIndCaseCohort = "Ind",
                                               VarStatforWeight = "statforweight_02",
                                               VarIdNCC = "Id",
                                               VarIndCaseNCC = "Ind")
  } else { # NN case has no returned risk set if Typical NCC
    km_matched <- compute_weights_KM_Typical(cohort, matching = TRUE,
                                             closest = TRUE)
                                             #,  tolMatching = TRUE)
  }

  cohort <- cohort %>% 
    left_join(gam_unmatched %>%
                rename(WeightsGAM_noMatching = WeightsGAM), by = "Id") %>%
    left_join(gam_matched %>%
                rename(WeightsGAM_Matching = WeightsGAM), by = "Id") %>%
    left_join(gam_matchedinter %>%
                rename(WeightsGAM_MatchingInter = WeightsGAM), by = "Id")  %>%
    left_join(km_matched %>%
                rename(WeightsKM_Matching = WeightsKM), by = "Id") %>%
    mutate(
      WeightsGAM_noMatching_thresh = pmin(WeightsGAM_noMatching, 100),
      WeightsGAM_Matching_thresh   = pmin(WeightsGAM_Matching, 100),
      WeightsCohort                = 1
    )
  return(cohort)
}