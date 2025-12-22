simul_cohort <- function(n = 1e5, R2SNPBMI, logHRSNP, logHRBMI,
                         logHRM, logHRM1M2, logHRSNP_M1M2,
                         MeanBMI, SDBMI, MeanSNP, SDSNP, 
                         MeanM1, SDM1, MeanM2, SDM2,
                         MeanM1M2, SDM1M2,
                         MeanSNP_M1M2, SDSNP_M1M2,
                         MAF, agemincens, agemaxcens) {
  ###### Exposures
  if (prod(is.na(MeanBMI),is.na(SDBMI),is.na(MeanSNP),is.na(SDSNP), is.na(MeanM1), is.na(SDM1), is.na(MeanM2), is.na(SDM2))==1) {
  # Generate Means and SDs for first simulation
    # Matching variables
    ntest         <- 1e6
    M1test        <- runif(ntest)
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
    MAFtest       <- rep(MAF, ntest)
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
    sigma2bmitest <- alpha2test*(1 - R2SNPBMI) / R2SNPBMI
    if (R2SNPBMI != 0) {
      BMIstest    <- rnorm(ntest, mean = snpbmitest, sd = sqrt(sigma2bmitest))
      MeanBMI     <- mean(BMIstest)
      SDBMI       <- sd(BMIstest)
    }
  }
  # Second Simulation onwards
    # Simulate Matching Variables
    M1              <- runif(n)
    M2              <- runif(n)
    M2              <- ifelse(M2 > 0.5, 1, 0)
    M1M2            <- M1 * M2
    M1M2            <- (M1M2 - MeanM1M2) / SDM1M2
    # Simulate SNPs
    MAFoo           <- rep(MAF, n)
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
      BMIs <- (rnorm(n, mean= snpbmi, sd = sqrt(sigma2bmi)) - MeanBMI) / SDBMI
    }
    BMIs        <- pmin(pmax(BMIs, -3), 3)

  # Compute Survival (https://www.econstor.eu/bitstream/10419/31002/1/483869465.PDF)
  COVS              <- cbind(BMIs, SNPs, M1, M2, M1M2, SNP_M1M2)
  linear_pred       <- COVS %*% matrix(c(logHRBMI, logHRSNP, logHRM, logHRM, logHRM1M2, logHRSNP_M1M2), nrow=6) 
  linear_pred       <- linear_pred - mean(linear_pred)
  lam0              <- - log(1 - 0.01/median(exp(linear_pred))) / 40 #0.000005
  survival_times    <- rexp(n, rate = 1) / (lam0 * exp(linear_pred))
  censoring_times   <- runif(n, min = agemincens, max = agemaxcens)
  observed_times    <- pmin(survival_times, censoring_times)
  delta             <- 1 * (survival_times <= censoring_times)
  if (mean(delta) > 0.07 | mean(delta) < 0.02) {
    n_iter <- 0
    if (mean(delta) > 0.07) {
      while (mean(delta) > 0.07 & n_iter < 10) {
        n_iter          <- n_iter + 1
        lam0            <- lam0 / 5
        survival_times  <- rexp(n, rate = 1) / (lam0 * exp(linear_pred))
        censoring_times <- runif(n, min = agemincens, max = agemaxcens)
        observed_times  <- pmin(survival_times, censoring_times)
        delta           <- 1 * (survival_times <= censoring_times)
      }
    }
    if (mean(delta) < 0.02) {
      while (mean(delta) < 0.02 & n_iter < 10) {
        n_iter         <- n_iter + 1
        lam0           <- lam0 * 5
        survival_times <- rexp(n, rate = 1) / (lam0 * exp(linear_pred))
        censoring_times <- runif(n, min = agemincens, max = agemaxcens)
        observed_times  <- pmin(survival_times, censoring_times)
        delta           <- 1 * (survival_times <= censoring_times)
      }
    }
  }
  # Store cohort data
  CohortData  <- tibble(Id = 1:n, SNP = SNPs, BMI = BMIs,
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
    # Select Risk Set by Caliper Matching
    RiskSet_noMatching <- cohort %>% filter(AgeOut >= ageout_case)
    RiskSet_Matching   <- cohort %>% filter(AgeOut >= ageout_case,
                                            abs(M1 - M1_case) < 0.1 * sqrt(12),
                                            abs(M2 - M2_case) < 0.1 * sqrt(12))
    # Make sure at least one control is in the Risk Set
    # a) For Typical NCC
    if (nrow(RiskSet_Matching) == 1) {
      eps1 <- eps2 <- 0.1 * sqrt(12)
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


####### Univ M 

simul_cohort_univM <- function(n=1e5,R2SNPBMI, logHRSNP, logHRBMI, logHRMatching, 
                               MeanBMI, SDBMI, MeanSNP, SDSNP, MAF, agemincens, agemaxcens,
                               asso_M_SNP = T)
{
  
  ###### Exposures
  if (prod(is.na(MeanBMI),is.na(SDBMI),is.na(MeanSNP),is.na(SDSNP))==1)
  {
    ntest         <- 1e6
    Matchingstest <- scale(runif(ntest))
    
    if(asso_M_SNP)
    {
      MAFtest     <- MAF + 0.2*(Matchingstest <= -1 )
      forSNPstest <- t(sapply(1:ntest, function(i){rmultinom(1, 1, c(1-2*MAFtest[i]+MAFtest[i]^2, 2*MAFtest[i]*(1-MAFtest[i]), MAFtest[i]^2))}))#t(rmultinom(ntest, 1, c(1-2*MAF+MAF^2, 2*MAF*(1-MAF), MAF^2)))
    }else{
      forSNPstest   <- t(rmultinom(ntest, 1, c(1-2*MAF+MAF^2, 2*MAF*(1-MAF), MAF^2)))  
    }
    
    SNPstest      <- apply(forSNPstest, 1, function(v){which(v==1)-1})
    MeanSNP       <- mean(SNPstest)
    SDSNP         <- sd(SNPstest)
    SNPstest      <- scale(SNPstest)
    snpbmitest    <- SNPstest * 0.1
    alpha2test    <- var(snpbmitest) #sum(snpbmitest^2)
    sigma2bmitest <- alpha2test*(1-R2SNPBMI)/R2SNPBMI #alpha2test*(1-R2SNPBMI)/(n*R2SNPBMI)
    if (R2SNPBMI !=0) 
    {
      BMIstest       <- rnorm(ntest, mean= snpbmitest, sd = sqrt(sigma2bmitest))
      MeanBMI        <- mean(BMIstest)
      SDBMI          <- sd(BMIstest)
    }
  }
  Matchings  <- scale(runif(n))
  if(asso_M_SNP)
  {
    MAFoo    <- MAF + 0.2*(Matchings <= -1 )
    forSNPs  <- t(sapply(1:n, function(i){rmultinom(1, 1, c(1-2*MAFoo[i]+MAFoo[i]^2, 2*MAFoo[i]*(1-MAFoo[i]), MAFoo[i]^2))}))#t(rmultinom(ntest, 1, c(1-2*MAF+MAF^2, 2*MAF*(1-MAF), MAF^2)))
  }else{
    forSNPs   <- t(rmultinom(n, 1, c(1-2*MAF+MAF^2, 2*MAF*(1-MAF), MAF^2)))  
  }
  SNPs       <- (apply(forSNPs, 1, function(v){which(v==1)-1}) - MeanSNP)/SDSNP
  snpbmi     <- SNPs * 0.1
  alpha2     <- var(snpbmi) #sum(snpbmi^2)
  sigma2bmi  <- alpha2*(1-R2SNPBMI)/R2SNPBMI # alpha2*(1-R2SNPBMI)/(n*R2SNPBMI)
  if (R2SNPBMI ==0) {
    BMIs      <- rnorm(n)
    # MeanBMI   <- 0
    # SDBMI     <- 1
  }else{
    BMIs      <- (rnorm(n, mean= snpbmi, sd = sqrt(sigma2bmi)) - MeanBMI)/SDBMI
  }
  
  BMIs        <- pmin(pmax(BMIs,-3), 3)
  COVS        <- cbind(BMIs, SNPs, Matchings)
  
  linear_pred <- COVS %*% matrix(c(logHRBMI, logHRSNP, logHRMatching), nrow=3)
  
  survival_times  <- 70 * (-log(runif(n))*exp(-linear_pred))^{1/3} #; hist(survival_times)
  # https://www.econstor.eu/bitstream/10419/31002/1/483869465.PDF
  #Generating Survival Times to Simulate Cox Proportional Hazards Models Ralf Bender1, Thomas Augustin2, Maria Blettner
  
  censoring_times <- runif(n, min=agemincens, max=agemaxcens) #pmin(rweibull(n, scale = (agemincens + agemaxcens)/2, shape = 2), agemaxcens) 
  observed_times  <- pmin(survival_times, censoring_times) 
  delta           <- 1*(survival_times<=censoring_times)
  #coxph(Surv(observed_times, delta)~ BMIs + SNPs + Matchings)
  
  CohortData    <- tibble(Id = 1:n, SNP = SNPs, BMI = BMIs, Matching = as.numeric(Matchings),
                          AgeOut = as.numeric(observed_times), Ind = as.numeric(delta))
  
  return(list(CohortData       = CohortData, 
              MeanBMI          = MeanBMI,
              SDBMI            = SDBMI,
              MeanSNP          = MeanSNP,
              SDSNP            = SDSNP))
}



simul_cohort_interMatching <- function(n=1e5,R2BMI, logHRSNP, logHRBMI, logHRMatching, 
                                       MeanBMI, SDBMI, MeanSNP, SDSNP, MAF, agemincens, agemaxcens)
{
  
  
  
  ###### Exposures
  if (prod(is.na(MeanBMI),is.na(SDBMI),is.na(MeanSNP),is.na(SDSNP))==1)
  {
    ntest         <- 1e6
    Matchingstest <- scale(runif(ntest))
    
    MAFtest     <- MAF + 0.2*(Matchingstest <= -1 )
    forSNPstest <- t(sapply(1:ntest, function(i){rmultinom(1, 1, c(1-2*MAFtest[i]+MAFtest[i]^2, 2*MAFtest[i]*(1-MAFtest[i]), MAFtest[i]^2))}))#t(rmultinom(ntest, 1, c(1-2*MAF+MAF^2, 2*MAF*(1-MAF), MAF^2)))
    SNPstest      <- apply(forSNPstest, 1, function(v){which(v==1)-1})
    MeanSNP       <- mean(SNPstest)
    SDSNP         <- sd(SNPstest)
    SNPstest      <- scale(SNPstest)
    
    Matchingstest <- scale(runif(ntest))
    SNPstest      <- scale(SNPstest)
    snpbmitest    <- SNPstest * 0.1 * Matchingstest
    alpha2test    <- var(snpbmitest) #sum(snpbmitest^2)
    sigma2bmitest <- alpha2test*(1-R2BMI)/R2BMI #alpha2test*(1-R2SNPBMI)/(n*R2SNPBMI)
    if (R2BMI !=0) 
    {
      BMIstest       <- rnorm(ntest, mean= snpbmitest, sd = sqrt(sigma2bmitest))
      MeanBMI        <- mean(BMIstest)
      SDBMI          <- sd(BMIstest)
    }
  }
  
  Matchings  <- scale(runif(n))
  MAFoo    <- MAF + 0.2*(Matchings <= -1 )
  forSNPs  <- t(sapply(1:n, function(i){rmultinom(1, 1, c(1-2*MAFoo[i]+MAFoo[i]^2, 2*MAFoo[i]*(1-MAFoo[i]), MAFoo[i]^2))}))#t(rmultinom(ntest, 1, c(1-2*MAF+MAF^2, 2*MAF*(1-MAF), MAF^2)))
  SNPs       <- (apply(forSNPs, 1, function(v){which(v==1)-1}) - MeanSNP)/SDSNP
  snpbmi     <- SNPs * 0.1 * Matchings
  alpha2     <- var(snpbmi) #sum(snpbmi^2)
  sigma2bmi  <- alpha2*(1-R2BMI)/R2BMI # alpha2*(1-R2SNPBMI)/(n*R2SNPBMI)
  if (R2BMI ==0) {
    BMIs      <- rnorm(n)
    # MeanBMI   <- 0
    # SDBMI     <- 1
  }else{
    BMIs      <- (rnorm(n, mean= snpbmi, sd = sqrt(sigma2bmi)) - MeanBMI)/SDBMI
  }
  
  BMIs        <- pmin(pmax(BMIs,-3), 3)
  COVS        <- cbind(BMIs, SNPs, Matchings)
  
  linear_pred <- COVS %*% matrix(c(logHRBMI, logHRSNP, logHRMatching), nrow=3)
  
  survival_times  <- 70 * (-log(runif(n))*exp(-linear_pred))^{1/3} #; hist(survival_times)
  # https://www.econstor.eu/bitstream/10419/31002/1/483869465.PDF
  #Generating Survival Times to Simulate Cox Proportional Hazards Models Ralf Bender1, Thomas Augustin2, Maria Blettner
  
  censoring_times <- runif(n, min=agemincens, max=agemaxcens) #pmin(rweibull(n, scale = (agemincens + agemaxcens)/2, shape = 2), agemaxcens) 
  observed_times  <- pmin(survival_times, censoring_times) 
  delta           <- 1*(survival_times<=censoring_times)
  #coxph(Surv(observed_times, delta)~ BMIs + SNPs + Matchings)
  
  CohortData    <- tibble(Id = 1:n, SNP = SNPs, BMI = BMIs, Matching = as.numeric(Matchings),
                          AgeOut = as.numeric(observed_times), Ind = as.numeric(delta))
  
  return(list(CohortData       = CohortData, 
              MeanBMI          = MeanBMI,
              SDBMI            = SDBMI,
              MeanSNP          = MeanSNP,
              SDSNP            = SDSNP))
}


sample_case_control_univM_sample_NNMatching <- function(cohort, pselCases)
{
  NTotCases    <- sum(cohort$Ind)
  indCases     <- which(cohort$Ind==1)
  NSelCases    <- pselCases*NTotCases
  indCasessel  <- sample(indCases, round(NSelCases), replace=F)
  NCC          <- cohort %>% mutate(Pair = "aa") %>% slice(0)
  if (pselCases < 1){ListRiskSets_noMatching <- vector("list", NSelCases); ListRiskSets_Matching <- vector("list", NSelCases)}else{ListRiskSets_noMatching <- NULL; ListRiskSets_Matching <- NULL}
  for (ii in 1:NSelCases)
  {
    icase          <- indCasessel[ii] 
    idcase         <- cohort$Id[icase]
    matching_case  <- cohort$Matching[icase]
    RiskSet        <- which(cohort$AgeOut >= cohort$AgeOut[icase])
    icontrol       <- cohort %>% slice(RiskSet) %>% 
                                  filter(Id != idcase) %>% 
                                  mutate(diff_match = abs(Matching - matching_case)) %>% 
                                  arrange(diff_match) %>% slice(1) %>% select(Id) %>% pull
   
    if (pselCases < 1)
    {
      ListRiskSets_noMatching[[ii]] <- RiskSet
      ListRiskSets_Matching[[ii]]   <- icontrol
    }
    
    NCC            <- bind_rows(NCC, cohort %>% slice(c(icase, icontrol)) %>% mutate(IndCaseCtrlNCC = c(1, 0), Pair = paste0("PairCase_", ii)))
  }
  #cat("not here anymore, right \n")
  return(list(NCC=NCC, ListRiskSets_noMatching=ListRiskSets_noMatching, ListRiskSets_Matching=ListRiskSets_Matching))
}

sample_case_control_univM_sample_CaliperMatching <- function(cohort, pselCases, tolMatching=0.1*sqrt(12), returnRiskSet=F)
{
  NTotCases    <- sum(cohort$Ind)
  indCases     <- which(cohort$Ind==1)
  NSelCases    <- pselCases*NTotCases
  indCasessel  <- sample(indCases, round(NSelCases), replace=F)
  NCC          <- cohort %>% mutate(Pair = "aa") %>% slice(0)
  if (pselCases < 1 | returnRiskSet){ListRiskSets_noMatching <- ListRiskSets_Matching <- vector("list", NSelCases)}else{ListRiskSets_noMatching <- ListRiskSets_Matching <- NULL}
  for (ii in 1:NSelCases)
  {
    icase              <- indCasessel[ii] 
    idcase             <- cohort$Id[icase]
    matching_case      <- cohort$Matching[icase]
    ageout_case        <- cohort$AgeOut[icase]
    RiskSet_noMatching <- cohort %>% filter(AgeOut >= ageout_case)
    RiskSet_Matching   <- cohort %>% filter(AgeOut >= ageout_case & abs(Matching - matching_case)<= tolMatching)
    if (pselCases < 1 | returnRiskSet)
    {
      ListRiskSets_noMatching[[ii]] <- RiskSet_noMatching$Id
      ListRiskSets_Matching[[ii]]   <- RiskSet_Matching$Id
    }
    #icontrol       <- sample(RiskSet_Matching$Id, 1)
    tempsample     <- RiskSet_Matching %>% filter(Id != idcase) %>% select(Id) %>% pull
    if (length(tempsample)>1) 
    {
      icontrol       <- sample(tempsample, 1)
    }else{icontrol <- tempsample}
    NCC            <- bind_rows(NCC, cohort %>% slice(c(icase, icontrol)) %>% mutate(IndCaseCtrlNCC = c(1, 0), Pair = paste0("PairCase_", ii)))
  }
  #cat("not here anymore, right \n")
  return(list(NCC=NCC, ListRiskSets_noMatching=ListRiskSets_noMatching, ListRiskSets_Matching=ListRiskSets_Matching)) # 
}

prepare_for_weights <- function(cohort, NCC, univar=FALSE) {
  ncc         <- NCC %>% arrange(Id, - IndCaseCtrlNCC) %>% distinct(Id, .keep_all = T)
  ind_ncccase <- ncc %>% filter(IndCaseCtrlNCC==1) %>% select(Id) %>% pull
  ind_nccctrl <- ncc %>% filter(IndCaseCtrlNCC==0) %>% select(Id) %>% pull
  cohort      <- cohort %>%
    mutate(statforweight_03 = case_when(
                                        (Id %in% ind_ncccase & Id %in% ind_nccctrl) ~ 3,
                                        Id %in% ind_ncccase ~ 2,
                                        Id %in% ind_nccctrl ~ 1,
                                        TRUE ~ 0),
           statforweight_02 = case_when(Id %in% ind_ncccase ~ 2,
                                        Id %in% ind_nccctrl ~ 1,
                                        TRUE ~ 0),
           indic_included = 1 * (statforweight_02 > 0))
  
if(!univar){ cohort      <- cohort %>% mutate(M1M2 = M1 * M2)}
  return(cohort)
}


compute_weights <- function(cohort, NCC, pselCases = 1,
                            ListRiskSets_Matching = NULL,
                            univariate = FALSE, nnmatch=FALSE, 
                            withTDontop = FALSE) { ### if withTDontop is TRUE, the function also returns GAM weights wo T,D (used in the simulation study of the Supp Mat. when focusing on the asso. between exposures)
  TypicalNCC_value <- (pselCases == 1)  # TRUE if pselCases==1, otherwise FALSE
  if (TypicalNCC_value){cases_var <- NULL}else{cases_var <-"1"}
  # If univariate, only use M1 for matching
  if (!univariate) {
    cohort <- cohort %>% mutate(M2fac = as.factor(M2))
    match_vars <- "s(M1) + M2fac"
  } else {
    match_vars <- "s(Matching)"
  }
  gam_matched   <- compute_GAM_weights(cohort,
                                       id.var       = "Id",
                                       status.var   = "statforweight_03",
                                       case.ind     = "Ind",
                                       samplestatus = "indic_included",
                                       survtime.var = "s(AgeOut)",
                                       match.vars   = match_vars,
                                       TD           = TRUE,
                                       cases.vars   = cases_var,
                                       typical.ncc  = TypicalNCC_value)
  
  if(withTDontop){
    gam_noTD <- compute_GAM_weights(cohort,
                                    id.var       = "Id",
                                    status.var   = "statforweight_03",
                                    case.ind     = "Ind",
                                    samplestatus = "indic_included",
                                    survtime.var = "s(AgeOut)",
                                    match.vars   = match_vars,
                                    TD           = FALSE,
                                    cases.vars   = cases_var,
                                    typical.ncc  = TypicalNCC_value)
  }else{gam_noTD <- NULL}
  
  
  gam_unmatched <- compute_GAM_weights(cohort,
                                       id.var       = "Id",
                                       status.var   = "statforweight_03",
                                       case.ind     = "Ind",
                                       samplestatus = "indic_included",
                                       survtime.var = "s(AgeOut)",
                                       match.vars   = "1",
                                       TD           = TRUE,
                                       cases.vars   = cases_var,
                                       typical.ncc  = TypicalNCC_value)
  # Only compute gam_matchedinter if univariate is FALSE
  if (!univariate) {
    gam_matchedinter <- compute_GAM_weights(cohort,
                                            id.var       = "Id",
                                            status.var   = "statforweight_03",
                                            case.ind     = "Ind",
                                            samplestatus = "indic_included",
                                            survtime.var = "s(AgeOut)",
                                            match.vars   = "s(M1) + M2fac + s(M1, by=M2fac)",
                                            TD           = TRUE,
                                            cases.vars   = cases_var,
                                            typical.ncc  = TypicalNCC_value)
  } else {
    gam_matchedinter <- NULL  # No need for interaction model in the univariate case
  }
  if(TypicalNCC_value)
  {
    km_matched   <- compute_weights_KM_Typical(cohort, matching=T, closest = nnmatch)
    km_unmatched <- compute_weights_KM_Typical(cohort, matching=F, closest = nnmatch)
  }else{
  km_matched       <- compute_weights_KM_Untypical(cohort, NCC, pselCases,
                                                   ListRiskSets_Matching,
                                                   VarIdCohort = "Id",
                                                   VarIndCaseCohort = "Ind",
                                                   VarStatforWeight = "statforweight_02",
                                                   VarIdNCC = "Id",
                                                   VarIndCaseNCC = "Ind")
  km_unmatched    <- NULL
  }
  # Join results, conditionally joining gam_matchedinter if it's not NULL
  cohort <- cohort %>% 
    left_join(gam_unmatched %>%
                rename(WeightsGAM_noMatching = WeightsGAM), by = "Id") %>%
    left_join(gam_matched %>%
                rename(WeightsGAM_Matching = WeightsGAM), by = "Id") %>%
    # Only join gam_matchedinter if it's not NULL
    left_join(if (!is.null(gam_matchedinter)) {
      gam_matchedinter %>% rename(WeightsGAM_MatchingInter = WeightsGAM)
    } else {
      tibble(Id = integer(0))
    }, by = "Id") %>%
    # Only join gam_noTD if it's not NULL
    left_join(if (!is.null(gam_noTD)) {
      gam_noTD %>% rename(WeightsGAM_noTD = WeightsGAM)
    } else {
      tibble(Id = integer(0))
    }, by = "Id") %>%
    left_join(km_matched %>%
                rename(WeightsKM_Matching = WeightsKM), by = "Id") %>%
    left_join(if (!is.null(km_unmatched)) {
      km_unmatched %>% rename(WeightsKM_noMatching = WeightsKM)
    } else {
      tibble(Id = integer(0))
    }, by = "Id") %>%
    mutate(WeightsGAM_noMatching_thresh = pmin(WeightsGAM_noMatching, 100),
           WeightsGAM_Matching_thresh   = pmin(WeightsGAM_Matching, 100),
           WeightsCohort                = 1)
  return(cohort)
}






