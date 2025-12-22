
compute_GAM_weights <- function(cohort,
                                id.var       = "Id",
                                status.var   = "statforweight_03",
                                case.ind     = "Ind",
                                samplestatus = "indic_included",
                                survtime.var = "s(LengthVV)",
                                match.vars   = "s(Age_Blood)",
                                TD           = TRUE,
                                cases.vars   = "s(LengthVV) + CenterFac",
                                typical.ncc  = FALSE) {
  # Select required vars                                
  IdVar              <- cohort %>% select(all_of(id.var)) %>% pull
  Selvar             <- cohort %>% select(all_of(status.var)) %>% pull
  IndCasevar         <- cohort %>% select(all_of(case.ind)) %>% pull
  forgamweights      <- cohort %>% mutate(Id = IdVar,
                                          VarSel = Selvar,
                                          IndCase = IndCasevar)

  # Here we select the variables that enter the GAM model
  if (TD) {
    VarsModCtls <- paste0(c(survtime.var, match.vars), collapse = " + ")
  } else {
    VarsModCtls <- paste(match.vars, collapse = " + ")
  }
  # Here we append them
  formula_GamWeigths <- paste0(samplestatus, " ~ ", VarsModCtls)
  formula_Cases      <- paste0("I(VarSel >= 2) ~ ", cases.vars)
  modGAM_ctls        <- gam(as.formula(formula_GamWeigths),
                            family = "binomial",
                            data = subset(forgamweights, VarSel != 2))
  # Computation for Untypical NCC
  if (!typical.ncc) {
      # Probs for cases
    modGAM_cases   <- gam(as.formula(formula_Cases), method="REML", family="binomial", data=subset(forgamweights, IndCase==1))
    pGAM_cases     <- as.numeric(predict(modGAM_cases, newdata = forgamweights, type='response'))
      # Probs for controls
    pGAM_ctls      <- as.numeric(predict(modGAM_ctls, newdata = forgamweights, type='response'))
    forgamweights  <- forgamweights %>% mutate(Pred_Interest = case_when(IndCase==0 ~ pGAM_ctls, TRUE ~ pGAM_cases + (1-pGAM_cases)*pGAM_ctls))
  # Computation for Typical NCC
  } else {
      # Probs for controls
    pGAM_ctls          <- as.numeric(predict(modGAM_ctls, newdata = forgamweights[forgamweights$VarSel !=2,], type='response'))
    forgamweights_ctrls <- forgamweights %>% filter(VarSel!=2) %>% mutate(Pred_Interest = pGAM_ctls) 
      # Probs for cases
    forgamweights_cases <- forgamweights %>% filter(VarSel==2) %>% mutate(Pred_Interest = 1) 
    forgamweights_all   <- bind_rows(forgamweights_ctrls, forgamweights_cases)
    forgamweights       <- forgamweights %>% left_join(forgamweights_all %>% select(Id, Pred_Interest))
  }
  res <- forgamweights %>%
          mutate(WeightsGAM = case_when(VarSel == 0 ~ 0, TRUE ~ 1 / Pred_Interest), 
                 gamforall  = 1 / Pred_Interest) %>%
          select(Id, WeightsGAM, gamforall)
  return(res)
}


compute_weights_KM_Typical <- function(cohort, matching = T, closest = F, tolMatching=ifelse(closest==T, 0.1, 0.1*sqrt(12)))#, interMonSurv=F)
{
  # if(!interMonSurv)
  # {
    if(!matching) 
    {
      KMWeights <- 1/KMprob(cohort$AgeOut, cohort$statforweight_02, m=1)
    }else
    {
      if (closest==F)
      {
        KMWeights <- 1/KMprob(cohort$AgeOut, cohort$statforweight_02, m=1, left.time = 0, match.var = cohort$Matching, match.int = c(-tolMatching, tolMatching))
      }else
      {
        KMWeights <- ifelse(cohort$statforweight_02==0, 0, 1)
      }
    }
  # }else{
  #   if(!matching) 
  #   {
  #     KMWeights <- 1/KMprob(cohort$AgeOut, cohort$statforweight_02, m=1)
  #   }else
  #   {
  #     KMWeights <- 1/KMprob(cohort$AgeOut, cohort$statforweight_02, m=1, left.time = 0, match.var = cbind(cohort$M1, cohort$M2), match.int = rep(0, 4))
  #   }
  #   
  # }
  res          <- cohort %>% mutate(WeightsKM = case_when(statforweight_02==0~0, 
                                                             statforweight_02==2~1,
                                                             TRUE~KMWeights)) %>%
                             select(Id, WeightsKM)
  return(res)
}

compute_weights_KM_Untypical <- function(cohort, NCC, pselCases,
                                         ListRiskSets, Nties = rep(1, length(ListRiskSets)), 
                                         VarIdCohort = "Id", VarIndCaseCohort = "Ind", VarStatforWeight = "statforweight", 
                                         VarIdNCC = "Id", VarIndCaseNCC = "Ind")
{
  
  if (is.null(ListRiskSets)){stop("KM-type weights will not be computed as ListRiskSets is required for their computation when pselCases < 1")} #| is.null(ListRiskSets_noMatching)
  
  IdVar              <- cohort %>% select(all_of(VarIdCohort)) %>% pull
  Selvar             <- cohort %>% select(all_of(VarStatforWeight)) %>% pull
  IndCasevar         <- cohort %>% select(all_of(VarIndCaseCohort)) %>% pull
  dataforWeightsVV   <- cohort %>% mutate(Id = IdVar,
                                          VarSel = Selvar,
                                          IndCase = IndCasevar)
  
  DatawithWeights   <- NCC %>% mutate(Ind = get(VarIndCaseNCC), Id= get(VarIdNCC)) %>% arrange(desc(Ind)) %>% distinct(Id, .keep_all = T)
  
  ## VV KM weights
  if (!(is.null(ListRiskSets))) #
  {
    fun_Probs_SelasControls <- function(integ, ListRiskSets)
    {
      prob <- 1
      for (ii in 1:length(ListRiskSets))
      {
        RiskSet <- ListRiskSets[[ii]]
        nTies   <- Nties[[ii]]
        if (integ %in% RiskSet)
        {
          prob <- prob * (1 - (nTies/(length(RiskSet) - 1)))  # Horwitz-Thompson weights (see Letter to the editor regarding the paper, Biom. J.)
          # \delta_j*V_{1,j} = 1 for all RiskSets, by construction:
          # the risk sets are only for the cases selected in the NCC)
        }
      }
      return(1- prob)
    }
    
    Ind_inNCC        <- DatawithWeights %>%  dplyr::select(Id) %>% pull %>% unique  # %>% arrange(Id, -Ind) %>% distinct(Id, .keep_all = T)
    Prob_SelAsCtrls  <- pbsapply(Ind_inNCC, function(integ){fun_Probs_SelasControls(integ, ListRiskSets)})
    
    
    dataforKMWeights   <- DatawithWeights %>% mutate(WeightsKM = case_when(
                                                        Ind==0 ~ 1/Prob_SelAsCtrls,
                                                        Ind==1 ~ 1/(pselCases + (1-pselCases)*Prob_SelAsCtrls)))
    
    res                <- left_join(dataforWeightsVV,
                                    dataforKMWeights %>% dplyr::select(Id, WeightsKM), by="Id") %>% #
                          mutate(WeightsKM = case_when(VarSel==0 ~ 0,
                                                       TRUE ~ WeightsKM)) %>%
                          select(Id, WeightsKM)
    return(res)
  }
  
}

# DEPRECATED: UNUSED IN THE INTERACTION CASE. ONLY CENSORING or No T, No D, which is not in the paper)
compute_weights_GAM_Typical <- function(cohort, matching = TRUE, TD = TRUE, quantileC = NA,
                                closest = FALSE,
                                tolMatching = ifelse(closest == TRUE, 0.1,
                                                     0.1 * sqrt(12)),
                                interMonSurv = FALSE,
                                interforMwheninteronSurv = FALSE) {
  if (!interMonSurv) {
    if (TD == TRUE)
    {
      if (!matching) {
        wgam <- 1 / GAMprob(cohort$AgeOut, cohort$statforweight_02)
      } else {
        wgam <- 1 / GAMprob(cohort$AgeOut, cohort$statforweight_02, left.time = 0, match.var = cohort$Matching, match.int = c(-tolMatching, tolMatching))  
      }
    }
    if (TD=="no_TD") {
      if(!matching) {
        modGLM  <- glm(indic_included ~ 1, family=binomial, data=cohort)
        wgam    <- 1/as.vector(predict(modGLM, newdata=cohort, type='response'))
      } else {
        modGAM <- gam(indic_included ~ s(Matching), method="REML", family=binomial, data=cohort)
        wgam   <- 1/as.vector(predict(modGAM, newdata=cohort, type='response')) 
      }
    }
    if (TD=="no_T")
    {
      if(!matching) 
      {
        modGLM  <- glm(indic_included ~ 1, family=binomial, data=subset(cohort, statforweight_02!=2))
        wgam    <- 1/as.vector(predict(modGLM, newdata=cohort, type='response'))
      }else{
        modGAM <- gam(indic_included ~ s(Matching), method="REML", family=binomial, data=subset(cohort, statforweight_02!=2))
        wgam   <- 1/as.vector(predict(modGAM, newdata=cohort, type='response')) 
      }
    }
    if(TD == "no_low_C") {
      if(!matching) {
        modGLM1  <- glm(indic_included ~ 1, family=binomial, data=subset(cohort, statforweight_02!=2))
        wgam1    <- 1/as.vector(predict(modGLM1, newdata=cohort, type='response'))
        modGAM2  <- gam(indic_included ~ s(AgeOut), method="REML", family=binomial, data=subset(cohort, statforweight_02!=2))
        wgam2    <- 1/as.vector(predict(modGLM2, newdata=cohort, type='response'))
        wgam     <- wgam2
        wgam[cohort$AgeOut < quantile(cohort$AgeOut[cohort$statforweight_02==2], probs = quantileC)] <- 
                                      wgam1[cohort$AgeOut < quantile(cohort$AgeOut[cohort$statforweight_02==2], probs = quantileC)]
      } else {
        modGAM1 <- gam(indic_included ~ s(Matching), method="REML", family=binomial, data=subset(cohort, statforweight_02!=2))
        wgam1   <- 1/as.vector(predict(modGAM1, newdata=cohort, type='response')) 
        modGAM2 <- gam(indic_included ~ s(AgeOut) + s(Matching), method="REML", family=binomial, data=subset(cohort, statforweight_02!=2))
        wgam2   <- 1/as.vector(predict(modGAM2, newdata=cohort, type='response')) 
        wgam    <- wgam2
        wgam[cohort$AgeOut < quantile(cohort$AgeOut[cohort$statforweight_02==2], probs = quantileC)]<- wgam1[cohort$AgeOut < quantile(cohort$AgeOut[cohort$statforweight_02==2], probs = quantileC)]
      }
    }
      
    if (TD!="no_TD")
    {
      res          <- cohort %>% mutate(WeightsGAM = case_when(statforweight_02==0~0, 
                                                             statforweight_02==2~1,
                                                             TRUE~wgam)) %>%
                               select(Id, WeightsGAM)
    }else{
      res          <- cohort %>% mutate(WeightsGAM = case_when(statforweight_02==0~0, 
                                                               TRUE~wgam)) %>%
        select(Id, WeightsGAM)
    }
  }else{
    if(!matching) 
    {
      wgam <- 1/GAMprob(cohort$AgeOut, cohort$statforweight_02)
      res          <- cohort %>% mutate(WeightsGAM = case_when(statforweight_02==0~0, 
                                                               statforweight_02==2~1,
                                                               TRUE~wgam)) %>%
                                  select(Id, WeightsGAM)
    } else {
      cohort <- cohort %>% mutate(M2fac = as_factor(round(M2, 2)))
      if (!interforMwheninteronSurv) {
        res <- compute_weights_GAM_Untypical(cohort, matching = TRUE, TD = TRUE,
                                      VarId            = "Id",
                                      VarStatforWeight = "statforweight_03",
                                      VarIndCase       = "Ind",
                                      VarIndicModCtrls = "indic_included",
                                      VarLenghtFup     = "s(AgeOut)",
                                      VarAdjCtls       = "s(M1) + M2fac",
                                      VarAdjCases      = "1",
                                      AllCases         = TRUE)
      } else {
        res <- compute_weights_GAM_Untypical(cohort, matching = TRUE, TD = TRUE,
                                             VarId            = "Id",
                                             VarStatforWeight = "statforweight_03",
                                             VarIndCase       = "Ind",
                                             VarIndicModCtrls = "indic_included",
                                             VarLenghtFup     = "s(AgeOut)",
                                             VarAdjCtls       = "s(M1) + M2fac + s(M1, by=M2fac)",
                                             VarAdjCases      = "1",
                                             AllCases         = TRUE)
      }
    }
  }
  return(res)
}