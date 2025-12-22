#### Extraction of the originating cohort from EPIC, that is the list of EPIC participants that fulfill the inclusion criteria
fun_eligible_subpop_interest <- function(StudyInterest, CncrInterest, TypeTumoInterest,
                                         Country_toexclude, Cond_FilterInterest,
                                         pmaxVV)
{
  id_study_interest        <- Biocrates_Normalized %>% filter(Study==StudyInterest)
  id_cases_Study_interest  <- Biocrates_Normalized %>% filter(Study==StudyInterest, Cncr_Caco ==1)
  
  #### exclusion criteria Endometrial cancer
  EligiblePopulationInterestBiocrat <- EPICeligible %>% 
                                          mutate(Age_maxBldRecr = pmaxVV(Age_Blood, Age_Recr)) %>%
                                          filter(!Country %in% Country_toexclude,  
                                                               Cancer_Prev==0,  
                                                               Baseline_Data_Avail==1, 
                                                               Fup_Vital_Status==1, 
                                                               Bld_Sample_Avail==1, 
                                                               Fup_Cncr_Blood == 1, 
                                                               Age_Exit >= Age_maxBldRecr)
                                                         
  filter_Spec <- paste0("EligiblePopulationInterestBiocrat <- EligiblePopulationInterestBiocrat %>% filter(", 
                        Cond_FilterInterest, ")")
  
  eval(parse(text = filter_Spec))
  
  ### this is the study-specific part... to be "functionalized"...
  
  if (mean(id_study_interest$Idepic %in% EligiblePopulationInterestBiocrat$Idepic)!=1){stop("EPICBiocratesks like inclusion/exclusion criteria are not correct")}
  
  
  EligiblePopulationInterestBiocrat <- EligiblePopulationInterestBiocrat %>% 
                                              left_join(EPICendpoints %>% select_if(grepl("Idepic|Cncr_|Typ_Tumo", names(.)))) %>%
                                              mutate(Cncr_Mal_Interest = get(CncrInterest), 
                                                     Typ_Tumo_Interest = get(TypeTumoInterest)) %>%
                                              mutate(Cncr_Mal_Interest = replace_na(Cncr_Mal_Interest, 0))
  
  ################################
  ####
  #### Checking a few things
  ####
  ################################
  
  EligiblePopulationInterestBiocrat %>% filter(Idepic %in% id_cases_Study_interest$Idepic) %>% select(Cncr_Mal_Interest) %>% pull %>% summary %>% print()
  EligiblePopulationInterestBiocrat %>% count(Cncr_Mal_Interest) %>% print()
  
  ### Problem
  # EligiblePopulationInterestBiocrat %>% filter(Age_Exit > Age_Csr_Cancer2010 + 0.001 & Cncr_Mal_Interest==1, Idepic %in% id_cases_Study_interest$Idepic)
  ### => 0 individuals d'Asturias (+68 si on n'ajoute pas le "+0.001"..mais c'est sans doute dû au calcul approché de Age_Csr_Cancer2010) ... 
  ################################
  ################################
  # Cas particulier du 41____41156061 (par exemple) qui a une tumeur avec une date de diag a missing, on exclut car on ne sait pas la date du 1er cancer => à virer donc 
  ################################
  ################################
  ################################
  EligiblePopulationInterestBiocrat %>% count(Typ_Tumo_Interest) %>% print()
  EligiblePopulationInterestBiocrat %>% filter(Idepic %in%id_cases_Study_interest$Idepic) %>% count(Typ_Tumo_Interest) %>% print()
  EligiblePopulationInterestBiocrat %>% filter(Cncr_Mal_Interest ==1, Idepic %in% id_cases_Study_interest$Idepic) %>% count(Typ_Tumo_Interest) %>% print()
  EligiblePopulationInterestBiocrat %>% filter(Cncr_Mal_Interest ==1, !(Idepic %in% id_cases_Study_interest$Idepic)) %>% count(Typ_Tumo_Interest) %>% print()
  EligiblePopulationInterestBiocrat %>% filter(Cncr_Mal_Interest ==1, Idepic %in% Biocrates_Normalized$Idepic) %>% count(Typ_Tumo_Interest) %>% print()
  EligiblePopulationInterestBiocrat %>% filter(is.na(Cncr_Mal_Interest), Idepic %in% Biocrates_Normalized$Idepic) %>% count(Typ_Tumo_Interest) %>% print()
  EligiblePopulationInterestBiocrat %>% filter(Cncr_Mal_Interest ==1, !(Idepic %in% id_cases_Study_interest$Idepic), Typ_Tumo_Interest!= "09-Excluded")  %>% print()
  
  
  weird <- Biocrates_Normalized %>% filter(Study==StudyInterest, !(Idepic%in%EligiblePopulationInterestBiocrat$Idepic)) %>% select(Idepic) %>% pull
  EPICeligible %>% filter(Idepic %in% weird) %>% print()
  
  TypTumo_Included <- EligiblePopulationInterestBiocrat %>% filter(Cncr_Mal_Interest ==1, Idepic %in% id_cases_Study_interest$Idepic) %>%
                          select(all_of(TypeTumoInterest)) %>% pull %>% unique
  
  #############################################
  ####
  #### Create the Cncr_Interest  
  #### and the Selection indicators
  ####
  #############################################
  EligiblePopulationInterestBiocrat_VV <- EligiblePopulationInterestBiocrat %>% 
    mutate(LengthVV = (Age_Exit - Age_maxBldRecr)*365.25,
           Cncr_Interest_VV = case_when(
             is.na(Cncr_Mal_Interest)~0, 
             ### this is debatable: NA means that this individual reported 2 tumors the same day.. 
             ### one corresponding to the cancer of interest, another one to another cancer.. not clear whether he/she should be considered as a case or not... 
             Cncr_Mal_Interest==0~0, 
             (Cncr_Mal_Interest ==1 & !(Typ_Tumo_Interest %in%  TypTumo_Included) )~0, 
             TRUE~1)) %>%
    left_join(Biocrates_Normalized %>% select(Idepic, Cncr_Caco) %>% arrange(desc(Cncr_Caco)) %>% distinct(Idepic, .keep_all = T) %>% mutate(InBiocrates = 1)) %>%
    left_join(Biocrates_Normalized %>% filter(Study == StudyInterest) %>% select(Idepic, Cncr_CacoInterest = Cncr_Caco) %>% arrange(desc(Cncr_CacoInterest)) %>% distinct(Idepic, .keep_all = T) %>% mutate(InBiocratesInterest = 1)) %>%
    mutate(InBiocrates         = replace_na(InBiocrates, 0), 
           InBiocratesInterest = replace_na(InBiocratesInterest, 0)) %>% 
    mutate(
      SexFac                   = as_factor(Sex),
      IndFac                   = as_factor(Cncr_Interest_VV),
      TbldMissing              = 1*(is.na(T_Bld_Coll)), 
      forFourierTbld           = case_when(is.na(T_Bld_Coll)~0, T~minute(T_Bld_Coll) + 60 *hour(T_Bld_Coll)), 
      FourierTblood1           = sin(2*pi*forFourierTbld/(24*60)),
      FourierTblood2           = cos(2*pi*forFourierTbld/(24*60)),
      FourierTblood3           = sin(4*pi*forFourierTbld/(24*60)),
      FourierTblood4           = cos(4*pi*forFourierTbld/(24*60)),
      DbldMissing              = 1*(is.na(D_Bld_Coll)),
      dayofyear_bdcol          = case_when(is.na(D_Bld_Coll)~0, T~yday(D_Bld_Coll)), 
      FourierDblood1           = sin(2*pi*dayofyear_bdcol/365),
      FourierDblood2           = cos(2*pi*dayofyear_bdcol/365),
      FourierDblood3           = sin(4*pi*dayofyear_bdcol/365),
      FourierDblood4           = cos(4*pi*dayofyear_bdcol/365),
      forTBldVV                = T_Bld_Coll,
      Menop_BldNum             = as.numeric(as.character(Menop_Bld)), 
      CenterFac                = as_factor(Center), 
      CountryFac               = as_factor(Country), 
      Fasting_Fac              = as_factor(replace_na(Fasting_C, 99)),
      Menop_BldFac             = as_factor(case_when(
                                      Sex==1~0,
                                      TRUE~Menop_Bld)),
      Phrt_BldFac              = as_factor(replace_na(case_when(
                                      Sex==1~0, 
                                      TRUE~Phrt_Bld), 99)), 
      Phase_MnsCyVV            = as_factor(replace_na(ifelse(Menop_Bld %in% c(1, 3), 0, Phase_Mnscycle), 99))
    )
  
  EligiblePopulationInterestBiocrat_VV %>% count(Cncr_Interest_VV, Cncr_Mal_Interest)%>% print()
  # EligiblePopulationInterestBiocrat_VV %>% count(Cncr_Interest_VV, InBiocrates)%>% print()
  # EligiblePopulationInterestBiocrat_VV %>% count(Cncr_Interest_VV, InBiocrates, Cncr_Caco)%>% print()
  # EligiblePopulationInterestBiocrat_VV %>% count(Cncr_Interest_VV, InBiocrates, Cncr_CacoInterest)%>% print()
  EligiblePopulationInterestBiocrat_VV %>% count(Cncr_Interest_VV, InBiocratesInterest, Cncr_CacoInterest)%>% print()
  if (mean(id_study_interest$Idepic %in% EligiblePopulationInterestBiocrat_VV$Idepic)!=1) {stop("Something got wrong here...")}
  return(EligiblePopulationInterestBiocrat_VV)
}



fun_computation_Weights_EPIC <- function(IndexMets, StudyInterest, EligiblePopulationInterestBiocrat_VV, 
                                         StrictMatchingInterest, FlexMatchingInterest, listToleranceInterest, 
                                         Adj_GamWeights_Interest, ## should be a list, whose elements are vectors of factors accounted for in the computation of GAM-weights
                                         Adj_GamWeights_All=NULL, ## one single vector can be given in the current version; could be extended eventually; but for now, I'm not really considering the "All" thing anymore... 
                                         FlexMatchVarTechnicalInterest=NULL)
{

  AllMets          <- colnames(Biocrates_Normalized)[IndexMets]
  temp_Interest    <- Biocrates_Normalized %>% filter(Study == StudyInterest) #%>% arrange(desc(Cncr_Caco)) %>% distinct(Idepic, .keep_all = T) #%>%
  #left_join(selections %>% select(Idepic, which(grepl("Cncr_Mal", colnames(selections)))))
  casesBio_Interest <- temp_Interest %>% filter(Cncr_Caco == 1) %>% select(Idepic) %>% pull
  ctrlsBio_Interest <- temp_Interest %>% filter(Cncr_Caco == 0) %>% select(Idepic) %>% pull
  Bio_all_nocaseInt <- Biocrates_Normalized %>% filter(Study != StudyInterest | Cncr_Caco == 0)%>% select(Idepic) %>% pull
  
  EligiblePopulationInterestBiocrat_VV     <- EligiblePopulationInterestBiocrat_VV %>% 
                                                      mutate(Sel_Interest         = case_when ((Idepic %in% casesBio_Interest & Idepic %in% ctrlsBio_Interest) ~ 3,
                                                                                       Idepic %in% casesBio_Interest ~2,
                                                                                       Idepic %in% ctrlsBio_Interest ~1, 
                                                                                       TRUE ~ 0), 
                                                         Sel_all                  = case_when ((Idepic %in% casesBio_Interest & Idepic %in% Bio_all_nocaseInt) ~3, 
                                                                                       Idepic %in% casesBio_Interest ~2,
                                                                                       Idepic %in% Bio_all_nocaseInt ~1, 
                                                                                       TRUE ~ 0),
                                                         indic_inluded_all        = 1*(Sel_all>0), 
                                                         indic_inluded_Interest   = 1*(Sel_Interest>0)) 
  
    ####### Construction of the gam weights
    
    forgamweights <-    EligiblePopulationInterestBiocrat_VV %>% left_join(Biocrates_Normalized %>% distinct(Idepic, .keep_all = T) %>% 
                                                                 select(Idepic, all_of(AllMets)))  %>% #Bmi_C is now in the eligible population already ;) 
                                                                 # replace_na(list(T_Bld_Coll = 9999999, Fasting_C = 999, forTBldVV=0)) %>% ### missing values for T_Bld will be considered as 0:00 (or most T_Bld=0 will be considered as missing in the modeling part...)
                                                                 # mutate(D_Blood_VV      = as_factor(year(D_Bld_Coll)), 
                                                                 #        Fasting_VV    = as_factor(Fasting_C), 
                                                                 #        T_Blood_VV    = as_factor(hour(forTBldVV))) %>% 
                                                                 mutate_if(grepl("Cncr_Mal_", names(.)), ~replace_na(., 0)) 
 
  forweights_all <- forgamweights
  for (le in 1:length(Adj_GamWeights_Interest))
  {
    adj_gamweights       <- Adj_GamWeights_Interest[[le]]
    suff_adj_gam         <- names(Adj_GamWeights_Interest)[le]
    WeightsGAM_suff      <- paste0("WeightsGAM_", suff_adj_gam)
    ProbInclGAM_suff     <- paste0("ProbInclGAM_", suff_adj_gam)
    gam_matched_interest <- compute_weights_GAM_Untypical(cohort=forgamweights, matching = T, TD=T, 
                                                        VarId            = "Idepic",
                                                        VarStatforWeight = "Sel_Interest",
                                                        VarIndCase       = "IndFac",
                                                        VarIndicModCtrls = "indic_inluded_Interest", 
                                                        VarLenghtFup     = "s(LengthVV)", 
                                                        VarAdjCtls       = c("s(Age_Blood)", adj_gamweights), 
                                                        VarAdjCases      = "s(LengthVV)") 
    
    forweights_all <- forweights_all %>% left_join(gam_matched_interest %>% rename(Idepic = Id, !!WeightsGAM_suff:=WeightsGAM, !!ProbInclGAM_suff:=ProbInclGAM)) 
  }
  
  if(!is.null(Adj_GamWeights_All))
  {
    gam_matched_all <- compute_weights_GAM_Untypical(cohort = forgamweights, matching = T, TD=T, 
                                                   VarId            = "Idepic",
                                                   VarStatforWeight = "Sel_all",
                                                   VarIndCase       = "IndFac",
                                                   VarIndicModCtrls = "indic_inluded_all", 
                                                   VarLenghtFup     = "s(LengthVV)", 
                                                   VarAdjCtls       = c("s(Age_Blood)", Adj_GamWeights_All), #paste0(Adj_GamWeights_Interest, collapse=" + ")
                                                   VarAdjCases      = "s(LengthVV)") 
    
    forweights_all <- forweights_all %>% left_join(gam_matched_all %>% rename(Idepic = Id, WeightsGAM_all=WeightsGAM, ProbInclGAM_all=ProbInclGAM)) 
    
  }
  
  # gam_matched_interest_noTD <- compute_weights_GAM_Untypical(cohort=forgamweights, matching = T, TD=F, 
  #                                                       VarId            = "Idepic",
  #                                                       VarStatforWeight = "Sel_Interest",
  #                                                       VarIndCase       = "IndFac",
  #                                                       VarIndicModCtrls = "indic_inluded_Interest", 
  #                                                       VarLenghtFup     = "s(LengthVV)", 
  #                                                       VarAdjCtls       = c("s(Age_Blood)", Adj_GamWeights_Interest), #paste0(Adj_GamWeights_Interest, collapse=" + ")
  #                                                       VarAdjCases      = "1") 
  # 
  # 
  # 
  # gam_matched_all_noTD <- compute_weights_GAM_Untypical(cohort = forgamweights, matching = T, TD=F, 
  #                                                  VarId            = "Idepic",
  #                                                  VarStatforWeight = "Sel_all",
  #                                                  VarIndCase       = "IndFac",
  #                                                  VarIndicModCtrls = "indic_inluded_all", 
  #                                                  VarLenghtFup     = "s(LengthVV)", 
  #                                                  VarAdjCtls       = c("s(Age_Blood)", Adj_GamWeights_All), #paste0(Adj_GamWeights_Interest, collapse=" + ")
  #                                                  VarAdjCases      = "1") 
  
  #  left_join(gam_matched_all_noTD %>% rename(Idepic = Id, WeightsGAM_all_noTD=WeightsGAM, ProbInclGAM_all_noTD=ProbInclGAM)) %>%
  #  left_join(gam_matched_interest_noTD %>% rename(Idepic = Id, WeightsGAM_int_noTD=WeightsGAM, ProbInclGAM_int_noTD=ProbInclGAM)) %>%
 
  ####### Computation of KM weights
  selCases <- as.character(casesBio_Interest)
  ActualControls <- lapply(selCases, function(charcase)
  {
    caseset00               <- Biocrates_Normalized %>% filter(Study == StudyInterest) %>% filter(Cncr_Caco == 1, Idepic == charcase) %>% dplyr::select(CaseSet) %>% pull
    if (length(caseset00)>1){print(icase, "\t", Dataindexcase$Idepic, "\t", caseset00, "\n")} # 
    idcaseset00             <- Biocrates_Normalized %>% filter(Study == StudyInterest) %>% filter(CaseSet == caseset00) %>% dplyr::select(Idepic) %>% pull
    return(Biocrates_Normalized %>% filter(Study == StudyInterest) %>% filter(Idepic %in% idcaseset00, Idepic != charcase) %>% dplyr::select(Idepic) %>% pull %>% as.character)
  })
  
  FlexMatchVarTechnical <- FlexMatchVarTechnical[FlexMatchVarTechnical%in%FlexMatchingInterest]
  if (length(FlexMatchVarTechnical)==0){FlexMatchVarTechnical = NULL}
  
  forKMweights <- EligiblePopulationInterestBiocrat_VV %>% mutate(across(ends_with("_Bld"), as.numeric), 
                                                                  Fasting_C = as.numeric(Fasting_C), 
                                                                  Phase_Mnscycle_3 = case_when(Phase_Mnscycle %in% c(1, 2)~1,
                                                                                               Phase_Mnscycle %in% c(4,5,6)~3, 
                                                                                               !is.na(Phase_Mnscycle)~2), 
                                                                  Phase_Mnscycle = as.numeric(Phase_Mnscycle), 
                                                                  Length_Bld2015 = LengthVV) %>%
                                                           replace_na(list(Menop_Bld = 999, Phrt_Bld= 999)) #, Phase_Mnscycle=999
  
  forStudyInterest       <- computation_risk_sets_foractualNCCs(forKMweights, selCases, ActualControls, 
                                                                StrictMatchVar=StrictMatchingInterest, FlexMatchVar=FlexMatchingInterest, 
                                                                listTol=listToleranceInterest, FlexMatchVarTechnical=FlexMatchVarTechnicalInterest)
  
  PselCases <- forKMweights  %>% filter(IndFac==1)  %>% mutate(forpsel = 1*(Sel_all>=2)) %>% select(forpsel) %>% pull %>% mean
  NCC       <- Biocrates_Normalized %>% filter(Study == StudyInterest) 
  
  km_matched <- compute_weights_KM_Untypical(forKMweights, NCC, pselCases = PselCases,
                               ListRiskSets = forStudyInterest$EligibleSets, 
                               VarIdCohort = "Idepic", VarIndCaseCohort = "IndFac", VarStatforWeight = "Sel_Interest", 
                               VarIdNCC = "Idepic", VarIndCaseNCC = "Cncr_Caco")
  
  forweights_all <- forweights_all %>% left_join(km_matched %>% rename(Idepic = Id))
  
  if (!is.null(FlexMatchVarTechnicalInterest))
  {
    km_matched_woTek <- compute_weights_KM_Untypical(forKMweights, NCC, pselCases = PselCases,
                                             ListRiskSets = forStudyInterest$EligibleSetswoTechnical, 
                                             VarIdCohort = "Idepic", VarIndCaseCohort = "IndFac", VarStatforWeight = "Sel_Interest", 
                                             VarIdNCC = "Idepic", VarIndCaseNCC = "Cncr_Caco")
  
    forweights_all <- forweights_all %>% left_join(km_matched_woTek %>% rename(Idepic = Id, WeightsKM_woTek= WeightsKM, ProbInclKM_woTek=ProbInclKM))
  }
  
  forweights_all <- forweights_all %>%  mutate(Weights_FullCohort = 1)
  
  # InterestwithProbs      <- forProbInclu(selCases, ActualControls= ActualControls, probsel=1, 
  #                                        EligibleSets= forStudyInterest$EligibleSets, NbmatchedCtrls= rep(1, length(selCases)), 
  #                                        EligibleSetswoTek = forStudyInterest$EligibleSetswoTechnical)
  
  # saveRDS(InterestwithProbs, file= paste0("RES/DataWithWeights_Biocrates_", StudyInterest,"_2025.rds"))
  
  all_eligible_KM       <- unique(unlist(forStudyInterest$EligibleSets))
  all_eligible_KM_woTek <- unique(unlist(forStudyInterest$EligibleSetswoTechnical))
  length(all_eligible_KM); length(all_eligible_KM_woTek)
  
  NoIn_woTek <- EligiblePopulationInterestBiocrat_VV %>% filter(!Idepic %in% all_eligible_KM_woTek) %>%
                                                         select(Idepic, all_of(c(FlexMatchingInterest[1:7], StrictMatchingInterest)), LengthVV)

  comprehensive_res <- list(
    forweights_all        = forweights_all, 
    all_eligible_KM       = all_eligible_KM, 
    all_eligible_KM_woTek = all_eligible_KM_woTek
  )
  return(comprehensive_res)
}



fun_KM_estimate <- function(forweights_all, StudyInterest, WeightsToKeep = c("WeightsGAM_NoTek", "WeightsGAM_WithTek", "WeightsKM", "WeightsKM_woTek", "Weights_FullCohort", "Weights_Selectable"), 
                            levelsWeights = c("Weights_FullCohort", "Weights_Selectable", "WeightsGAM_WithTek", "WeightsGAM_NoTek", "WeightsKM", "WeightsKM_woTek"), 
                            legendWeights = c("Originating Cohort", "Eligible", paste0("GAM-Weights", c("", " wo Tech.")), paste0("KM-Weights", c("", " wo Tech."))))
{
  ####### Survival
  cohort_all <- forweights_all %>% pivot_longer(cols = which(grepl("Weight", colnames(forweights_all))), names_to = "TypeAnalysis", values_to = "Weights")
  cohort_all <- cohort_all %>% filter(TypeAnalysis %in% WeightsToKeep) %>%
                               mutate(TypeAnalysis = as_factor(TypeAnalysis))                
  cohort_all$TypeAnalysis = factor(cohort_all$TypeAnalysis, levels=levelsWeights)
    #"WeightsGAM_all", 
  fitKM      <- survfit(Surv(LengthVV/365.25, 1*(IndFac==1)) ~ TypeAnalysis, weights = Weights, data=cohort_all) 
  
  keep <- c("Weights_FullCohort","Weights_Selectable","WeightsGAM_WithTek","WeightsKM")
  labels_keep <- c("Originating Cohort","Eligible","GAM-Weights","KM-Weights")
  lt          <- c("solid","longdash","dashed","dotted")
  cols <- c("#F8766D", "#7CAE00", "#007C80", "#C77CFF")

  cohort_plot <- droplevels(cohort_all[cohort_all$TypeAnalysis %in% keep, ])
  cohort_plot$TypeAnalysis <- factor(cohort_plot$TypeAnalysis, levels = keep)

  fitKM_plot <- survfit(Surv(LengthVV/365.25, 1*(IndFac==1)) ~ TypeAnalysis,
                        weights = Weights, data = cohort_plot)

  kmplot <- ggsurvplot(
    fitKM_plot, data = cohort_plot, ylim = c(0.55, 1), censor = FALSE, conf.int = TRUE,
    palette = cols,                 # <- colors applied
    linetype = lt,                  # <- linetypes applied
    size = 1.1,
    legend.title = "", legend.labs = labels_keep,
    xlab = "Time (years)"
  )

  kmplot$plot <- kmplot$plot +
    guides(color = guide_legend(override.aes = list(linetype = lt, size = 1.3)),
          linetype = "none")


  kmplot
  
  cohortred              <- cohort_all %>% filter(!grepl("oTek", TypeAnalysis))
  cohortred$TypeAnalysis = factor(cohortred$TypeAnalysis, levels=levelsWeights[which(!grepl("oTek", levelsWeights))])
  
  fitKMred      <- survfit(Surv(LengthVV/365.25, 1*(IndFac==1)) ~ TypeAnalysis, weights = Weights, data=cohortred) 
  kmplotred     <- ggsurvplot(fitKMred, data=cohort_all , ylim =c(0.55, 1), censor=FALSE, conf.int = T, 
                              legend.title = "",
                              size = "strata",
                              legend.labs = legendWeights[which(!grepl("wo Tech", legendWeights))],
                              xlab = "Time (years)")
  kmplotred$plot <- kmplotred$plot + 
    scale_size_manual(values = c(2, rep(1, length(levels(cohortred$TypeAnalysis))-1)), guide = "none") +
    scale_linetype_manual(values = c("solid","dashed","dotdash","dotted","longdash","twodash"))
                      # scale_color_manual(values = c(hue_pal()(2)[2], hue_pal()(2)[1])) + 
                      # scale_fill_manual(values = c(hue_pal()(2)[2], hue_pal()(2)[1]))
                      # 
  ggsave(paste0("Figures/KM_", StudyInterest,"_June2025.png"),
         plot=print(kmplot$plot), scale = 1, width = 8, height= 5, units = c("in"), dpi = 300)
  
  ggsave(paste0("Figures/KM_red_", StudyInterest,"_June2025.png"),
         plot=print(kmplotred$plot), scale = 1, width = 8, height= 5, units = c("in"), dpi = 300)
  
  # ggsave(paste0("Figures/KM_", StudyInterest,"_June2025.eps"), device="eps",
  #        plot=print(kmplot$plot), scale = 1, width = 8, height= 5, units = c("in"), dpi = 300)

  grDevices::cairo_ps(
    filename = paste0("Figures/KM_", StudyInterest, "_June2025.eps"),
    width = 8, height = 5, pointsize = 12,
    onefile = FALSE,
    fallback_resolution = 800
  )
  print(kmplot$plot)
  dev.off()

  grDevices::cairo_ps(
    filename = paste0("Figures/KM_red_", StudyInterest, "_June2025.eps"),
    width = 8, height = 5, pointsize = 12,
    onefile = FALSE,
    fallback_resolution = 800
  )
  print(kmplotred$plot)
  dev.off()
  
  saveRDS(fitKM, file=paste0("RES/KM_", StudyInterest,"_June2025.rds"))
  saveRDS(fitKMred, file=paste0("RES/KM_red_", StudyInterest,"_June2025.rds"))
}

fun_HR_std <- function(forweights_all, EligiblePopulationInterestBiocrat_VV, StudyInterest, AdjCLogit, AdjWeightedCox, 
                       Weightstotest = c("WeightsGAM_WithTek", "WeightsGAM_NoTek", "WeightsKM", "WeightsKM_woTek"))
{
  if(AdjCLogit[1] != AdjWeightedCox[1]){stop("AdjCLogit[1] should be the same as AdjWeightedCox[1]")}
  ### HR for the "standard risk factors"
  formula.clogit                       <- paste0("Cncr_Caco ~ ", paste0(AdjCLogit, collapse = " + "), 
                                                 "+ strata(CaseSet)")
  formula.weightedCox                  <- paste0("Surv(Age_maxBldRecr, Age_Exit, 1*(IndFac==1)) ~ ",  
                                                 paste0(AdjWeightedCox, collapse = " + "),  
                                                 " + strata(Center)")
  
  Recap_Asso_Bmi                       <- tibble(RiskFactor= character(), HR = numeric(), SD = numeric(), pval = numeric(), Analysis = character())
  ### C-logit
  NCC_Interest    <- Biocrates_Normalized %>% 
                        filter(Study ==StudyInterest) %>% 
                        mutate(Menop_BldFac = as_factor(case_when(
                                                Sex==1~0,
                                                TRUE~Menop_Bld)), 
                              Phrt_BldFac   = as_factor(replace_na(case_when(
                                                  Sex==1~0, 
                                                  TRUE~Phrt_Bld), 99)), 
                              Fasting_Fac              = as_factor(replace_na(Fasting_C, 99)),
                              
                              
                              TbldMissing              = 1*(is.na(T_Bld_Coll)), 
                              forFourierTbld           = case_when(is.na(T_Bld_Coll)~0, T~minute(T_Bld_Coll) + 60 *hour(T_Bld_Coll)), 
                              FourierTblood1           = sin(2*pi*forFourierTbld/(24*60)),
                              FourierTblood2           = cos(2*pi*forFourierTbld/(24*60)),
                              FourierTblood3           = sin(4*pi*forFourierTbld/(24*60)),
                              FourierTblood4           = cos(4*pi*forFourierTbld/(24*60)),
                              DbldMissing              = 1*(is.na(D_Bld_Coll)),
                              dayofyear_bdcol          = case_when(is.na(D_Bld_Coll)~0, T~yday(D_Bld_Coll)), 
                              FourierDblood1           = sin(2*pi*dayofyear_bdcol/365),
                              FourierDblood2           = cos(2*pi*dayofyear_bdcol/365),
                              FourierDblood3           = sin(4*pi*dayofyear_bdcol/365),
                              FourierDblood4           = cos(4*pi*dayofyear_bdcol/365),
                              forTBldVV                = T_Bld_Coll,
                              Menop_BldNum             = as.numeric(as.character(Menop_Bld)), 
                              Phase_MnsCyVV            = as_factor(replace_na(ifelse(Menop_Bld %in% c(1, 3), 0, Phase_Mnscycle), 99)),
                              Age_maxBldRecr           = pmax(Age_Blood, Age_Recr))
  
  mod_clogit      <- clogit(as.formula(formula.clogit), data= NCC_Interest)
  
  hrstemp   <- exp(c(summary(mod_clogit)$coef[1, 1]))
  sdstemp   <- c(summary(mod_clogit)$coef[1, 3]) #, summary(modCox_WeightsKM_woTek)$coef[1, 4])
  pvalstemp <- c(summary(mod_clogit)$coef[1, 5]) #, summary(modCox_WeightsKM_woTek)$coef[1, 6])
  risktemp  <- c(AdjCLogit[1])
  analytemp <- c("cLogit") 
  
  for (typewe in c("FullCohort", Weightstotest))
  {
    indcolsd <- 4; indcolpval <- 6
    if(!grepl("GAM_all", typewe))
    {
      if(!grepl("FullCohort", typewe)) 
      {
        expr_WeightedCox <- paste0("modCoxWeighted <- coxph(as.formula(formula.weightedCox), weights=", typewe, ", data = subset(forweights_all, Sel_Interest>0))")
        
      }else{
        expr_WeightedCox <- "modCoxWeighted <- coxph(as.formula(formula.weightedCox), data = EligiblePopulationInterestBiocrat_VV)"
        indcolsd <- 3; indcolpval <- 5  # for unweighted models, robust variance is not present in the output of the model :)
      }
    }else{
      expr_WeightedCox <- paste0("modCoxWeighted <- coxph(as.formula(formula.weightedCox), weights=", typewe, ", data = subset(forweights_all, Sel_all>0))")
    }
    eval(parse(text = expr_WeightedCox))
    
    hrstemp   <- c(hrstemp, exp(summary(modCoxWeighted)$coef[1, 1]))
    sdstemp   <- c(sdstemp, summary(modCoxWeighted)$coef[1, indcolsd]) 
    pvalstemp <- c(pvalstemp, summary(modCoxWeighted)$coef[1, indcolpval]) 
    risktemp  <- c(risktemp, AdjCLogit[1])
    analytemp <- c(analytemp, typewe)
  }
  
  Recap_Asso_Bmi <- bind_rows(Recap_Asso_Bmi, tibble(RiskFactor= risktemp, HR = hrstemp, SD = sdstemp, pval = pvalstemp, Analysis = analytemp))

}


fun_HR_Metabolites <- function(forweights_all, IndexMets, StudyInterest, AdjCLogit, AdjWeightedCox, 
                               Weightstotest = c("WeightsGAM_WithTek", "WeightsGAM_NoTek", "WeightsKM", "WeightsKM_woTek"), 
                               WeightsOrder = c("WeightsKM", "WeightsKM_woTek", "WeightsGAM_NoTek", "WeightsGAM_WithTek"), 
                               WeightsToPrint = c("KM weights", "KM weights wo Technical M", "GAM weights", "GAM weights wo Technical M"))
{
  AllMets          <- colnames(Biocrates_Normalized)[IndexMets]
  #### PCA for ENT...
  scaleVV <- function(x){return( (x-mean(x))/sd(x))}
  MatrixMeta <- forweights_all %>% filter(Sel_all!=0) %>% select(all_of(AllMets))
  MatrixMeta <- MatrixMeta %>% mutate_all(scaleVV)
  resPCA     <- PCA(MatrixMeta, graph=F)
  saveRDS(PCA, file = paste0("RES/PCA_", StudyInterest, "_June2025.rds"))
  
  

  
  
  ######### HR of each Metabolite in StudyInterest
  Recap_Asso                     <- tibble(Metabolite= character(), HR = numeric(), SD = numeric(), pval = numeric(), Analysis = character())
  foranalysisweigthed            <- forweights_all
  Biocrates_Normalized           <- Biocrates_Normalized %>% mutate(Menop_BldFac = as_factor(case_when(
                                                                                    Sex==1~0,
                                                                                    TRUE~Menop_Bld)), 
                                                                                    Phrt_BldFac   = as_factor(replace_na(case_when(
                                                                                      Sex==1~0, 
                                                                                      TRUE~Phrt_Bld), 99)), 
                                                                    Fasting_Fac              = as_factor(replace_na(Fasting_C, 99)),
                                                                    
                                                                    
                                                                    TbldMissing              = 1*(is.na(T_Bld_Coll)), 
                                                                    forFourierTbld           = case_when(is.na(T_Bld_Coll)~0, T~minute(T_Bld_Coll) + 60 *hour(T_Bld_Coll)), 
                                                                    FourierTblood1           = sin(2*pi*forFourierTbld/(24*60)),
                                                                    FourierTblood2           = cos(2*pi*forFourierTbld/(24*60)),
                                                                    FourierTblood3           = sin(4*pi*forFourierTbld/(24*60)),
                                                                    FourierTblood4           = cos(4*pi*forFourierTbld/(24*60)),
                                                                    DbldMissing              = 1*(is.na(D_Bld_Coll)),
                                                                    dayofyear_bdcol          = case_when(is.na(D_Bld_Coll)~0, T~yday(D_Bld_Coll)), 
                                                                    FourierDblood1           = sin(2*pi*dayofyear_bdcol/365),
                                                                    FourierDblood2           = cos(2*pi*dayofyear_bdcol/365),
                                                                    FourierDblood3           = sin(4*pi*dayofyear_bdcol/365),
                                                                    FourierDblood4           = cos(4*pi*dayofyear_bdcol/365),
                                                                    forTBldVV                = T_Bld_Coll,
                                                                    Menop_BldNum             = as.numeric(as.character(Menop_Bld)), 
                                                                    Phase_MnsCyVV            = as_factor(replace_na(ifelse(Menop_Bld %in% c(1, 3), 0, Phase_Mnscycle), 99)),
                                                                    Age_maxBldRecr           = pmax(Age_Blood, Age_Recr))
  for (individualmet in AllMets)
  {
    cat("******************\t", individualmet, "\n")
    
    formula.clogit                       <- paste0("Cncr_Caco ~ ", individualmet, " + ", paste0(AdjCLogit, collapse = " + "), 
                                                      "+ strata(CaseSet)")
    formula.weightedCox                  <- paste0("Surv(Age_maxBldRecr, Age_Exit, 1*(IndFac==1)) ~ ", individualmet,  " + ", 
                                                   paste0(AdjWeightedCox, collapse = " + "),  
                                                   " + strata(Center)")
    
    sdmet                                <- sqrt(wtd.var(forweights_all[, individualmet], weights = forweights_all$WeightsGAM_WithTek,na.rm=T)) ### WeightsGAM_int would give very similar results...
    foranalysisweigthed[, individualmet] <- foranalysisweigthed[, individualmet]/sdmet
    
    ### C-logit
    NCC_Interest    <- Biocrates_Normalized %>% filter(Study ==StudyInterest) 
    
    NCC_Interest[, individualmet]      <- NCC_Interest[, individualmet]/sdmet
    mod_clogit      <- clogit(as.formula(formula.clogit), data= NCC_Interest)
    
    hrstemp   <- exp(c(summary(mod_clogit)$coef[1, 1]))
    sdstemp   <- c(summary(mod_clogit)$coef[1, 3]) #, summary(modCox_WeightsKM_woTek)$coef[1, 4])
    pvalstemp <- c(summary(mod_clogit)$coef[1, 5]) #, summary(modCox_WeightsKM_woTek)$coef[1, 6])
    metstemp  <- c(individualmet)
    analytemp <- c("cLogit") 
    
    for (typewe in Weightstotest)
    {
      if(!grepl("GAM_all", typewe))
      {
        expr_WeightedCox <- paste0("modCoxWeighted <- coxph(as.formula(formula.weightedCox), weights=", typewe, ", data = subset(foranalysisweigthed, Sel_Interest>0))")
      }else{
        expr_WeightedCox <- paste0("modCoxWeighted <- coxph(as.formula(formula.weightedCox), weights=", typewe, ", data = subset(foranalysisweigthed, Sel_all>0))")
      }
      eval(parse(text = expr_WeightedCox))
      
      hrstemp   <- c(hrstemp, exp(summary(modCoxWeighted)$coef[1, 1]))
      sdstemp   <- c(sdstemp, summary(modCoxWeighted)$coef[1, 4]) 
      pvalstemp <- c(pvalstemp, summary(modCoxWeighted)$coef[1, 6]) 
      metstemp  <- c(metstemp, individualmet)
      analytemp <- c(analytemp, typewe)
    }

    Recap_Asso <- bind_rows(Recap_Asso, tibble(Metabolite= metstemp, HR = hrstemp, SD = sdstemp, pval = pvalstemp, Analysis = analytemp))
    
  }
  
  Recap_Asso_wide <- Recap_Asso %>% pivot_wider(names_from = "Analysis", values_from = c("HR", "pval", "SD"))
  Recap_Asso_wide <- Recap_Asso_wide %>% mutate(# FDR_cLogit = p.adjust(pval_cLogit, method = "BH"),
                                                across(starts_with("pval_"), ~ p.adjust(.x, , method = "BH"), .names = "{sub('pval', 'FDR', col)}"))
                                                # FDR_Weighted_Interest = p.adjust(pval_Weighted_Interest),
                                                # FDR_Weighted_all = p.adjust(pval_Weighted_all, method = "BH"),
                                                # FDR_WeightedKM = p.adjust(pval_WeightedKM, method = "BH")) #,
                                                #FDR_WeightedKM_woTek = p.adjust(pval_WeightedKM_woTek, method = "BH"))
  
  
  saveRDS(Recap_Asso_wide, file= paste0("RES/Recap_Asso_wide_", StudyInterest, "_June2025.rds"))
  
  Recap_HRs <- Recap_Asso_wide %>% select(Metabolite) %>%
                bind_cols(Recap_Asso_wide %>% select_if(grepl("HR_", names(.))) %>%
                              mutate_all(log)) %>%
                pivot_longer(-c(1,2), names_to = "Weights", values_to = "logHR") %>%
                rename(logHR_cLogit = 2) %>%
                mutate(Weights = str_replace(Weights, "HR_", ""))
  Recap_FDRs <- Recap_Asso_wide %>% select(Metabolite) %>%
                    bind_cols(Recap_Asso_wide %>% select_if(grepl("FDR_", names(.)))) %>%
                    pivot_longer(-c(1,2), names_to = "Weights", values_to = "FDR") %>%
                    rename(FDR_cLogit = 2)%>%
                     mutate(Weights = str_replace(Weights, "FDR_", ""))
  
  Recap_all <- Recap_HRs %>% left_join(Recap_FDRs) %>% 
                  mutate(Signif = case_when(
                      FDR_cLogit < 0.05 & FDR < 0.05 ~ "Both",
                      FDR_cLogit >= 0.05 & FDR >= 0.05 ~ "None", 
                      FDR_cLogit >= 0.05 & FDR < 0.05 ~ "Only Weighted",
                      FDR_cLogit < 0.05 & FDR >= 0.05 ~ "Only cLogit",
                  ))
  
  Recap_all <- Recap_all %>% 
                  mutate(Weights = as_factor(Weights)) %>%
                  mutate(Weights = factor(Weights, levels = WeightsOrder))
  
  suppWeights <- WeightsToPrint
  names(suppWeights) <- WeightsOrder
  
  # PlotHRs <- ggplot(Recap_HRs, aes(x=logHR_cLogit, y=logHR, colour=Weights)) + 
  #   geom_point() +
  #   geom_abline(intercept = 0, slope=1) +  
  #   geom_smooth(method = "loess", se = TRUE) + 
  #   facet_wrap(~Weights, labeller = labeller(Weights = suppWeights))+ 
  #   xlab("Log HR (cond. logistic)") + ylab("Log HR (Weighted analysis)") + 
  #   theme(legend.position = "bottom") + 
  #   scale_colour_discrete(name = "", labels = c("GAM-Weights ALL", "GAM-Weights ENDO", "KM-Weights ENDO"))
  # 
  # PlotHRs <- ggplot(Recap_all, aes(x=logHR_cLogit, y=logHR, colour=Signif)) + 
  #   geom_point() +
  #   geom_abline(intercept = 0, slope=1) +  
  #   #geom_smooth(method = "loess", se = TRUE) + 
  #   facet_wrap(~Weights) + #, labeller = labeller(Weights = suppWeights)) + 
  #   xlab("Log HR (cond. logistic)") + ylab("Log HR (Weighted analysis)") + 
  #   theme(legend.position = "bottom") # + 
  #   #scale_colour_discrete(name = "", labels = c("GAM-Weights ALL", "GAM-Weights ENDO", "KM-Weights ENDO"))
  # ggsave(paste0("Figures/HRs_", StudyInterest,"_June2025.png"),
  #        plot=PlotHRs, scale = 1, width = 14, height= 5, units = c("in"), dpi = 300)
    

  bias_km         <- mean((log(Recap_Asso_wide$HR_cLogit) - log(Recap_Asso_wide$HR_WeightsKM)))
  bias_km_notek   <- mean((log(Recap_Asso_wide$HR_cLogit) - log(Recap_Asso_wide$HR_WeightsKM_woTek)))
  bias_gam        <- mean((log(Recap_Asso_wide$HR_cLogit) - log(Recap_Asso_wide$HR_WeightsGAM_NoTek)))
  bias_gam_notek  <- mean((log(Recap_Asso_wide$HR_cLogit) - log(Recap_Asso_wide$HR_WeightsGAM_WithTek)))
  
  rmse_km         <- format(mean((log(Recap_Asso_wide$HR_cLogit) - log(Recap_Asso_wide$HR_WeightsKM))^2), scientific = T, digits = 3)
  rmse_km_notek   <- format(mean((log(Recap_Asso_wide$HR_cLogit) - log(Recap_Asso_wide$HR_WeightsKM_woTek))^2), scientific = T, digits = 3)
  rmse_gam        <- format(mean((log(Recap_Asso_wide$HR_cLogit) - log(Recap_Asso_wide$HR_WeightsGAM_NoTek))^2), scientific = T, digits = 3)
  rmse_gam_notek  <- format(mean((log(Recap_Asso_wide$HR_cLogit) - log(Recap_Asso_wide$HR_WeightsGAM_WithTek))^2), scientific = T, digits = 3)
  
  for_ggtext      <- tibble(Text    = paste0("    RMSE = ", c(rmse_km, rmse_km_notek, rmse_gam, rmse_gam_notek)), 
                            Weights = WeightsOrder, 
                            x       = -0.15, #min(Recap_all$logHR_cLogit), 
                            y       = 0.29) #max(Recap_all$logHR)
    
  PlotHRs_onlyInt <- ggplot(Recap_all %>% filter(!grepl("_all", Weights))) + #Recap_all
    geom_ellipse(aes(x0 = logHR_cLogit, y0=logHR, a=1.96*SD_cLogit, b = 1.96*SD, angle=0), fill="grey90", alpha=0.1, colour="grey90") + 
    geom_point(aes(x=logHR_cLogit, y=logHR)) +
    geom_abline(intercept = 0, slope=1, col="grey30") +  
    #geom_smooth(method = "loess", se = TRUE) + 
    facet_wrap(~Weights, labeller = labeller(Weights = suppWeights), nrow=1) + 
    geom_text(data= for_ggtext, aes(x=x, y=y, label=Text)) +
    xlab("Log HR (conditional logistic reg.)") + ylab("Log HR (Weighted analysis)") + 
    theme_minimal() + 
    coord_fixed()
  #theme(legend.position = "bottom") # + 
  

  
  ggsave(paste0("Figures/HRs_", StudyInterest,"_noAll_June2025.png"),
         plot=PlotHRs_onlyInt, scale = 1, width = 18, height= 4, units = c("in"), dpi = 300)
  
  return(Recap_Asso_wide)
}


fun_Correlation_Exposures <- function(Eligible_withWeights, 
                                      fullEPICData, 
                                      Vec_Cont_Variables = c("Age_Recr", "QE_ENERGY", "Bmi_C", "Whr_C", 
                                                             "MRMed_Score", "Alc_Re", "Height_C", "Hli_Smoke", 
                                                             "Pa_Work", "Pa_Vig", "Pa_Recr", "Pa_Mets_C"), 
                                      WhichAnalyses = c("_FullCohort", "NCC_int_Total", "NCC_int_Ctls", "KM", "KM_woTek", "GAM_NoTek", "GAM_WithTek", "GAM_thresh"))
  
  # "_FullCohort", "KM", "GAM_int", "GAM_all", "GAM_int_thresh", "GAM_all_thresh", "GAM_int_noTD", "GAM_all_noTD", "NCC_int_Total", "NCC_int_Ctls", "NCC_all_Total", "NCC_all_Ctls")
{
  if (mean(Eligible_withWeights$Idepic %in% fullEPICData$Idepic)<1){stop("The Full EPIC Data is supposed to include all the particiapnts of the Eligible sub-study...")}
  red_cohort <- fullEPICData %>% select(Idepic, all_of(Vec_Cont_Variables))
  red_eligib <- Eligible_withWeights %>% select_if(names(.) =="Idepic" | grepl("Weight", names(.)))
  red_forcor <- red_cohort %>% inner_join(red_eligib) %>%
                               mutate(WeightsGAM_thresh = pmin(WeightsGAM_NoTek, 100))#, 
                               #        WeightsGAM_all_thresh = pmin(WeightsGAM_all, 100))
  red_forcor <- red_forcor %>% modify_if(names(.)%in%Vec_Cont_Variables, ~scale(.))
  colSums(is.na(red_forcor %>% select(all_of(Vec_Cont_Variables))))
  
  Res_all_as <- tibble(Var1 = character(), Var2= character(), Analysis = character(), CoefLM = double(), SD =double())
  for (typeAnalysis in WhichAnalyses)
  {
    eval(parse(text=paste0("MatCor_", typeAnalysis, "<- matrix(1, length(Vec_Cont_Variables), length(Vec_Cont_Variables))")))
    for (ivar1 in 1:(length(Vec_Cont_Variables)-1))
    {
      for (ivar2 in ((ivar1+1):length(Vec_Cont_Variables)))
      {
        var1      <- Vec_Cont_Variables[ivar1] 
        var2      <- Vec_Cont_Variables[ivar2]

        spec_weights <- paste0("Weights", typeAnalysis)
        tempdata     <- red_forcor %>% 
                              mutate(vvar1= get(var1), vvar2= get(var2)) %>%
                              filter(spec_weights > 0 & ((!is.na(vvar1))*(!is.na(vvar2))==1))

        Covar1       <- tempdata %>% select(any_of(var1)) %>% pull
        Covar2       <- tempdata %>% select(any_of(var2)) %>% pull
        weightsss    <- tempdata %>% select(any_of(spec_weights)) %>% pull
        #wCorrsss     <- weightedCorr(Covar1, Covar2, weights = weightsss, method ="Pearson")
        
        wCorrsss     <- wtd.cor(x= Covar1, y=Covar2, weight=weightsss, mean1=F, collapse=TRUE, bootse=T,
                                    bootp=FALSE, bootn=1000)
        wcorrrr      <- wCorrsss[1, "correlation"]
        sdwcor       <- wCorrsss[1, "std.err"]
        # mod_lm     <- lm(get(var1)~ get(var2), weights = get(paste0("Weights", typeAnalysis)), data= subset(red_forcor, get(paste0("Weights", typeAnalysis))>0))        
        temptib      <- tibble(Var1 = var1, Var2= var2, Analysis = typeAnalysis, CoefLM =wcorrrr, SD =sdwcor)
        Res_all_as   <- bind_rows(Res_all_as, temptib)
        eval(parse(text=paste0("MatCor_", typeAnalysis, "[ivar1, ivar2] <- MatCor_", typeAnalysis, "[ivar2, ivar1] <- wcorrrr")))
      }
    }
    ### Now set the upper tri. part of the Correlation matrix to that of FullCohort (for ease of visual comparison)
    eval(parse(text = paste0("MatCor_", typeAnalysis, "[upper.tri(MatCor_", typeAnalysis,")] <- MatCor__FullCohort[upper.tri(MatCor__FullCohort)]")))
  }
  Res_wide <- Res_all_as %>% 
                pivot_wider(names_from = Analysis, values_from = c(CoefLM, SD)) %>%
                rename(FullCohort = 3)
  saveRDS(Res_wide, file=paste0("RES/MatCorrExpossures_", StudyInterest,"_June2025_ForTableRossella.rds"))
  Res_wide  <- Res_wide %>% 
                    select_if(grepl("Var|Full|Coef", names(.)) & !(grepl("SD", names(.)))) 
  Res_wide2 <- Res_wide %>%
                mutate(across(4: (2 + length(WhichAnalyses)), ~ .x - FullCohort))
  colnames(Res_wide2)[4: (2 + length(WhichAnalyses))]<- paste0("Diff", colnames(Res_wide2)[4: (2 + length(WhichAnalyses))])
  Res_wide  <- Res_wide %>% left_join(Res_wide2)

  saveRDS(Res_wide, file=paste0("RES/MatCorrExpossures_", StudyInterest,"_June2025.rds"))
  return(Res_wide)
  
  # par(mfrow=c(4,2))
  # 
  # corrplot(MatCor__FullCohort)
  # corrplot(MatCor_KM)
  # corrplot(MatCor_GAM_int)
  # corrplot(MatCor_GAM_all)
  # corrplot(MatCor_NCC_all_Total)
  # corrplot(MatCor_NCC_all_Ctls)
  # corrplot(MatCor_NCC_int_Total)
  # corrplot(MatCor_NCC_int_Ctls)
  # 
  # 
  # sum(Res_wide$DiffKM^2)
  # sum(Res_wide$DiffGAM_int^2)
  # sum(Res_wide$DiffGAM_all^2)
  # sum(Res_wide$DiffNCC_int_Total^2)
  # sum(Res_wide$DiffNCC_all_Total^2)
  # sum(Res_wide$DiffNCC_int_Ctls^2)
  # sum(Res_wide$DiffNCC_all_Ctls^2)
  # 
  # 
  # max(Res_wide$DiffKM^2)
  # max(Res_wide$DiffGAM_int^2)
  # max(Res_wide$DiffGAM_all^2)
  # max(Res_wide$DiffNCC_int_Total^2)
  # max(Res_wide$DiffNCC_all_Total^2)
  # max(Res_wide$DiffNCC_int_Ctls^2)
  # max(Res_wide$DiffNCC_all_Ctls^2)
  
}


fun_Corr_Exposures_Metabo <- function(Eligible_withWeights, 
                                      IndexMets, 
                                      Vec_Cont_Variables = c("Bmi_C"), 
                                      WhichAnalyses = c("NCC_int_Total", "NCC_int_Ctls", "KM", "KM_woTek", "GAM_NoTek", "GAM_WithTek", "GAM_thresh"))
{
  Vec_names_Metabo <- colnames(Eligible_withWeights)[IndexMets]
  red_forcor       <- Eligible_withWeights %>% 
                        modify_if(names(.)%in%c(Vec_names_Metabo, Vec_Cont_Variables), ~scale(.)) %>%
                        mutate(WeightsGAM_thresh = pmin(WeightsGAM_NoTek, 100)) #, 
                               # WeightsGAM_all_thresh = pmin(WeightsGAM_all, 100))
  
  
  Res_all_as <- tibble(Var1 = character(), Var2= character(), Analysis = character(), CoefLM = double())
  AllMat     <- NULL 
  for (typeAnalysis in WhichAnalyses)
  {
    eval(parse(text=paste0("MatCor_", typeAnalysis, "<- matrix(1, length(IndexMets), length(Vec_Cont_Variables))")))
    for (ivar1 in 1:length(IndexMets))
    {
      for (ivar2 in 1:length(Vec_Cont_Variables))
      {
        
        var1      <- Vec_names_Metabo[ivar1] 
        var2      <- Vec_Cont_Variables[ivar2]
        
        spec_weights <- paste0("Weights", typeAnalysis)
        tempdata     <- red_forcor %>% 
          mutate(vvar1= get(var1), vvar2= get(var2)) %>%
          filter(spec_weights > 0 & (!is.na(vvar1))*(!is.na(vvar2))==1)
        
        Covar1       <- tempdata %>% select(any_of(var1)) %>% pull
        Covar2       <- tempdata %>% select(any_of(var2)) %>% pull
        weightsss    <- tempdata %>% select(any_of(spec_weights)) %>% pull
        wCorrsss     <- weightedCorr(Covar1, Covar2, weights = weightsss, method ="Pearson")
        # mod_lm     <- lm(get(var1)~ get(var2), weights = get(paste0("Weights", typeAnalysis)), data= subset(red_forcor, get(paste0("Weights", typeAnalysis))>0))        
        temptib    <- tibble(Var1 = var1, Var2= var2, Analysis = typeAnalysis, CoefLM =wCorrsss)
        Res_all_as <- bind_rows(Res_all_as, temptib)
        eval(parse(text=paste0("MatCor_", typeAnalysis, "[ivar1, ivar2] <- wCorrsss")))
      }
    }
    eval(parse(text = paste0("AllMat <- cbind(AllMat, MatCor_", typeAnalysis, ")")))
  }
  # AllMat <- cbind(MatCor_NCC_int_Ctls, MatCor_NCC_int_Total, MatCor_NCC_all_Ctls, MatCor_NCC_all_Total,
  #                 MatCor_KM, 
  #                 MatCor_GAM_int, MatCor_GAM_all, 
  #                 MatCor_GAM_int_thresh, MatCor_GAM_all_thresh, 
  #                 MatCor_GAM_int_noTD, MatCor_GAM_all_noTD)
  colnames(AllMat) <- WhichAnalyses
  # c( "NCC_int_Ctls", "NCC_int_Total",
  #                        "NCC_all_Ctls", "NCC_all_Total", 
  #                        "KM", "GAM_int", "GAM_all", "GAM_int_thresh", "GAM_all_thresh", 
  #                        "GAM_int_noTD", "GAM_all_noTD")
  AllMat <- as_tibble(AllMat) %>% mutate(Diff_int = abs(NCC_int_Ctls - NCC_int_Total)) #, 
                                         # Diff_all = abs(NCC_all_Ctls - NCC_all_Total))
  
  saveRDS(AllMat, file=paste0("RES/MatCorr_BMI_Metabo_", StudyInterest,"_June2025.rds"))
  return(AllMat)
}
  


##### old... and hopefully not useful anymore
# #forgamweights <- forgamweights %>% mutate(across(c(Alc_Re), ~replace_na(., median(., na.rm=TRUE)))) #%>% left_join(selections %>% select(Idepic, which(grepl("Cncr_Mal", colnames(selections)))))
# 
# formula_GamWeigths_Interest_ctls  <- paste0("indic_inluded_Interest ~ s(LengthVV) + s(Age_Blood) + ", paste0(Adj_GamWeights_Interest, collapse=" + "))
# #formula_GamWeigths_Interest_cases <- paste0("indic_inluded_Interest ~ ", paste0(Adj_GamWeights_Interest, collapse=" + "))
# formula_GamWeigths_all_ctls  <- paste0("indic_inluded_all ~ s(LengthVV) + s(Age_Blood) + ", paste0(Adj_GamWeights_All, collapse=" + "))
# #formula_GamWeigths_all_cases <- paste0("indic_inluded_all ~ ", paste0(Adj_GamWeights_Interest, collapse=" + "))
# 
# modGAM_Interest_red_ctls     <- gam(as.formula(formula_GamWeigths_Interest_ctls), 
#                                     method="REML", family="binomial", data=subset(forgamweights, Sel_Interest!=2)) 
# modGAM_Interest_red_cases    <- gam(I(Sel_Interest==2) ~ s(LengthVV) + CenterFac, method="REML", family="binomial", data=subset(forgamweights, IndFac==1)) # 
# pGAM_Interest_red_ctls       <- as.numeric(predict(modGAM_Interest_red_ctls, newdata = forgamweights, type='response'))
# pGAM_Interest_red_cases      <- as.numeric(predict(modGAM_Interest_red_cases, newdata = forgamweights, type='response'))#mean(1*(forgamweights$Sel_Interest[forgamweights$IndFac==1]==2))#as.numeric(predict(modGAM_Interest_red_cases, type='response'))
# 
# 
# modGAM_all_red_ctls          <- gam(as.formula(formula_GamWeigths_all_ctls), 
#                                     method="REML", family=binomial, data=subset(forgamweights, Sel_all!=2))
# modGAM_all_red_cases         <- gam(I(Sel_all==2) ~ s(LengthVV) + CenterFac, method="REML", family="binomial", data=subset(forgamweights, IndFac==1))
# pGAM_all_red_ctls            <- as.numeric(predict(modGAM_all_red_ctls, newdata = forgamweights, type='response'))
# pGAM_all_red_cases           <- as.numeric(predict(modGAM_all_red_cases, newdata = forgamweights, type='response'))#mean(1*(forgamweights$Sel_all[forgamweights$IndFac==1]==2))# as.numeric(predict(modGAM_all_red_cases, type='response'))
# 
# # forgamweights_cases          <- forgamweights %>% filter(IndFac==1) %>% mutate(Pred_Interest = pGAM_Interest_red_cases, 
# #                                                                                Pred_all      = pGAM_all_red_cases) 
# # 
# # forgamweights_ctls           <- forgamweights %>% filter(IndFac==0) %>% mutate(Pred_Interest = pGAM_Interest_red_ctls, 
# #                                                                                Pred_all      = pGAM_all_red_ctls) 
# 
# # forgamweights_ctls  %>%  group_by(Center) %>% summarise(MeanPred_Interest = mean(Pred_Interest), MeanPred_all = mean(Pred_all), MeanIndic_Interest = mean(indic_inluded_Interest), MeanIndic_all = mean(indic_inluded_all)) %>% print(n=100)
# # forgamweights_cases %>%  group_by(Center) %>% summarise(MeanPred_Interest = mean(Pred_Interest), MeanPred_all = mean(Pred_all), MeanIndic = mean(indic_inluded_Interest)) %>% print(n=100)
# # 
# 
# forgamweights             <- forgamweights %>% mutate(Pred_Interest = case_when(IndFac==0~ pGAM_Interest_red_ctls, 
#                                                                                 TRUE~ pGAM_Interest_red_cases + (1-pGAM_Interest_red_cases)*pGAM_Interest_red_ctls), 
#                                                       Pred_all      = case_when(IndFac==0~pGAM_all_red_ctls, 
#                                                                                 TRUE~ pGAM_all_red_cases + (1-pGAM_all_red_cases)*pGAM_all_red_ctls)) 
# 
# forgamweights_final <-forgamweights %>% # bind_rows(forgamweights_ctls, forgamweights_cases) %>% 
#   mutate(Weight_Interest = case_when(
#     indic_inluded_Interest==0~0, 
#     TRUE~1/Pred_Interest), 
#     Weight_all = case_when(
#       indic_inluded_all==0~0, 
#       TRUE~1/Pred_all) 
#   )

