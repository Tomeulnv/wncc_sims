rm(list=ls())


setwd("/data/Epic/subprojects/BiocratesGWAS/work/WeightedNCC_Final/Simulations/")

library(tidyverse)
library(ggplot2)



remove_first_bit <- Vectorize(function(char){paste0(unlist(strsplit(char, "_"))[-1], collapse="_")})
sizecoho      <- 1e5

##### Simulation Supp Asso: 
##### Illustration of the limitations of weighted analyses
##### when estimating associations between exposures


## Usual  & Inter

pselCases     <- 1 
MAF           <- 0.25
R2SNPBMI      <- R2BMI <- 0.1
logHRMatching <- log(2)


for (Inter in c(0, 1))
{
  RECAP_RES_ASSO <- NULL
  RECAPALL0 <- NULL
  
  for (Collider in c("No collider bias", "Moderate collider bias", "Strong collider bias")) for(Censoring in c("Low", "Moderate", "High")) 
  {
    if (Collider == "No collider bias") {logHRBMI   <- log(1);  logHRSNP   <- log(2)}
    if (Collider == "Moderate collider bias") {logHRBMI   <- log(2);  logHRSNP   <- log(1)}
    if (Collider == "Strong collider bias") {logHRBMI   <- log(2);  logHRSNP   <- log(2)}
    
    if (Censoring == "Low") {agemincens   <- 20;  agemaxcens   <- 50}
    if (Censoring == "Moderate") {agemincens   <- 20;  agemaxcens   <- 30}
    if (Censoring == "High") {agemincens   <- 0;  agemaxcens   <- 30}
    
    if (Inter==0)
    {
      suffsuffix <- paste0("_pselCases_", pselCases, "_ageCens_", agemincens, "_", agemaxcens, "_NNmatching_0_AssoMsnp_1") 
      suffix     <- paste0("Sim_Supp_Asso_CohortSize_", sizecoho, "_MAF_", MAF, "_R2SNPBMI_", R2SNPBMI, "_logHRSNP_", round(logHRSNP,2),
                         "_logHRBMI_", round(logHRBMI,2), "_logHRMatching_", round(logHRMatching,2), 
                         suffsuffix)
      suffixFig  <- "Usual"
    }else{
      suffix <- paste0("Sim_Supp_Asso_Inter_CohortSize_", sizecoho, "_MAF_", MAF, "_R2BMI_", R2BMI, "_logHRSNP_", round(logHRSNP,2),
                       "_logHRBMI_", round(logHRBMI,2), "_logHRMatching_", round(logHRMatching,2), "_pselCases_", pselCases, "_ageCens_", agemincens, "_", agemaxcens)
      suffixFig  <- "Inter"
    }
    namefile <- paste0("Res/Res_", suffix,".rds")
    if (file.exists(namefile))
    {
      recap    <- readRDS(namefile)$RECAP_RES_ASSO
      RECAPALL0 <- rbind(RECAPALL0, recap)
    }
  }
  
  RECAPALL   <- as_tibble(RECAPALL0)%>% rename(agemincens = Agemincens, agemaxcens=Agemaxcens)
  RECAPALL   <- RECAPALL %>% mutate(Collider = as_factor(case_when(logHRBMI==log(1)&logHRSNP==log(2) ~ "No Collider Bias", 
                                                                   logHRBMI==log(2)&logHRSNP==log(1)~ "Moderate Collider Bias", 
                                                                   logHRBMI==log(2)&logHRSNP==log(2)~ "Strong Collider Bias")),
                                    Censoring = as_factor(case_when(agemincens==0 & agemaxcens == 30 ~ "High censoring",
                                                                    agemincens==20 & agemaxcens == 30 ~ "Moderate censoring",
                                                                    agemincens==20 & agemaxcens == 50 ~ "Low censoring")),
                                    across(2:7, ~as_factor(round(.x,2))))
  #levels(RECAPALL$logHRMatching) <- paste0("alpha[M] == ", levels(RECAPALL$logHRMatching))
  
  
  Betas      <- RECAPALL %>% dplyr::select(Censoring, Collider, starts_with(c("beta"))) %>% 
    mutate(across(starts_with("beta"), ~.x - beta_allcohort)) %>% 
    select(-beta_allcohort) %>% 
    gather(Type, beta, starts_with("beta"), factor_key=TRUE) %>% 
    mutate(Typec = as.character(Type), Type = as_factor(remove_first_bit(Typec))) %>% 
    select(-Typec) 
  
  
  # Betas      <- Betas %>% filter(!grepl("VV", Type)) %>% droplevels() 
  ALL_CRIT   <- Betas %>% rename(Value = beta) %>% 
    mutate(Type = fct_recode(Type, "GAM-weights" = "NCC_WeightedGAM_matching", 
                             "KM-weights" = "NCC_WeightedKM_matching", 
                             "GAM-weights thresh." = "NCC_WeightedGAM_matching_thresh", 
                             "GAM-weights no (T,D)" = "NCC_WeightedGAM_noTD", 
                             "Random Subsample" = "allcohort_subsample", 
                             "Unweighted" = "NCC", 
                             "Unweighted (controls only)" = "NCC_ctrl"))
  
  
  ALL_CRIT$Type <- factor(ALL_CRIT$Type, levels = levels(ALL_CRIT$Type)[c(1, 3, 2, 4, 5, 11, 6, 7)]) 
  
  colVVs <- c("Random Subsample" = "gray27", "Unweighted (controls only)" = "seagreen4", "Unweighted" = "palegreen", 
              "KM-weights wo M" = "sienna1", "KM-weights" ="sienna1", 
              "GAM-weights wo M" = "dodgerblue1", "GAM-weights" ="dodgerblue1", 
              "GAM-weights wo M thresh." = "skyblue", "GAM-weights thresh." ="skyblue", 
              "GAM-weights no (T,D)" = "cyan", "GAM-weights no T" = "cadetblue2", 
              "GAM-weights no low C (0.1)" = "darkorchid1", 
              "GAM-weights no low C (0.25)" = "darkorchid3", 
              "GAM-weights no low C (0.5)" = "darkorchid4")
  
  lineVVs <- c("Random Subsample" = "solid", "Unweighted (controls only)" = "solid", "Unweighted" = "solid", 
               "KM-weights wo M" = "dotted", "KM-weights" ="solid", 
               "GAM-weights wo M" = "dotted", "GAM-weights" ="solid", 
               "GAM-weights wo M thresh." = "dotted", "GAM-weights thresh." ="solid", 
               "GAM-weights no (T,D)" = "dashed", "GAM-weights no T" = "dashed", 
               "GAM-weights no low C (0.1)" = "dashed",
               "GAM-weights no low C (0.25)" = "dashed",
               "GAM-weights no low C (0.5)" = "dashed")
  
  
  forggplot_meanSD <- ALL_CRIT %>% group_by(Censoring, Collider, Type) %>% summarise(mean = mean(Value, na.rm=T), sd =sd(Value, na.rm=T))
  
  gplotallcrit <- ggplot(ALL_CRIT, aes(y=Value, fill=Type, linetype = Type)) + 
    facet_grid(Collider~Censoring, scales = "free_y") +
    geom_boxplot(position=position_dodge(1)) + 
    geom_hline(yintercept=0, linetype="longdash", color = "red") +
    scale_fill_manual(values = colVVs) + 
    scale_linetype_manual(values= lineVVs) + 
    theme_light()+ 
    theme(legend.position = "bottom", 
          legend.title    = element_blank(),
          legend.text     = element_text(size=30),
          axis.text.x     = element_blank(),
          axis.ticks.x    = element_blank(),
          axis.text.y     = element_text(size=25), 
          strip.text      = element_text(size=30), 
          axis.title      = element_text(size=25), 
          legend.key.height = unit(2, 'cm'),
          legend.key.width  = unit(4, 'cm')) +
    labs(y = "", x="", fill = "Method", linetype= "Method")
  
  
  
  ggsave(paste0("Figures/DistrX_", suffixFig, ".pdf"), plot=gplotallcrit, height=20, width=30)

  
  
  ##### Performance of Gam-weights without (T,D) vs Diff between NCC all and NCC ctrls
  
  Betas      <- RECAPALL %>% dplyr::select(Censoring, Collider, starts_with(c("beta"))) %>%
    mutate(Impact_ColliderBias = beta_NCC_ctrl - beta_NCC,
           beta_withoutTD = beta_NCC_WeightedGAM_noTD - beta_allcohort) %>%
    select(1:2, Impact_ColliderBias, beta_withoutTD)
  
  
  gplot_withoutTD<- ggplot(Betas, aes(y=beta_withoutTD, x= Impact_ColliderBias)) +
    facet_grid(~Censoring) +
    geom_point(aes(color=Collider)) +
    geom_smooth(method = lm, se = FALSE, linewidth=1.5, aes(color = Collider))  +#aes(group=Collider)
    geom_hline(yintercept=0, linetype="longdash", color = "grey60") +
    theme_light()+
    theme(legend.position = "bottom",
          legend.title    = element_blank(),
          legend.text     = element_text(size=20),
          axis.text.x     = element_text(size=20),
          #axis.ticks.x    = element_blank(),
          axis.text.y     = element_text(size=20),
          strip.text      = element_text(size=30),
          axis.title      = element_text(size=25),
          legend.key.height = unit(1, 'cm'),
          legend.key.width  = unit(2, 'cm')) +
    labs(y = "Bias of analyses based on GAM weights wo (T,D)", x="Diff. between estimates from the unweighted analyses of the controls and of the full NCC", color= "Collider Bias")
  
  ggsave(paste0("Figures/DistrX_", suffixFig, "_FocusGAMWeights_Wo_TD.pdf"), plot = gplot_withoutTD,height=9, width=15)
}