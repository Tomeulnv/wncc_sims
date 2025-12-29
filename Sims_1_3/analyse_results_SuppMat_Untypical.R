rm(list=ls())

 
ROOT <- Sys.getenv("Epic_simulations")
stopifnot(nzchar(ROOT))

setwd(ROOT)


library(tidyverse)
library(ggplot2)



sizecoho      <- 1e5
MAF           <- 0.25
R2SNPBMI      <- 0.1
logHRSNP      <- log(2) #0
logHRBMI      <- log(2)
agemincens    <- 20
agemaxcens    <- 50


######## Illustration of the performance of GAM-weights for untypical NCC studies

nnmatching <- nnmatch <- 0

RECAP_RES_ASSO <- NULL
RECAPALL0 <- NULL
for(logHRMatching in c(0, log(2)))  
{
  if(logHRMatching==0) {
    AllAssoMSnp <- c(0)
  }else{AllAssoMSnp <- c(1)}
  for (assoMsnp_num in AllAssoMSnp)
  {
    for (pselCases in c(1, 0.5))
    {
      CohortSize   <- round(sizecoho/pselCases)
      suffix <- paste0("Sim_Supp_Untypical_CohortSize_", CohortSize, "_MAF_", MAF, "_R2SNPBMI_", R2SNPBMI, "_logHRSNP_", round(logHRSNP,2),
                       "_logHRBMI_", round(logHRBMI,2), "_logHRMatching_", round(logHRMatching,2), 
                       "_pselCases_", pselCases, "_ageCens_", agemincens, "_", agemaxcens, 
                       "_NNmatching_", nnmatching, "_AssoMsnp_", assoMsnp_num)
      namefile <- paste0("Res/Res_", suffix,".rds")
      if (file.exists(namefile))
      {
        recap    <- readRDS(namefile)    
        RECAPALL0 <- rbind(RECAPALL0, cbind(NNmatch = nnmatch, recap$RECAP_RES_ASSO))
      }
    }
  }
}

RECAPALL   <- as_tibble(RECAPALL0) 
RECAPALL   <- RECAPALL %>% mutate(across(c(1, 3:9), ~as_factor(round(.x,2))))
RECAPALL   <- RECAPALL %>% mutate(AssoMsnp = as_factor(case_when(
                                    logHRM==0 & AssoMsnp==0~"M ~ Z",           
                                    logHRM==0.69 & AssoMsnp==1~"M ~ W")))

# RECAPALL$AssoMsnp  <- factor(RECAPALL$AssoMsnp, levels=levels(RECAPALL$AssoMsnp))

remove_first_bit <- Vectorize(function(char){paste0(unlist(strsplit(str_replace(char, "cond_", ""), "_"))[-1], collapse="_")})

Betas      <- RECAPALL %>% dplyr::select(1:9, starts_with(c("beta"))) %>% mutate(across(starts_with("beta"), ~.x - beta_allcohort)) %>% 
  select(-beta_allcohort) %>% 
  gather(Type, beta, starts_with("beta"), factor_key=TRUE) %>% 
  mutate(Typec = as.character(Type), Type = as_factor(remove_first_bit(Typec))) %>% select(-Typec) %>% 
  mutate(noMatchingInWeights = as_factor(case_when( grepl("Weighted", Type) & grepl("nomatching", Type) ~ 1, TRUE~0)))
logHRs_BMI <- RECAPALL %>% dplyr::select(1:9,starts_with(c("alphaBMI"))) %>% mutate(across(starts_with("alphaBMI"), ~.x - alphaBMI_allcohort)) %>% 
  select(-alphaBMI_allcohort) %>% 
  gather(Type, logHR, starts_with("alpha"), factor_key=TRUE) %>%
  mutate(Typec = as.character(Type), Type = as_factor(remove_first_bit(Typec))) %>% select(-Typec) %>% 
  mutate(noMatchingInWeights = as_factor(case_when( grepl("Weighted", Type) & grepl("nomatching", Type) ~ 1, TRUE~0)))
# logHRs_SNP <- RECAPALL %>% dplyr::select(1:7,starts_with(c("alphaSNP"))) %>% gather(Type, logHR, starts_with("alpha"), factor_key=TRUE)
KMs        <- RECAPALL %>% dplyr::select(1:9,starts_with(c("KM_cond"))) %>% mutate(across(starts_with("KM_cond"), ~.x - KM_cond_allcohort)) %>% 
  select(-KM_cond_allcohort) %>% 
  gather(Type, km, starts_with("KM"), factor_key=TRUE) %>%
  mutate(Typec = as.character(Type), Type = as_factor(remove_first_bit(Typec))) %>% select(-Typec) %>% 
  mutate(noMatchingInWeights = as_factor(case_when( grepl("Weighted", Type) & grepl("nomatching", Type) ~ 1, TRUE~0)))



Betas            <- Betas %>% filter(!grepl("VV", Type)) %>% droplevels() 
logHRs_BMI       <- logHRs_BMI %>% filter(!grepl("VV", Type)) %>% mutate(Type = case_when(Type=="NCC"~"Cond. Log. Reg.", T~Type))
KMs              <- KMs %>% filter(!grepl("VV", Type)) 
ALL_CRIT_Untyp   <- bind_rows(Betas %>% rename(Value = beta) %>% mutate(Crit = "Linear~Reg.~Coef~ (theta)"), 
                             logHRs_BMI %>% rename(Value = logHR) %>% mutate(Crit = "log-HR~ (alpha[b])"), 
                             KMs %>% rename(Value = km) %>% mutate(Crit = "Cond.~Survival")) 


ALL_CRIT_Untyp   <- ALL_CRIT_Untyp  %>% 
                      mutate(Crit = as_factor(Crit), Crit = factor(Crit, levels=levels(Crit)[c(2,3,1)])) %>% 
                      mutate(Type = fct_recode(Type, "GAM-weights" = "NCC_WeightedGAM_matching", 
                                               "KM-weights" = "NCC_WeightedKM_matching", 
                                               "GAM-weights wo M" = "NCC_WeightedGAM_nomatching", 
                                               "Cond. Log. Reg." = "Cond. Log. Reg.", 
                                               "Random Subsample" = "allcohort_subsample", 
                                               "Unweighted" = "NCC", 
                                               "Unweighted (controls only)" = "NCC_ctrl"))

ALL_CRIT_Untyp$Type <- factor(ALL_CRIT_Untyp$Type, levels = c("Random Subsample", "Unweighted (controls only)", "Unweighted",  
                                                            "Cond. Log. Reg.", "GAM-weights", "KM-weights",  "GAM-weights wo M", "KM-weights wo M"))

ALL_CRIT_Untyp <- ALL_CRIT_Untyp %>% mutate(pselCases = paste0("pi[1] == ", pselCases))

colVVs <- c("Random Subsample" = "gray27", "Unweighted (controls only)" = "seagreen4", "Unweighted" = "palegreen",
            "Cond. Log. Reg." = "purple", 
            "KM-weights" ="sienna1", 
            "GAM-weights" ="dodgerblue1", 
            "KM-weights wo M" ="sienna1", 
            "GAM-weights wo M" ="dodgerblue1")

lineVVs <- c("Random Subsample" = "solid", "Unweighted (controls only)" = "solid", "Unweighted" = "solid",
             "Cond. Log. Reg." = "solid", 
             "KM-weights wo M" = "dotted", "KM-weights" ="solid", 
             "GAM-weights wo M" = "dotted", "GAM-weights" ="solid")




for (leAsso in levels(ALL_CRIT_Untyp$AssoMsnp))
{
  
  gplot_red_untyp_1 <- ggplot(ALL_CRIT_Untyp %>% filter(Type %in% c("Cond. Log. Reg.", 
                                                              "KM-weights", "GAM-weights", 
                                                              "GAM-weights wo M"), AssoMsnp==leAsso) %>% droplevels(), aes(y=Value, fill=Type, linetype=Type)) + 
  facet_grid(Crit~pselCases, scales = "free_y", labeller = labeller(pselCases = label_parsed, Crit = label_parsed))  + #ncol=2, 
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

  ggsave(paste0("Figures/Sim_Untyp_red_", str_replace(leAsso, " ~ ", "_"), ".pdf"), plot=gplot_red_untyp_1, height=20, width=30)
}
