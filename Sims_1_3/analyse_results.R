rm(list = ls())
library(tidyverse)
library(ggplot2)

ROOT <- Sys.getenv("WNCC_ROOT")
stopifnot(nzchar(ROOT))

setwd(file.path(ROOT, "wncc_sims", "Sims_1_3"))

##### Plots comparing GAM and Weights in the NN matching case
oo <- readRDS("./Res/Data_test_Sim_1_3_CohortSize_1e+05_MAF_0.25_R2SNPBMI_0.1_logHRSNP_0.69_logHRBMI_0.69_logHRMatching_0.69_pselCases_1_ageCens_20_50_NNmatching_1_AssoMsnp_1.rds")

quantilesMatch <- quantile(
  oo$cohort_withweights$Matching, probs = c(0.01, 0.25, 0.26, 0.5, 0.51, 0.99)
) 

oo_long <- oo$cohort_withweights %>%
  filter(AgeOut > 48, Ind==0) %>%
  dplyr::select(Matching, gamforall.y , WeightsKM_Matching) %>%
  pivot_longer(cols=-1, names_to = "Weights") %>%
  mutate(
    value = case_when(
    grepl("GAM|gam", Weights)~1/value, T ~ ifelse(value==0, 0, 1/value)
  ))


oo_long <- oo_long %>%
  mutate(
    Matching_range  = case_when(
      Matching < quantilesMatch[1] ~ "< P1",
      Matching > quantilesMatch[6] ~ ">P99",
      Matching < quantilesMatch[3] & Matching > quantilesMatch[2] ~ "P25-P26",
      Matching < quantilesMatch[5] & Matching > quantilesMatch[4] ~ "P50-P51"
    )
  ) %>%
  filter(!is.na(Matching_range)) %>%
  mutate(
    Matching_range = factor(
      as_factor(Matching_range),
      levels = c("< P1", "P25-P26", "P50-P51", ">P99"))
  )

CompWeights_NN <- ggplot(
  oo_long %>%
  filter(!is.na(Matching_range)), aes(x = Matching, y = value, color = Weights)) +
  geom_line() +
  facet_wrap(~Matching_range, scales = "free_x") +
  scale_color_discrete(name= "Inclusion Probability", labels = c("GAM", "KM")) + ylab("")

ggsave(paste0("Figures/Comp_prob_incl_GAM_KM_NNMatching.pdf"), plot=CompWeights_NN, height=6, width=9)


remove_first_bit <- Vectorize(function(char){paste0(unlist(strsplit(char, "_"))[-1], collapse="_")})
sizecoho      <- 1e5
pselCases     <- 1
MAF           <- 0.25
R2SNPBMI      <- 0.1
logHRSNP      <- log(2) #0
logHRBMI      <- log(2)
agemincens    <- 20
agemaxcens    <- 50

##### Simulation 1: NN versus Caliper matching; 
##### Illustration of the limitations of KM-weights 
##### when inclusion probabilities are null for some 
##### participants of the originating cohort


# logHRMatching <- log(2)
# assoMsnp_num  <- 1

RECAP_RES_ASSO <- NULL
RECAPALL0 <- NULL

for (nnmatch in c(0, 1)) {
  if (nnmatch == 1) {
    Match_effect <- list(c(0, 0), c(log(2), 1))  ## for NNmatching, we want to illustrate the fact that ignoring M in the computation of KM-weights, when M ~ Z, makes inference based on KM-weights ok (even if selection effectively depends on S)
  } else {
    Match_effect <- list(c(log(2), 1)) ## for caliper matching we don't need to illustrate this at this stage (this will be taken care of in sim_3)
  }

  for (le in 1:length(Match_effect)){
    logHRMatching <- Match_effect[[le]][1]
    assoMsnp_num  <- Match_effect[[le]][2]
    suffsuffix    <- paste0("_pselCases_", pselCases, "_ageCens_", agemincens, "_", agemaxcens, "_NNmatching_", nnmatch, "_AssoMsnp_", assoMsnp_num) 
    suffix        <- paste0("Sim_1_3_CohortSize_", sizecoho, "_MAF_", MAF, "_R2SNPBMI_", R2SNPBMI, "_logHRSNP_", round(logHRSNP,2),
                         "_logHRBMI_", round(logHRBMI,2), "_logHRMatching_", round(logHRMatching,2), suffsuffix)
    namefile <- paste0("Res/Res_", suffix,".rds")
    if (file.exists(namefile)) {
      recap    <- readRDS(namefile)    
      RECAPALL0 <- rbind(RECAPALL0, cbind(NNmatch = nnmatch, recap$RECAP_RES_ASSO)) # here, assoMsnp_num=1 <=> logHRMatching=log(2); so assoMsnp_num is enough to tell whether M~W or M~Z
    }
  }
}

RECAPALL   <- as_tibble(RECAPALL0) 
RECAPALL   <- RECAPALL %>% mutate(across(c(3:9), ~as_factor(round(.x,2))))
RECAPALL   <- RECAPALL %>% mutate(NNmatch = as_factor(case_when(
  NNmatch == 1 & logHRM==0 & AssoMsnp==0~"NN Matching, M ~ Z",  
  NNmatch == 1 & logHRM==0.69 & AssoMsnp==1~"NN Matching, M ~ W",  
  NNmatch == 0 & logHRM==0.69 & AssoMsnp==1~"Caliper Matching, M ~ W")))
RECAPALL$NNmatch <- factor(RECAPALL$NNmatch, levels=sort(levels(RECAPALL$NNmatch)))
remove_first_bit <- Vectorize(function(char){paste0(unlist(strsplit(str_replace(char, "cond_", ""), "_"))[-1], collapse="_")})

Betas      <- RECAPALL %>% dplyr::select(1:9, starts_with(c("beta"))) %>% mutate(across(starts_with("beta"), ~.x - beta_allcohort)) %>% 
  select(-beta_allcohort) %>% # select_if(!grepl("nomatching", names(.))) %>% #
  gather(Type, beta, starts_with("beta"), factor_key=TRUE) %>% 
  mutate(Typec = as.character(Type), Type = as_factor(remove_first_bit(Typec))) %>% select(-Typec) %>% 
  mutate(noMatchingInWeights = as_factor(case_when( grepl("Weighted", Type) & grepl("nomatching", Type) ~ 1, TRUE~0)))
logHRs_BMI <- RECAPALL %>% dplyr::select(1:9,starts_with(c("alphaBMI"))) %>% mutate(across(starts_with("alphaBMI"), ~.x - alphaBMI_allcohort)) %>% 
  select(-alphaBMI_allcohort) %>% # select_if(!grepl("nomatching", names(.))) %>% 
  gather(Type, logHR, starts_with("alpha"), factor_key=TRUE) %>%
  mutate(Typec = as.character(Type), Type = as_factor(remove_first_bit(Typec))) %>% select(-Typec) %>% 
  mutate(noMatchingInWeights = as_factor(case_when( grepl("Weighted", Type) & grepl("nomatching", Type) ~ 1, TRUE~0)))
# logHRs_SNP <- RECAPALL %>% dplyr::select(1:7,starts_with(c("alphaSNP"))) %>% gather(Type, logHR, starts_with("alpha"), factor_key=TRUE)
KMs        <- RECAPALL %>% dplyr::select(1:9,starts_with(c("KM_cond"))) %>% mutate(across(starts_with("KM_cond"), ~.x - KM_cond_allcohort)) %>% 
  select(-KM_cond_allcohort) %>% # select_if(!grepl("nomatching", names(.))) %>% 
  gather(Type, km, starts_with("KM"), factor_key=TRUE) %>%
  mutate(Typec = as.character(Type), Type = as_factor(remove_first_bit(Typec))) %>% select(-Typec) %>% 
  mutate(noMatchingInWeights = as_factor(case_when( grepl("Weighted", Type) & grepl("nomatching", Type) ~ 1, TRUE~0)))



Betas      <- Betas %>% filter(!grepl("VV", Type)) %>% droplevels() 
logHRs_BMI <- logHRs_BMI %>% filter(!grepl("VV", Type)) %>% mutate(Type = case_when(Type=="NCC"~"Cond. Log. Reg.", T~Type))
KMs        <- KMs %>% filter(!grepl("VV", Type)) 
ALL_CRIT_Sim1   <- bind_rows(
  Betas %>% 
  rename(Value = beta) %>% 
  mutate(Crit = "Linear~Reg.~Coef~ (theta)"),
  logHRs_BMI %>% rename(Value = logHR) %>% mutate(Crit = "log-HR~ (alpha[b])"),
  KMs %>% rename(Value = km) %>% mutate(Crit = "Cond.~Survival"))

ALL_CRIT_Sim1   <- ALL_CRIT_Sim1  %>%
  mutate(Crit = as_factor(Crit), Crit = factor(Crit, levels=levels(Crit)[c(2,3,1)])) %>% 
  mutate(
    Type = fct_recode(
      Type,
      "GAM-weights"                 = "NCC_WeightedGAM_matching",
      "KM-weights"                  = "NCC_WeightedKM_matching",
      "GAM-weights wo M"            = "NCC_WeightedGAM_nomatching",
      "KM-weights wo M"             = "NCC_WeightedKM_nomatching",
      "Cond. Log. Reg."             = "Cond. Log. Reg.",
      "Random Subsample"            = "allcohort_subsample",
      "Unweighted"                  = "NCC",
      "Unweighted (controls only)"  = "NCC_ctrl"
    )
  )

ALL_CRIT_Sim1$Type <- factor(ALL_CRIT_Sim1$Type, levels = c("Random Subsample", "Unweighted (controls only)", "Unweighted",  "Cond. Log. Reg.", "GAM-weights", "KM-weights",  "GAM-weights wo M", "KM-weights wo M"))


######## Simulation 3
######## Some matching factors can be safely ignored when 
######## computing weights (and inclusion probabilities)

nnmatch <- 0

RECAP_RES_ASSO <- NULL
RECAPALL0 <- NULL

for(logHRMatching in c(0, log(2)))  
{
  if(logHRMatching==0) {
    AllAssoMSnp <- c(0)
  }else{AllAssoMSnp <- c(0, 1)}
  for (assoMsnp_num in AllAssoMSnp)
  {
    suffsuffix   <- paste0("_pselCases_", pselCases, "_ageCens_", agemincens, "_", agemaxcens, "_NNmatching_", nnmatch, "_AssoMsnp_", assoMsnp_num) 
    suffix       <- paste0("Sim_1_3_CohortSize_", sizecoho, "_MAF_", MAF, "_R2SNPBMI_", R2SNPBMI, "_logHRSNP_", round(logHRSNP,2),
                         "_logHRBMI_", round(logHRBMI,2), "_logHRMatching_", round(logHRMatching,2), suffsuffix)
    namefile <- paste0("Res/Res_", suffix,".rds")
    if (file.exists(namefile))
    {
      recap    <- readRDS(namefile)    
      RECAPALL0 <- rbind(RECAPALL0, cbind(NNmatch = nnmatch, recap$RECAP_RES_ASSO))
    }
  }
}

RECAPALL   <- as_tibble(RECAPALL0) 
RECAPALL   <- RECAPALL %>% mutate(across(c(1, 3:9), ~as_factor(round(.x,2))))
RECAPALL   <- RECAPALL %>% mutate(AssoMsnp = as_factor(case_when(
                                    logHRM==0 & AssoMsnp==0~"M ~ Z",           
                                    logHRM==0.69 & AssoMsnp==0~"M ~ V",
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



Betas      <- Betas %>% filter(!grepl("VV", Type)) %>% droplevels() 
logHRs_BMI <- logHRs_BMI %>% filter(!grepl("VV", Type)) %>% mutate(Type = case_when(Type=="NCC"~"Cond. Log. Reg.", T~Type))
KMs        <- KMs %>% filter(!grepl("VV", Type)) 
ALL_CRIT_Sim3   <- bind_rows(Betas %>% rename(Value = beta) %>% mutate(Crit = "Linear~Reg.~Coef~ (theta)"), 
                             logHRs_BMI %>% rename(Value = logHR) %>% mutate(Crit = "log-HR~ (alpha[b])"), 
                             KMs %>% rename(Value = km) %>% mutate(Crit = "Cond.~Survival")) 


ALL_CRIT_Sim3   <- ALL_CRIT_Sim3  %>%
  mutate(
    Crit = as_factor(Crit), Crit = factor(Crit, levels=levels(Crit)[c(2,3,1)])
  ) %>%
  mutate(Type = fct_recode(
    Type, "GAM-weights" = "NCC_WeightedGAM_matching",
    "KM-weights" = "NCC_WeightedKM_matching",
    "GAM-weights wo M" = "NCC_WeightedGAM_nomatching",
    "KM-weights wo M" = "NCC_WeightedKM_nomatching",
    "Cond. Log. Reg." = "Cond. Log. Reg.",
    "Random Subsample" = "allcohort_subsample",
    "Unweighted" = "NCC",
    "Unweighted (controls only)" = "NCC_ctrl"
  ))

ALL_CRIT_Sim3$Type <- factor(
  ALL_CRIT_Sim3$Type,
  levels = c(
    "Random Subsample",
    "Unweighted (controls only)",
    "Unweighted",
    "Cond. Log. Reg.",
    "GAM-weights",
    "KM-weights",
    "GAM-weights wo M",
    "KM-weights wo M"
  )
)

colVVs <- c(
  "Random Subsample"           = "gray27",
  "Unweighted (controls only)" = "seagreen4",
  "Unweighted"                 = "palegreen",
  "Cond. Log. Reg."            = "#F8D77A",
  "KM-weights"                 = "#0D2A35",
  "GAM-weights"                = "#5A8FE6",
  "KM-weights wo M"            = "#5D7D87",
  "GAM-weights wo M"           = "#1E5AA8"
)

colVVs_bwcheck <- colorspace::desaturate(colVVs)

lineVVs <- c(
  "Random Subsample"           = "solid",
  "Unweighted (controls only)" = "solid",
  "Unweighted"                 = "solid",
  "Cond. Log. Reg."            = "solid",
  "KM-weights"                 = "dashed",
  "KM-weights wo M"            = "dotted",
  "GAM-weights"                = "solid",
  "GAM-weights wo M"           = "dotted"
)


dat_sim1_red <- ALL_CRIT_Sim1 %>%
  filter(
    NNmatch != "NN Matching, M ~ Z",
    Type %in% c("Cond. Log. Reg.", "KM-weights", "GAM-weights")
  ) %>% droplevels() %>%
  mutate(oo = as_factor(str_replace(as.character(NNmatch), ", M ~ W", "")))

gplot_red_sim1 <- ggplot(dat_sim1_red, aes(y = Value, fill = Type)) +
  facet_wrap(
    Crit ~ oo, scales = "free", ncol = 2, labeller = labeller(Crit = label_parsed)
  ) +
  geom_boxplot(position = position_dodge(1)) +
  geom_hline(yintercept = 0, linetype = "longdash", color = "red") +
  scale_fill_manual(values = colVVs) +
  #  scale_linetype_manual(values= lineVVs) +
  theme_light() +
  theme(legend.position = "bottom",
        legend.title    = element_blank(),
        legend.text     = element_text(size = 30),
        axis.text.x     = element_blank(),
        axis.ticks.x    = element_blank(),
        axis.text.y     = element_text(size = 25),
        strip.text      = element_text(size = 30),
        axis.title      = element_text(size = 25),
        legend.key.height = unit(2, 'cm'),
        legend.key.width  = unit(4, 'cm')) +
  labs(y = "", x = "", fill = "Method", linetype = "Method")

ggsave(paste0("Figures/Sim_1_red_", sizecoho, ".pdf"), plot=gplot_red_sim1, height=20, width=20)
ggsave(
  paste0("Figures/Sim_1_red_", sizecoho, ".eps"),
  plot = gplot_red_sim1,
  width = 20, height = 20,
  dpi = 850
)





gplot_red_sim1_suppMat <- ggplot(
  ALL_CRIT_Sim1 %>% filter(
    NNmatch != "Caliper Matching, M ~ W",
     Type %in% c(
      "Cond. Log. Reg.", "KM-weights", "GAM-weights", "KM-weights wo M")
    ) %>% droplevels(), aes(y=Value, fill=Type, linetype=Type)) + 
  facet_wrap(Crit~NNmatch, scales = "free", ncol=2, labeller = labeller(Crit = label_parsed)) +
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

ggsave(paste0("Figures/Sim_1_red_", sizecoho, "_suppMat.pdf"), plot=gplot_red_sim1_suppMat, height=20, width=20)

gplot_red_sim3 <- ggplot(
  ALL_CRIT_Sim3 %>% 
    filter(Type %in% c(
      "Cond. Log. Reg.","KM-weights", "GAM-weights", 
      "KM-weights wo M", "GAM-weights wo M")
    ) %>% droplevels(), aes(y = Value, fill = Type, linetype = Type)) +
  facet_wrap(
    Crit~AssoMsnp, scales = "free", ncol=3,
    labeller = labeller(Crit = label_parsed)
  ) +
  geom_boxplot(position = position_dodge(1)) +
  geom_hline(yintercept = 0, linetype ="longdash", color = "red") +
  scale_fill_manual(values = colVVs) +
  scale_linetype_manual(values = lineVVs) +
  theme_light() +
  theme(legend.position = "bottom",
        legend.title    = element_blank(),
        legend.text     = element_text(size = 30),
        axis.text.x     = element_blank(),
        axis.ticks.x    = element_blank(),
        axis.text.y     = element_text(size = 25),
        strip.text      = element_text(size = 30),
        axis.title      = element_text(size = 25),
        legend.key.height = unit(2, 'cm'),
        legend.key.width  = unit(4, 'cm')) +
  labs(y = "", x ="", fill = "Method", linetype = "Method")

ggsave(paste0("Figures/Sim_3_red_", sizecoho, ".pdf"), plot=gplot_red_sim3, height=20, width=30)

ggsave(
  paste0("Figures/Sim_3_red_", sizecoho, ".eps"),
  plot = gplot_red_sim3,
  width = 20, height = 20,
  dpi = 850
)
