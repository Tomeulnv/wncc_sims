# --- CLEAR ENVIRONMENT & LOAD LIBRARIES ---
rm(list = ls())

library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggh4x)

# --- SET WORKING DIRECTORY ---
setwd("/home/PERSONALE/tomeu.lopeznieto/wNCC/Res")

# --- HELPER FUNCTION ---
remove_first_bit <- Vectorize(function(char) {
  paste0(
    unlist(strsplit(str_replace(char, "cond_", ""), "_"))[-1],
    collapse = "_"
  )
})

# --- PARAMETER SPECIFICATIONS ---
NSIM             <- 100
pselCases        <- 1
MAF              <- 0.25
R2SNPBMI         <- 0.1
logHRSNP         <- log(2)
logHRBMI         <- log(2)
agemincens       <- 20
agemaxcens       <- 50
Asso_M_SNP       <- 0.2
R2BMI            <- 0.1
ALLlogHRM        <- log(c(1, 2))
ALLlogHRM1M2     <- log(c(1, 2))
ALLlogHRSNP_M1M2 <- log(c(1, 2))
matchMethodLabel <- "Caliper"

# --- IMPORT SIMULATION RESULTS ---
RECAP_RES_ASSO <- NULL
RECAPALL0      <- NULL

for (logHRM in ALLlogHRM)
for (logHRM1M2 in ALLlogHRM1M2)
for (logHRSNP_M1M2 in ALLlogHRSNP_M1M2) {

  suffix <- paste0(
    "InterMatching_onSurvival_MAF_", MAF,
    "_Asso_M_SNP_", Asso_M_SNP,
    "_R2SNPBMI_", R2SNPBMI,
    "_logHRSNP_", round(logHRSNP, 2),
    "_logHRBMI_", round(logHRBMI, 2),
    "_logHRM_", round(logHRM, 2),
    "_logHRM1M2_", round(logHRM1M2, 2),
    "_logHRSNP_M1M2_", round(logHRSNP_M1M2, 2),
    "_pselCases_", pselCases,
    "_ageCens_", agemincens, "_", agemaxcens,
    "_", matchMethodLabel
  )

  namefile <- paste0("Res_", suffix, ".rds")

  if (file.exists(namefile)) {
    recap     <- readRDS(namefile)
    RECAPALL0 <- rbind(RECAPALL0, recap$RECAP_RES_ASSO)
  }
}

# --- CLEAN AND FORMAT DATA ---
RECAPALL <- as_tibble(RECAPALL0) %>%
  mutate(across(2:9, ~ as_factor(round(.x, 2))))

levels(RECAPALL$logHRM)        <- paste0("alpha[M] == ",     levels(RECAPALL$logHRM))
levels(RECAPALL$logHRM1M2)     <- paste0("alpha[M1M2] == ",  levels(RECAPALL$logHRM1M2))
levels(RECAPALL$logHRSNP_M1M2) <- paste0("alpha[MSNP] == ",  levels(RECAPALL$logHRSNP_M1M2))

###################################################################
## OBSERVED DISTRIBUTION PLOTS (DIFF TO ALL COHORT)
###################################################################
# --- BETAS ---
Betas <- RECAPALL %>%
  dplyr::select(1:9, starts_with("beta"), logHRSNP_M1M2) %>%
  mutate(across(starts_with("beta"), ~ .x - beta_allcohort)) %>%
  select(-any_of(c(
    "beta_allcohort",
    "beta_NCC_WeightedGAMth_nomatching",
    "beta_NCC_WeightedGAMth_matching",
    "beta_NCC_WeightedGAM_matching_noTD"
  ))) %>%
  gather(Type, beta, starts_with("beta"), factor_key = TRUE) %>%
  mutate(
    Typec = as.character(Type),
    Type = as_factor(remove_first_bit(Typec))
  ) %>%
  select(-Typec) %>%
  mutate(noMatchingInWeights = as_factor(case_when(
    grepl("Weighted", Type) & grepl("nomatching", Type) ~ 1,
    TRUE ~ 0
  )))

# --- LOGHRs BMI ---
logHRs_BMI <- RECAPALL %>%
  dplyr::select(1:9, starts_with("alphaBMI"), logHRSNP_M1M2) %>%
  mutate(across(starts_with("alphaBMI"), ~ .x - alphaBMI_allcohort)) %>%
  select(-any_of(c(
    "alphaBMI_allcohort",
    "alphaBMI_NCC_WeightedGAMth_nomatching",
    "alphaBMI_NCC_WeightedGAMth_matching"
  ))) %>%
  gather(Type, logHR, starts_with("alpha"), factor_key = TRUE) %>%
  mutate(
    Typec = as.character(Type),
    Type = as_factor(remove_first_bit(Typec))
  ) %>%
  select(-Typec) %>%
  mutate(noMatchingInWeights = as_factor(case_when(
    grepl("Weighted", Type) & grepl("nomatching", Type) ~ 1,
    TRUE ~ 0
  )))

# --- KM Estimates ---
KMs <- RECAPALL %>%
  dplyr::select(1:9, starts_with("KM_cond"), logHRSNP_M1M2) %>%
  mutate(across(starts_with("KM_cond"), ~ .x - KM_cond_allcohort)) %>%
  select(-any_of(c(
    "KM_cond_allcohort",
    "KM_cond_NCC_WeightedGAMth_nomatching",
    "KM_cond_NCC_WeightedGAMth_matching"
  ))) %>%
  gather(Type, km, starts_with("KM"), factor_key = TRUE) %>%
  mutate(
    Typec = as.character(Type),
    Type = as_factor(remove_first_bit(Typec))
  ) %>%
  select(-Typec) %>%
  mutate(noMatchingInWeights = as_factor(case_when(
    grepl("Weighted", Type) & grepl("nomatching", Type) ~ 1,
    TRUE ~ 0
  )))

# --- CLEAN & RELABEL ESTIMANDS ---
Betas      <- Betas %>% filter(!grepl("VV", Type)) %>% droplevels()
logHRs_BMI <- logHRs_BMI %>%
  filter(!grepl("VV", Type)) %>%
  mutate(Type = case_when(Type == "NCC" ~ "Cond. Log. Reg.", TRUE ~ Type))
KMs        <- KMs %>% filter(!grepl("VV", Type))

# --- JOIN ESTIMANDS ---
ALL_CRIT <- bind_rows(
  Betas      %>% rename(Value = beta)  %>% mutate(Crit = "Linear~Reg.~Coef~ (theta)"),
  logHRs_BMI %>% rename(Value = logHR) %>% mutate(Crit = "log-HR~ (alpha[b])"),
  KMs        %>% rename(Value = km)    %>% mutate(Crit = "Cond.~Survival")
) %>%
  mutate(Crit = as_factor(Crit)) %>%
  mutate(Type = fct_recode(Type,
    "GAM-weights"             = "NCC_WeightedGAM_matching",
    "GAM-weights wo M"        = "NCC_WeightedGAM_nomatching",
    "KM-weights"              = "NCC_WeightedKM_matching",
    "Cond. Log. Reg."         = "Cond. Log. Reg.",
    "GAM-weights with InterM" = "NCC_WeightedGAM_matchingInter"
  ))

ALL_CRIT$Type <- factor(ALL_CRIT$Type, levels = c(
  "Cond. Log. Reg.",
  "GAM-weights",
  "KM-weights",
  "GAM-weights wo M",
  "GAM-weights with InterM"
))

# --- PLOT COLORS & LINE STYLES ---
colVVs <- c(
  "Cond. Log. Reg."              = "#F8D77A",  # L≈82
  "KM-weights"                   = "#0D2A35",  # L≈16
  "GAM-weights"                  = "#5A8FE6",  # L≈63 (mid blue)
  "GAM-weights wo M"             = "#D8B4F3",  # L≈78 (if used elsewhere)
  "GAM-weights with InterM"      = "#9B2F2A",  # L≈35 (much darker, warm)
  "NCC_WeightedGAMth_nomatching" = "#8A8A8A",  # L≈58
  "NCC_WeightedGAMth_matching"   = "#3A3A3A"   # L≈23
)

lineVVs <- c(
  "Cond. Log. Reg."              = "solid",
  "KM-weights"                   = "dashed",
  "GAM-weights"                  = "solid",
  "GAM-weights wo M"             = "dotted",
  "GAM-weights with InterM"      = "dotdash",
  "NCC_WeightedGAMth_nomatching" = "dotted",
  "NCC_WeightedGAMth_matching"   = "solid"
)

# Check if Printed in Black and White
colVVs_bwcheck <- colorspace::desaturate(colVVs)


######## PLOT OF ALL ESTIMANDS ########
make_alpha_plot <- function(alpha_label) {
  ALL_CRIT %>%
    filter(Type %in% c("Cond. Log. Reg.",
                       "KM-weights",
                       "GAM-weights",
                       "GAM-weights with InterM")) %>%
    filter(logHRM == alpha_label) %>%
    droplevels() %>%
    mutate(
      Type = factor(Type, levels = c("Cond. Log. Reg.",
                                     "KM-weights",
                                     "GAM-weights",
                                     "GAM-weights with InterM")),
      Crit = factor(Crit, levels = c(
        "log-HR~ (alpha[b])",
        "Cond.~Survival",
        "Linear~Reg.~Coef~ (theta)"
      ))
    ) %>%
    ggplot(aes(y = Value, fill = Type, linetype = Type)) +
    facet_nested(
      Crit ~ logHRSNP_M1M2 + logHRM1M2,
      labeller = labeller(
        Crit = label_parsed,
        logHRSNP_M1M2 = label_parsed,
        logHRM1M2 = label_parsed
      ),
      scales = "free_y"
    ) +
    geom_boxplot(position = position_dodge(1)) +
    geom_hline(yintercept = 0, linetype = "longdash", color = "red") +
    scale_fill_manual(values = colVVs) +
    scale_linetype_manual(values = lineVVs) +
    theme_light() +
    theme(
      legend.position     = "bottom",
      legend.title        = element_blank(),
      legend.text         = element_text(size = 24),
      axis.text.x         = element_blank(),
      axis.ticks.x        = element_blank(),
      axis.text.y         = element_text(size = 25),
      strip.text          = element_text(size = 30),
      axis.title          = element_text(size = 25),
      legend.key.height   = unit(3, "cm"),
      legend.key.width    = unit(5, "cm")
    ) +
    labs(y = "", x = "", fill = "Method", linetype = "Method")
}

gplotallcrit_0    <- make_alpha_plot("alpha[M] == 0")
gplotallcrit_069  <- make_alpha_plot("alpha[M] == 0.69")

ggsave("Estimands_diffs_Inters_logHRM_0.pdf",   gplotallcrit_0,   height = 20, width = 30)
ggsave("Estimands_diffs_Inters_logHRM_0.69.pdf", gplotallcrit_069, height = 20, width = 30)


# ###################################################################
# ## 2) MEAN; SD & CI TABLE
# ###################################################################
# Betas <- RECAPALL %>%
#   select(1:9, starts_with("beta"), logHRSNP_M1M2) %>%
#   gather(Type, beta, starts_with("beta"), factor_key = TRUE) %>%
#   mutate(
#     Typec = as.character(Type),
#     Type = as_factor(remove_first_bit(Typec)),
#     Type = fct_recode(Type, "Full Cohort" = "allcohort"),
#     noMatchingInWeights = as_factor(case_when(
#       grepl("Weighted", Type) & grepl("nomatching", Type) ~ 1,
#       TRUE ~ 0
#     ))
#   ) %>%
#   select(-Typec)

# logHRs_BMI <- RECAPALL %>%
#   select(1:9, starts_with("alphaBMI"), logHRSNP_M1M2) %>%
#   gather(Type, logHR, starts_with("alphaBMI"), factor_key = TRUE) %>%
#   mutate(
#     Typec = as.character(Type),
#     Type = as_factor(remove_first_bit(Typec)),
#     Type = fct_recode(Type, "Full Cohort" = "allcohort"),
#     noMatchingInWeights = as_factor(case_when(
#       grepl("Weighted", Type) & grepl("nomatching", Type) ~ 1,
#       TRUE ~ 0
#     ))
#   ) %>%
#   select(-Typec)

# KMs <- RECAPALL %>%
#   select(1:9, starts_with("KM_cond"), logHRSNP_M1M2) %>%
#   gather(Type, km, starts_with("KM_cond"), factor_key = TRUE) %>%
#   mutate(
#     Typec = as.character(Type),
#     Type = as_factor(remove_first_bit(Typec)),
#     Type = fct_recode(Type, "Full Cohort" = "allcohort"),
#     noMatchingInWeights = as_factor(case_when(
#       grepl("Weighted", Type) & grepl("nomatching", Type) ~ 1,
#       TRUE ~ 0
#     ))
#   ) %>%
#   select(-Typec)

# # --- CLEAN & RELABEL ESTIMANDS ---
# Betas <- Betas %>%
#   filter(!grepl("VV", Type)) %>%
#   filter(Type %in% c(
#     "Full Cohort",
#     "NCC_WeightedKM_matching",
#     "NCC_WeightedGAM_matching",
#     "NCC_WeightedGAM_matchingInter"
#   )) %>%
#   droplevels()

# logHRs_BMI <- logHRs_BMI %>%
#   filter(!grepl("VV", Type)) %>%
#   filter(Type %in% c(
#     "Full Cohort",
#     "NCC",
#     "NCC_WeightedKM_matching",
#     "NCC_WeightedGAM_matching",
#     "NCC_WeightedGAM_matchingInter"
#   )) %>%
#   mutate(Type = case_when(Type == "NCC" ~ "Cond. Log. Reg.", TRUE ~ Type))

# KMs <- KMs %>%
#   filter(!grepl("VV", Type)) %>%
#   filter(Type %in% c(
#     "Full Cohort",
#     "NCC_WeightedKM_matching",
#     "NCC_WeightedGAM_matching",
#     "NCC_WeightedGAM_matchingInter"
#   )) %>%
#   droplevels()

# # Auxilliary Function
# compute_summary <- function(df, value_col, crit_label) {
#   df %>%
#     group_by(Type, logHRM, logHRM1M2, logHRSNP_M1M2) %>%
#     summarise(
#       mean = mean(.data[[value_col]], na.rm = TRUE),
#       sd = sd(.data[[value_col]], na.rm = TRUE),
#       lower_ci = mean - 1.96 * (sd / sqrt(max(RECAPALL$nsim))),
#       upper_ci = mean + 1.96 * (sd / sqrt(max(RECAPALL$nsim))),
#       .groups = "drop"
#     ) %>%
#     mutate(Crit = crit_label)
# }

# # Compute summary tables
# summary_beta <- compute_summary(Betas, "beta", "Correlation")
# summary_logHR <- compute_summary(logHRs_BMI, "logHR", "log-HR")
# summary_KM <- compute_summary(KMs, "km", "Cond. Survival")

# # Combine summaries with clean labels
# summary_all <- bind_rows(summary_beta, summary_logHR, summary_KM) %>%
#   mutate(Type = fct_recode(Type,
#     "GAM-weights"             = "NCC_WeightedGAM_matching",
#     "KM-weights"              = "NCC_WeightedKM_matching",
#     "Cond. Log. Reg."         = "Cond. Log. Reg.",
#     "GAM-weights with InterM" = "NCC_WeightedGAM_matchingInter"
#   )) %>%
#   select(Crit, Type, logHRM, logHRM1M2, logHRSNP_M1M2, mean, sd, lower_ci, upper_ci) %>%
#   arrange(Crit, Type, logHRM, logHRM1M2, logHRSNP_M1M2)

# summary_table <- summary_all %>%
#   mutate(
#     mean = round(mean, 4),
#     sd = round(sd, 4),
#     lower_ci = round(lower_ci, 4),
#     upper_ci = round(upper_ci, 4)
#   ) %>%
#   select(
#     Crit,
#     logHRM, logHRM1M2, logHRSNP_M1M2,
#     Type,
#     mean, sd, lower_ci, upper_ci
#   ) %>%
#   arrange(Crit, logHRM, logHRM1M2, logHRSNP_M1M2, Type)

# library(writexl)
# write_xlsx(summary_table, path = "../summary_table.xlsx")

# ###################################################################
# ## 2) CONFIDENCE INTERVALS PLOTS
# ###################################################################

# KMs        <- RECAPALL %>% dplyr::select(1:8,starts_with(c("KM_cond"))) %>%# select(-starts_with(c("KM_cond"))) %>%
#                           #mutate(across(starts_with("KM_cond"), ~.x - KM_cond_allcohort)) %>%
#                           select(-any_of(c(#"KM_cond_allcohort",
#                           "KM_cond_NCC_WeightedGAMth_nomatching",
#                           "KM_cond_NCC_WeightedGAMth_matching"))) %>%
#                           # mutate(across(starts_with("KM"), ~.x - KM_allcohort)) %>% 
#                           # select(-any_of(c("KM_allcohort", "KM_NCC_WeightedGAMth_nomatching", 
#                           #                  "KM_NCC_WeightedGAMth_matching"))) %>% 
#                           gather(Type, km, starts_with("KM"), factor_key=TRUE) %>%
#                           mutate(Typec = as.character(Type), Type = as_factor(remove_first_bit(Typec))) %>% 
#                           select(-Typec) %>% 
#                           mutate(noMatchingInWeights = 
#                                    as_factor(case_when( grepl("Weighted", Type) & grepl("nomatching", Type) ~ 1, TRUE~0))) %>%
#                           mutate(Type = fct_recode(Type,
#                                       "Full Cohort" = "allcohort",
#                                       "GAM-weights" = "NCC_WeightedGAM_matching", 
#                                       "GAM-weights wo M" = "NCC_WeightedGAM_nomatching", 
#                                       "KM-weights" = "NCC_WeightedKM_matching",
#                                       "GAM-weights with InterM" = "NCC_WeightedGAM_matchingInter"))



# KMs$Type <- factor(KMs$Type, levels = c("Full Cohort", "Cond. Log. Reg.",  "KM-weights",
#                                                   "GAM-weights", "GAM-weights wo M", "GAM-weights with InterM"))

# colVVs <- c(
#   "Full Cohort" = "palegreen",
#   "KM-weights" = "sienna1",
#   "GAM-weights" = "dodgerblue1",
#   "GAM-weights wo M" = "dodgerblue1",  
#   "GAM-weights with InterM" = "dodgerblue1"
# )

# lineVVs <- c(
#   "Full Cohort" = "solid",
#   "KM-weights" = "solid",  
#   "GAM-weights" = "solid", 
#   "GAM-weights wo M" = "dotted",
#   "GAM-weights with InterM" = "dashed"
# )


# forggplot_meanSD <- KMs %>% group_by(Type, logHRM, logHRM1M2) %>% 
#   summarise(mean = mean(km, na.rm=T),
#             sd = sd(km, na.rm=T), 
#             lower_ci = mean - 1.96*(sd/sqrt(max(RECAPALL$nsim))),
#             upper_ci = mean + 1.96*(sd/sqrt(max(RECAPALL$nsim))))

# #### CI PLOT
# gplotallcritIC <- ggplot(forggplot_meanSD, aes(x = Type,  linetype = Type)) + 
#   geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci, y = mean, color = Type), stat = "identity", linewidth = 1.2) +
#   geom_point(aes(y = mean, color = Type), size = 3) +
#   facet_grid(logHRM ~ logHRM1M2, 
#              labeller = labeller(logHRM = label_parsed, logHRM1M2 = label_parsed), 
#              scales = "free_y") +
#   scale_color_manual(values = colVVs) +
#   scale_linetype_manual(values = lineVVs) +
#   theme_light() +
#   theme(legend.position = "bottom", 
#         legend.title    = element_blank(),
#         legend.text     = element_text(size=15),
#         axis.text.x     = element_blank(),
#         axis.ticks.x    = element_blank(),
#         axis.text.y     = element_text(size=15), 
#         strip.text      = element_text(size=25), 
#         axis.title      = element_text(size=25), 
#         legend.key.height = unit(2, 'cm'),
#         legend.key.width  = unit(4, 'cm')) +
#   labs(y = "", x = "")

# ggsave(paste0("KM_IC_", suffix, ".pdf"), plot=gplotallcritIC,height=20, width=20)

# #### Boxplot of CIs
# CI_boxplot  <- ggplot(forggplot_meanSD, aes(x = Type, color = Type, linetype = Type)) + 
#   geom_boxplot(aes(lower = lower_ci, upper = upper_ci, middle = mean,
#   ymin = mean - 2.57*(sd/sqrt(max(RECAPALL$nsim))),
#   ymax = mean + 2.57*(sd/sqrt(max(RECAPALL$nsim))),
#   fill=Type), stat="identity", color="black") +
#   facet_grid(logHRM ~ logHRM1M2, 
#              labeller = labeller(logHRM = label_parsed, logHRM1M2 = label_parsed), 
#              scales = "free_y") +  
#   scale_color_manual(values = colVVs) +
#   scale_fill_manual(values = colVVs) +
#   scale_linetype_manual(values = lineVVs) +
#   theme_light() + 
#   theme(legend.position = "bottom", 
#         legend.title    = element_blank(),
#         legend.text     = element_text(size=15),
#         axis.text.x     = element_blank(),
#         axis.ticks.x    = element_blank(),
#         axis.text.y     = element_text(size=15), 
#         strip.text      = element_text(size=25), 
#         axis.title      = element_text(size=25), 
#         legend.key.height = unit(2, 'cm'),
#         legend.key.width  = unit(4, 'cm')) +
#   labs(y = "", x = "") 

# ggsave(paste0("KM_IC_Box", suffix, ".pdf"), plot=CI_boxplot,height=20, width=20)

# ###################################################################
# # OUTPUT TABLE OF RESULTS WITH CIs
# ###################################################################

# library(knitr)
# library(kableExtra)

# forggplot_meanSD_sorted <- forggplot_meanSD %>%
#   arrange(logHRM, logHRM1M2)

# html_table  <- kable(forggplot_meanSD_sorted, format = "html", caption = "Sorted Results by logHRM and logHRM1M2")%>%
#   kable_styling(bootstrap_options = c("striped", "hover")) %>%
#    save_kable(file = "sorted_table.html")

#kable(forggplot_meanSD_sorted, format = "latex", caption = "Sorted Results by logHRM and logHRM1M2")