#### Estimates of the survival and conditional survival functions
#### 

km_estim <- function(cohort, suffix_file_plotkm=NULL)
{
  cohort_all <- cohort %>% select(-any_of("WeightsGAM_Matching_noTD"))
  cohort_all <- cohort_all %>% pivot_longer(cols = which(grepl("Weights", colnames(cohort_all))), names_to = "TypeAnalysis", values_to = "Weights")
  fitKM      <- survfit(Surv(AgeOut, Ind) ~ TypeAnalysis, weights = Weights, data=cohort_all) 
  kmplot     <- ggsurvplot(fitKM , data=cohort_all , ylim =c(0.2, 1), censor=FALSE, conf.int = F)$plot + guides(colour = guide_legend(nrow = 2))#xlim = c(ageforcheck, 80),
  if(!is.null(suffix_file_plotkm))
  { 
    ggsave(paste0("Figures/KM_",  suffix_file_plotkm, ".png"),
           plot=print(kmplot$plot), scale = 1, width = 10, height= 8, units = c("in"), dpi = 300)
  }
  return(list(fitKM = fitKM, kmplot=kmplot))
}

km_cond_estim <- function(cohort, suffix_file_plotkm=NULL)
{
  cohort_all <- cohort %>% select(-any_of("WeightsGAM_Matching_noTD"))
  cohort_all <- cohort_all %>% pivot_longer(cols = which(grepl("Weights", colnames(cohort_all))), names_to = "TypeAnalysis", values_to = "Weights")
  fitKM      <- survfit(Surv(AgeOut, Ind) ~ TypeAnalysis, weights = Weights, data=subset(cohort_all, SNP == min(cohort_all$SNP))) 
  kmplot     <- ggsurvplot(fitKM , data=cohort_all , ylim =c(0.2, 1), censor=FALSE, conf.int = F)$plot + guides(colour = guide_legend(nrow = 2))#xlim = c(ageforcheck, 80),
  if(!is.null(suffix_file_plotkm))
  { 
    ggsave(paste0("Figures/KM_",  suffix_file_plotkm, ".png"),
           plot=print(kmplot$plot), scale = 1, width = 10, height= 8, units = c("in"), dpi = 300)
  }
  return(list(fitKM = fitKM, kmplot=kmplot))
}