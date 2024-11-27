library(connector)
library(dplyr)
library(survival)
library(survminer)

Data <- readRDS( system.file("Data","WholeData.RDs", package = "ShinyMCL") )

# names(WholeData) = c("MCL0208", "YoungerMCL" )
# saveRDS(WholeData, "inst/Data/WholeData.RDs")

SU_Younger = list()

SUgeneration = function(tissueMCL,WholeData=Data){
  
  CONNECTORList.FCM = WholeData$MCL0208[[tissueMCL]]$CONNECTORList.FCM
  classificationFCM = WholeData$YoungMCL[[tissueMCL]]$Dataset
  
  
  ClassMatrix_entropyInit= ClassificationNewCurves(classificationFCM,
                                                   CONNECTORList.FCM, Cores = 2)
  ClassMatrix_entropy =  ClassMatrix_entropyInit$ClassMatrix_entropy %>%
    filter(Cluster != "Unclassified")
  Info_tmp= merge(classificationFCM %>% select(ID,INIERG,rnd1, ttpevent, ttp) %>%
                    distinct() %>% mutate(ttpDiff = ttp*365.25 ) %>% as.data.frame(), 
                  as.data.frame(ClassMatrix_entropy), by = "ID")
  
  SmallNumberIDs = Info_tmp %>% group_by(Cluster) %>% dplyr::summarise(N = n()) %>% filter(N <3)
  if(dim(SmallNumberIDs)[1] > 0 )
  {
    Info_tmp = Info_tmp %>% filter(! Cluster %in% SmallNumberIDs$Cluster )
  }
  
  ############################ SURVIVAL ANALYSIS
  fit <- eval(parse(text = paste0("survfit(Surv(Info_tmp$ttpDiff,Info_tmp$ttpevent) ~ Cluster, data = Info_tmp)")))
  surv_median <- as.vector(summary(fit)$table[, "median"])
  df <- data.frame(x1 = surv_median, x2 = surv_median,
                   y1 = rep(0, length(surv_median)), y2 = rep(0.5, length(surv_median)))
  list2env(list(Info_tmp = Info_tmp), envir = .GlobalEnv)
  ggsurv=ggsurvplot(fit = fit, data = Info_tmp,
                    xlab = "Year",  ylab = "TTP",
                    size = 1, pval = TRUE, risk.table = TRUE,
                    risk.table.col="strata",ggtheme = theme_bw() )
  
  list(fit = fit,
       Info_tmp = Info_tmp,
       ggsurvData = ggsurv$data.survplot,
       Classification = ClassMatrix_entropyInit["ClassMatrix_entropy"])  
  
}


SU_Younger$BM = SUgeneration("BM")
SU_Younger$PB = SUgeneration("PB")


saveRDS(SU_Younger, "inst/Data/YoungerMCLsurvival.RDs")
