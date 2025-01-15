load("Data/Connector_dataset_ArmUniti.RData")
source("test_funct_MCL0208.R")

qual_vars<-c("PFS","OS","TTP","TP53_loss_or_mut","TP53","TP53_loss_array","TP53_disruption","Ki67_Classes","BM_INFILTRATION","HIGH_LDH\n","PS_ECOG","HISTOLOGY","AA_STAGE","MARKER (0=NO; 1=IGH; 2=BCL1; 3=BOTH)","Ki67_Classes_1","MIPI_Classes","CLINICAL_SIT_PRE_ASCT","CLINICAL_SIT_postASCT")
for(tissue in c("Cluster_BM","Cluster_PB")){
  if(tissue=="Cluster_BM"){
    palette<-c("#F1BB7B", "#FD6467", "#5B1A18", "#D67236")
  }
  else if(tissue=="Cluster_PB"){
    palette<-c("#440154", "#31688E", "#35B779", "#8FD744")
  }
  for(variable in qual_vars){
    test<-test_indipendence(info_MCL0208_connector,variable,tissue,palette)
    if(is.null(test)){next}
    plot_independence<-test$plot+
      theme(
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent')
      )
    ggplot2::ggsave(plot=plot_independence,
                    path = "/Users/danielavolpatto/Documents/Università/dottorato/Connector/Plot connector/MCL0208/ClusterArmUniti/4Cluster/VariabiliCliniche",
                    filename = paste(variable,"_",tissue,".png",sep=""),
                    width = 18,
                    height = 10,
                    device="png",
                    units = "cm",
                    bg ='transparent')
  }
}

qual_vars<-c(qual_vars,"MIPI")
for(tissue in c("Cluster_BM","Cluster_PB")){
  if(tissue=="Cluster_BM"){
    palette<-c("#F79071","#994627")
    tibble<-info_MCL0208_connector%>%
      mutate(Cluster_BM=ifelse(Cluster_BM=="A"|Cluster_BM=="B","A/B","C/D"))
  }
  else if(tissue=="Cluster_PB"){
    palette<-c("#3B3571", "#62C75F")
    tibble<-info_MCL0208_connector%>%
      mutate(Cluster_PB=ifelse(Cluster_PB=="A"|Cluster_PB=="B","A/B","C/D"))
  }
  for(variable in qual_vars){
    test<-test_indipendence(tibble,variable,tissue,palette)
    plot_independence_unified<-test$plot+
      theme(
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent')
      )
    ggplot2::ggsave(plot=plot_independence_unified,
                    path = "/Users/danielavolpatto/Documents/Università/dottorato/Connector/Plot connector/MCL0208/ClusterArmUniti/ClusterUniti/VariabiliCliniche",
                    filename = paste(variable,"_",tissue,".png",sep=""),
                    width = 18,
                    height = 10,
                    device="png",
                    units = "cm",
                    bg ='transparent')
  }
}


mut<-Condensed_mut_Moia%>%
      group_by(ID)%>%
      mutate(across(contains("mutated"),~any(unlist(.),na.rm = TRUE)))%>%
      select(ID,contains("mutated"))
info_tmp<-full_join(info_MCL0208_connector,mut)
qual_vars<-colnames(mut[,2:ncol(mut)])
for(tissue in c("Cluster_BM","Cluster_PB")){
    if(tissue=="Cluster_BM"){
      palette<-c("#F1BB7B", "#FD6467", "#5B1A18", "#D67236")
    }
    else if(tissue=="Cluster_PB"){
      palette<-c("#440154", "#31688E", "#35B779", "#8FD744")
    }
    info_tmp<-info_tmp%>%
      filter(!is.na(!!sym(tissue))&!!sym(tissue)!="Unclassified")
    num_TRUE<-info_tmp%>%
      filter(!is.na(!!sym(tissue)))%>%
      select(qual_vars)%>%
      colSums(na.rm = TRUE)%>%
      as.vector()
    qual_vars_tmp<-qual_vars[num_TRUE>0]
    for(variable in qual_vars_tmp){
      test<-test_indipendence(info_tmp,variable,tissue,palette)
      if(is.null(test)) {next}
      plot_independence<-test$plot+
        theme(
          panel.background = element_rect(fill='transparent'),
          plot.background = element_rect(fill='transparent', color=NA),
          legend.background = element_rect(fill='transparent'),
          legend.box.background = element_rect(fill='transparent')
        )
      
      path<-paste("/Users/danielavolpatto/Documents/Università/dottorato/Connector/Plot connector/MCL0208/ClusterArmUniti/4Cluster/NsMut/",type,sep="")
      ggplot2::ggsave(plot=plot_independence,
                      path = path,
                      filename = paste("NsMut_",variable,"_",tissue,".png",sep=""),
                      width = 18,
                      height = 10,
                      device="png",
                      units = "cm",
                      bg ='transparent')
    }
}

quant_var<-c("flow_at_dia","Age","flow PB","flow BM","%_BM_INFILTR\n","MIPI","HB_LEVEL\n","NEUTRO_COUNT","LYMPHO_COUNT","PLTS","KI_67")
for(tissue in c("Cluster_BM","Cluster_PB")){
  if(tissue=="Cluster_BM"){
    palette<-c("#F1BB7B", "#FD6467", "#5B1A18", "#D67236")
  }
  else if(tissue=="Cluster_PB"){
    palette<-c("#440154", "#31688E", "#35B779", "#8FD744")
  }
  for(variable in quant_var){
    test<-test_distribution(info_MCL0208_connector,variable,tissue,palette)
    if(is.null(test)){next}
    plot<-test$plot
    stat_test<-test[names(test)!="plot"]
    null_stat_test<-sapply(stat_test,is.null)
    labs<-c("Kruskal-Wallis test for equal medians: ",
            "ANOVA test for equal means: ",
            "Bartlett's for equal variances: ",
            "Levene's test for equal variance: ")[!null_stat_test]
    pvalues<-c(stat_test$kw_test$p.value,
               summary(test$anova_test)[[1]]$'Pr(>F)'[1],
               test$bartlett_test$p.value,
               test$levene_test$`Pr(>F)`[1])
    caption<-paste(paste0(labs,round(pvalues,4),"\n"),collapse="")
    
    plot<-test$plot+
      labs(caption = paste("Shapiro test for normality: * = pvalue<0.05, ** = p-value<0.01 \n",
                           caption))
    ggplot2::ggsave(plot=plot,
                    path = "/Users/danielavolpatto/Documents/Università/dottorato/Connector/Plot connector/MCL0208/ClusterArmUniti/4Cluster/VariabiliCliniche",
                    filename = paste(variable,"_",tissue,".png",sep=""),
                    width = 18,
                    height = 10,
                    device="png",
                    units = "cm",
                    bg ='transparent')
  }
}

for(tissue in c("Cluster_BM","Cluster_PB")){
  for(variable in quant_var){
    if(tissue=="Cluster_BM"){
      palette<-c("#F1BB7B", "#FD6467", "#5B1A18", "#D67236")
      filtered_tibble<-info_MCL0208_connector%>%
        filter(!is.na(Cluster_BM)&!is.na(variable))%>%
        mutate(Cluster_unified=ifelse(Cluster_BM=="A"|Cluster_BM=="B","A/B","C/D"))
    }
    else if(tissue=="Cluster_PB"){
      palette<-c("#440154", "#31688E", "#35B779", "#8FD744")
      filtered_tibble<-info_Arm0_MCL0208_connector%>%
        filter(!is.na(Cluster_PB)&!is.na(variable))%>%
        mutate(Cluster_unified=ifelse(Cluster_PB=="A"|Cluster_PB=="B","A/B","C/D"))
    }
    test<-test_distr_unified(filtered_tibble,variable,"Cluster_unified",palette)
    
    plot<-test$plot
    if(!is.null(test$kw_test)){
      plot<-plot+
        labs(caption = paste("Shapiro test for normality: * = pvalue<0.05, ** = p-value<0.01 \n",
                             "Kolmogorov-Smirnov test for equal distribution: ",round(test$kw_test$p.value,4),"\n",
                             "T-test for equal means:",round(test$t_test$p.value,4),"\n",
                             "Wilcoxon-Mann-Whitney test:",round(test$wilcox_test$p.value,4),"\n",
                             "F-test for equal variances:",round(test$f_test$p.value,4),"\n",
                             "Levene's test for equal variances:", round(test$levene_test$`Pr(>F)`[1],4),"\n")
        )
    }
    ggplot2::ggsave(plot=plot,
                    path = "/Users/danielavolpatto/Documents/Università/dottorato/Connector/Plot connector/MCL0208/ClusterArmUniti/ClusterUniti/VariabiliCliniche",
                    filename = paste(variable,"_",tissue,".png",sep=""),
                    width = 18,
                    height = 10,
                    device="png",
                    units = "cm",
                    bg ='transparent')
    
  }
}

pharmagen_muts<-c("ABCB1 1236 C>T","ABCB1 2677 G>T,A","ABCB1 2677 G>T,A _2","ABCB1 3435 C,T","Aplotype ABCB1","VEGFA -2055 A>C","VEGFA -2055 A>C _2","ABCG2 421 C>A","FCGR2A 497 A>G","NCF4 -368 G>A","GSTP1 313 A>G","CRBN_1_rs1714327","CRBN_2_rs1705814") 
for(tissue in c("Cluster_BM","Cluster_PB")){
  if(tissue=="Cluster_BM"){
    palette<-c("#F1BB7B", "#FD6467", "#5B1A18", "#D67236")
  }
  else if(tissue=="Cluster_PB"){
    palette<-c("#440154", "#31688E", "#35B779", "#8FD744")
  }
  for(variable in pharmagen_muts){
    test<-test_indipendence(info_MCL0208_connector,variable,tissue,palette)
    if(is.null(test)){next}
    plot_independence<-test$plot+
      theme(
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent')
      )
    ggplot2::ggsave(plot=plot_independence,
                    path = "/Users/danielavolpatto/Documents/Università/dottorato/Connector/Plot connector/MCL0208/ClusterArmUniti/4Cluster/FarmacogenMut",
                    filename = paste(variable,"_",tissue,".png",sep=""),
                    width = 18,
                    height = 10,
                    device="png",
                    units = "cm",
                    bg ='transparent')
  }
}
for(tissue in c("Cluster_BM","Cluster_PB")){
  if(tissue=="Cluster_BM"){
    palette<-c("#F79071","#994627")
    tibble<-info_Arm0_MCL0208_connector%>%
      mutate(Cluster_BM=ifelse(Cluster_BM=="A"|Cluster_BM=="B","A/B","C/D"))
  }
  else if(tissue=="Cluster_PB"){
    palette<-c("#3B3571", "#62C75F")
    tibble<-info_Arm0_MCL0208_connector%>%
      mutate(Cluster_PB=ifelse(Cluster_PB=="A"|Cluster_PB=="B","A/B","C/D"))
  }
  for(variable in pharmagen_muts){
    test<-test_indipendence(tibble,variable,tissue,palette)
    plot_independence_unified<-test$plot+
      theme(
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent')
      )
    ggplot2::ggsave(plot=plot_independence_unified,
                    path = "/Users/danielavolpatto/Documents/Università/dottorato/Connector/Plot connector/MCL0208/ClusterArmUniti/ClusterUniti/FarmacogenMut",
                    filename = paste(variable,"_",tissue,".png",sep=""),
                    width = 18,
                    height = 10,
                    device="png",
                    units = "cm",
                    bg ='transparent')
  }
}

