
{
  setwd("/Users/danielavolpatto/Documents/Università/dottorato/Connector")
  
  library(dplyr)
  library(ggthemr)
  library(corrplot)
  library(tidyr)
  library(purrr)
  library(ggmosaic)
  #library(ggsankeyfier)
  library(readxl)
  library(stats)
  library(readr)
  library("vcd")
  library(car)
  library(ggplot2)
  library(patchwork)
  
  source("./inst/Rfunctions/Daniela/test_funct_MCL0208.R")
  
WholeData<-readRDS("./inst/Rfunctions/Daniela/WholeData.RDs")
load("./inst/Rfunctions/Daniela/Moia_preprocessing.RData")
load("./inst/Rfunctions/Daniela/Pavia_preprocessing.RData")

Variant_calling_mut_Pavia<-read.table("./inst/Rfunctions/Daniela/mut_Pavia_long_filt.tsv",sep=" ")
Variant_calling_mut_Moia<-read.table("./inst/Rfunctions/Daniela/mut_Moia_long_filt.tsv",sep=" ")
Condensed_mut_Pavia<-Variant_calling_mut_Pavia%>%
  filter(ACMG%in%c("Pathogenic","Likely Pathogenic"))%>%
  select(gene,ID,time,tissue)%>%
  distinct()%>%
  pivot_wider(names_from = "gene",
              names_prefix = "mutated_",
              values_from = "gene",
              values_fn = function(x){!is.na(x)})%>%
  mutate(ID=as.character(ID))%>%
  full_join(df_Pavia)

Condensed_mut_Pavia[is.na(Condensed_mut_Pavia)]<-FALSE
Condensed_mut_Pavia$time<-factor(Condensed_mut_Pavia$time,levels = c("DIAG","FU3","M1","M2"))


Condensed_mut_Moia<-Variant_calling_mut_Moia%>%
  filter(ACMG%in%c("Pathogenic","Likely Pathogenic"))%>%
  select(gene,ID)%>%
  distinct()%>%
  pivot_wider(names_from = "gene",
              names_prefix = "mutated_",
              values_from = "gene",
              values_fn = function(x){!is.na(x)})%>%
  mutate(ID=as.character(ID))%>%
  full_join(df_Moia)%>%
  rowwise()
Condensed_mut_Moia[is.na(Condensed_mut_Moia)]<-FALSE

# Condensed_mut_all<-bind_rows(Condensed_mut_Pavia%>%select(c("ID","time","tissue")),
#                          Condensed_mut_Moia%>%mutate(time="DIAG")%>%select(c("ID","time","tissue")))%>%
#   pivot_wider(names_from = time,
#               values_from = tissue,
#               values_fn = list)%>%
#   rowwise()%>%
#   mutate(DIAG=paste(DIAG,collapse=", "),
#          FU3=paste(FU3,collapse=", "),
#          M1=paste(M1,collapse=", "),
#          M2=paste(M2,collapse=", "))%>%
#   merge(full_join(Condensed_mut_Pavia%>%select(-c(time,tissue)),Condensed_mut_Moia%>%select(-tissue),by="ID")%>%
#           distinct()%>%
#           mutate(across(ends_with(".x"), ~ coalesce(., get(sub(".x$", ".y", cur_column()))))) %>%
#           rename_with(~ sub(".x$", "", .), ends_with(".x")) %>%
#           select(-ends_with(".y"))%>%
#           #mutate(na_count = rowSums(is.na(.)))%>%
#           group_by(ID) %>%
#           summarise_all(~ any(. == TRUE, na.rm = TRUE))%>%
#           distinct(),by="ID")

Condensed_mut_all<-bind_rows(Condensed_mut_Pavia,Condensed_mut_Moia%>%mutate(time="DIAG"))%>%
  group_by(ID,time)%>%
  mutate(across(c(tissue,contains("mutated")),~list(.)))%>%
  ungroup()%>%
  distinct()
Condensed_mut_all[Condensed_mut_all==""]<-NA
Condensed_mut_all$time<-factor(Condensed_mut_all$time,levels = c("DIAG","FU3","M1","M2"))


# prova<-Condensed_mut_Pavia%>%
#   pivot_wider(names_from = time,
#               values_from = tissue,
#               values_fn = list)%>%
#   rowwise()%>%
#   mutate(DIAG=ifelse(length(DIAG)>0,paste(DIAG,collapse=","),NA),
#          FU3=ifelse(length(FU3)>0,paste(FU3,collapse=","),NA),
#          M1=ifelse(length(M1)>0,paste(M1,collapse=","),NA),
#          M2=ifelse(length(M2)>0,paste(M2,collapse=","),NA))%>%
#   group_by(ID)%>%
#   mutate(across(!contains("ID"), ~ ifelse(sum(!is.na(.)>0),list(.[!is.na(.)]),NA)))%>%
#   distinct()


dataBM = merge(WholeData$MCL0208$BM$Dataset,
               WholeData$MCL0208$BM$CONNECTORList.FCM$CONNECTORList$LabCurv %>% select(IDold,ID))%>% 
  select(-ID) %>% rename(ID = IDold)

dataPB = merge(WholeData$MCL0208$PB$Dataset,
               WholeData$MCL0208$PB$CONNECTORList.FCM$CONNECTORList$LabCurv %>% select(IDold,ID))%>% 
  select(-ID) %>% rename(ID = IDold)

ClusterBM<-dataBM%>%
  select(ID,Cluster,Arm,TTP)%>%
  rename("Cluster_BM"="Cluster")

ClusterPB<-dataPB%>%
  select(ID,Cluster,Arm,TTP)%>%
  rename("Cluster_PB"="Cluster")

Cluster<-full_join(ClusterBM,ClusterPB)%>%
  select(ID,Arm,TTP,Cluster_BM,Cluster_PB)
#intersect(intersect(c(ClusterPB$ID,ClusterBM$ID),df_Moia$ID),df_Pavia$ID)

ID_comuni<-intersect(c(ClusterPB$ID,ClusterBM$ID),df_Pavia$ID)

ID_chip_TP53mutdel_age <- read_excel("./inst/Rfunctions/Daniela/ID_chip_TP53mutdel_age.xlsx")
ID_chip_TP53mutdel_age[ID_chip_TP53mutdel_age=="NA"]<-NA

info_pz_MCL0208<-full_join(ID_chip_TP53mutdel_age,Cluster)
info_pz_MCL0208$ID<-as.character(info_pz_MCL0208$ID)

df_Pavia$tissue<-as.character(df_Pavia$tissue)
df_Moia$tissue<-as.character(df_Moia$tissue)

mut_info<-rbind(df_Pavia,df_Moia%>%mutate(time="DIAG"))%>%
  pivot_wider(names_from = time,
              values_from = tissue,
              values_fn = list)%>%
  rowwise()%>%
  mutate(DIAG=paste(DIAG,collapse=", "),
         FU3=paste(FU3,collapse=", "),
         M1=paste(M1,collapse=", "),
         M2=paste(M2,collapse=", "))
mut_info[mut_info==""]<-NA

info_pz_MCL0208<-full_join(mut_info,info_pz_MCL0208)

info_pz_MCL0208_1 <- read_csv("./inst/Rfunctions/Daniela/MCL0208_database- Gennaio 2022 BA.xlsx - Foglio1.csv")
info_pz_MCL0208_1<-info_pz_MCL0208_1 %>% 
  select(which(info_pz_MCL0208_1[1, ] == "YES"),
         -AGE_DIA)%>%
  rename("ID"="code")
info_pz_MCL0208_1<-info_pz_MCL0208_1[-1,]
info_pz_MCL0208_1<-info_pz_MCL0208_1[rowSums(is.na(info_pz_MCL0208_1)) != ncol(info_pz_MCL0208_1), ]

info_pz_MCL0208<-full_join(info_pz_MCL0208,info_pz_MCL0208_1)

npat<-info_pz_MCL0208%>%nrow()
npat_connector<-info_pz_MCL0208%>%filter(!is.na(Cluster_BM)|!is.na(Cluster_PB))%>%nrow()
npat_connector_PB<-info_pz_MCL0208%>%filter(!is.na(Cluster_PB))%>%nrow()
npat_connector_BM<-info_pz_MCL0208%>%filter(!is.na(Cluster_BM))%>%nrow()
#npat_seq<-info_pz_MCL0208%>%filter(!is.na(DIAG)|!is.na(FU3)|!is.na(M1)|!is.na(M2))%>%nrow()

info_pz_MCL0208_connector<-info_pz_MCL0208%>%filter(!is.na(Cluster_BM)|!is.na(Cluster_PB))

info_pz_MCL0208_connector<-info_pz_MCL0208_connector%>%
  mutate(Ki67_Classes_1=ifelse(as.numeric(KI_67)<30,0,1))
info_pz_MCL0208_connector<-info_pz_MCL0208_connector%>%
  mutate(MIPI_Classes=ifelse(MIPI< 5.70,"low", ifelse(MIPI<6.20,"intermediate", "high")))
info_pz_MCL0208_connector$MIPI_Classes<-as.factor(info_pz_MCL0208_connector$MIPI_Classes)
levels(info_pz_MCL0208_connector$MIPI_Classes)<-c("low","intermediate","high")

## save Informations

saveRDS(info_pz_MCL0208_connector,file = "inst/Data/info_pz_MCL0208_connector.RDs")
Condensed_mut = list(All = Condensed_mut_all, Moia = Condensed_mut_Moia, Pavia =Condensed_mut_Pavia )
saveRDS(Condensed_mut,file = "inst/Data/Condensed_mut.RDs")



qual_vars<-c("PFS","OS","TTP","TP53_loss_or_mut","TP53","TP53_loss_array","TP53_disruption","Arm","Ki67_Classes","BM_INFILTRATION","HIGH_LDH\n","PS_ECOG","HISTOLOGY","AA_STAGE","MARKER (0=NO; 1=IGH; 2=BCL1; 3=BOTH)","Ki67_Classes_1","MIPI_Classes")
for(tissue in c("BM","PB")){
  for(variable in qual_vars){
    test<-test_indipendence(info_pz_MCL0208_connector,variable,tissue)
    plot_independence<-test$plot+
      theme(
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent')
      )
    print(plot_independence)
    # ggplot2::ggsave(plot=plot_independence,
    #                 path = "/Users/danielavolpatto/Documents/Università/dottorato/MCLVariantCalling/Plot connector/Qualitative/Variabili Cliniche/4 Cluster",
    #                 filename = paste("plot_independence_",variable,"_",tissue,".png",sep=""),
    #                 width = 18,
    #                 height = 10,
    #                 device="png",
    #                 units = "cm",
    #                 bg ='transparent')
  }
}

qual_vars<-c(qual_vars,"MIPI")
for(tissue in c("BM","PB")){
  for(variable in qual_vars){
    test<-test_unified(info_pz_MCL0208_connector,variable,tissue)
    plot_independence_unified<-test$plot+
      theme(
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent')
      )
    ggplot2::ggsave(plot=plot_independence_unified,
                    path = "/Users/danielavolpatto/Documents/Università/dottorato/MCLVariantCalling/Plot connector/Qualitative/Variabili Cliniche/Cluster Uniti",
                    filename = paste("plot_independence_unified_",variable,"_",tissue,".png",sep=""),
                    width = 18,
                    height = 10,
                    device="png",
                    units = "cm",
                    bg ='transparent')
  }
}

cond_type<-c("_all_anyTimeAndTissue","Pavia_DIAG","_Pavia_anyTimeAndTissue","_Pavia_lost","_Pavia_gained","_Moia")
for(type in cond_type){
  if(type=="_all_anyTimeAndTissue"){
    mut<-Condensed_mut_all%>%
      rowwise()%>%
      mutate(across(contains("mutated"),~any(unlist(.),na.rm = TRUE)))%>%
      select(ID,contains("mutated"))
  }
  else if(type=="_Pavia_anyTimeAndTissue"){
    mut<-Condensed_mut_Pavia%>%
      rowwise()%>%
      mutate(time_tissue=paste(time,tissue,sep="_"))%>%
      select(-c(time,tissue))%>%
      group_by(ID)%>%
      mutate(across(c(time_tissue,contains("mutated")),~list(.)))%>%
      ungroup()%>%
      distinct()%>%
      rowwise()%>%
      mutate(across(contains("mutated"),~any(unlist(.),na.rm = TRUE)))%>%
      select(ID,contains("mutated"))
  }
  else if(type=="_Moia"){
    mut<-Condensed_mut_Moia%>%
      group_by(ID)%>%
      mutate(across(contains("mutated"),~any(unlist(.),na.rm = TRUE)))%>%
      select(ID,contains("mutated"))
  }
  else if(type=="_PaviaDIAG"){
    mut<-Condensed_mut_Pavia%>%
      filter(time=="DIAG")%>%
      select(-time)%>%
      group_by(ID)%>%
      mutate(across(contains("mutated"),~any(.,na.rm = TRUE)))%>%
      distinct()%>%
      select(ID,contains("mutated"))
  }
  else if(type=="_Pavia_lost"){
    mut<-Condensed_mut_Pavia%>%
      order(time)
      select(-tissue)%>%
      group_by(ID,time)%>%
      mutate(across(contains("mutated"),~any(.)))%>%
      distinct()%>%
      ungroup()%>%
      select(-time)%>%
      group_by(ID)%>%
      mutate(across(contains("mutated"),~ifelse(is.unsorted(.),TRUE,FALSE)))%>%
      distinct()
  }
  else if(type=="_Pavia_gained"){
    mut<-Condensed_mut_Pavia%>%
      select(-tissue)%>%
      group_by(ID,time)%>%
      mutate(across(contains("mutated"),~any(.)))%>%
      distinct()%>%
      ungroup()%>%
      select(-time)%>%
      group_by(ID)%>%
      mutate(across(contains("mutated"),~ifelse(!is.unsorted(.)&length(unique(.))>1,TRUE,FALSE)))%>%
      distinct()
  }
  info_tmp<-full_join(info_pz_MCL0208_connector,mut)
  qual_vars<-colnames(mut[,2:ncol(mut)])
  for(tissue in c("BM","PB")){
    num_TRUE<-info_tmp%>%
      filter(!is.na(get(paste("Cluster",tissue,sep = "_"))))%>%
      select(qual_vars)%>%
      colSums(na.rm = TRUE)%>%
      as.vector()
    qual_vars_tmp<-qual_vars[num_TRUE>0]
    for(variable in qual_vars_tmp){
      test<-test_indipendence(info_tmp,variable,tissue)
      print(paste(variable,tissue))
      print(test$table%>%sum())
      print(info_tmp%>%filter(!is.na(get(paste("Cluster",tissue,sep = "_"))),!is.na(get(variable)))%>%nrow())
      #prova<-info_pz_MCL0208_connector%>%filter(!is.na(get(paste("Cluster",tissue,sep = "_"))),!is.na(get(variable)))%>%select(Cluster_BM,variable)
      
      plot_independence<-test$plot+
        theme(
          panel.background = element_rect(fill='transparent'),
          plot.background = element_rect(fill='transparent', color=NA),
          legend.background = element_rect(fill='transparent'),
          legend.box.background = element_rect(fill='transparent')
        )
      if(grepl("Pavia",type)){folder<-"Solo Pavia/"}
      else if(grepl("Moia",type)){folder<-"Solo Moia/"}
      else{folder<-"Almeno uno/"}
      path<-paste("/Users/danielavolpatto/Documents/Università/dottorato/MCLVariantCalling/Plot connector/Qualitative/Geni Mutati/",folder,sep="")
      ggplot2::ggsave(plot=plot_independence,
                      path = path,
                      filename = paste("plot_independence_",variable,"_",tissue,type,".png",sep=""),
                      width = 18,
                      height = 10,
                      device="png",
                      units = "cm",
                      bg ='transparent')
    }
  }
}

info_TP53<-info_pz_MCL0208_connector%>%select(ID,DIAG,FU3,M1,M2,Cluster_BM,Cluster_PB,contains("TP53"))
write.csv(info_TP53,"info_TP53.csv")

info_TP53_table<-info_TP53%>%
  mutate(mut_info=any(!is.na(c(DIAG,FU3,M1,M2))),
         mut_info=ifelse(mut_info,mut_info,NA))%>%
  pivot_longer(c(Cluster_BM,Cluster_PB),names_to = "tissue",values_to = "Cluster")%>%
  select(-c(ID,DIAG,FU3,M1,M2))%>%
  group_by(tissue)%>%
  summarize(
    tot=sum(!is.na(Cluster)),
    across(c(mut_info,contains("TP53")),~sum(!is.na(.)&!is.na(Cluster))))

write.csv(info_TP53_table,"info_TP53_table.csv")

quant_var<-c("flow_at_dia","Age","flow PB","flow BM","%_BM_INFILTR\n","MIPI","HB_LEVEL\n","NEUTRO_COUNT","LYMPHO_COUNT","PLTS","KI_67")
  
for(tissue in c("BM","PB")){
  for(variable in quant_var){
    test<-test_distribution(info_pz_MCL0208_connector,variable,tissue)
    plot<-test$plot+
      labs(caption = paste("Shapiro test for normality: * = pvalue<0.05, ** = p-value<0.01 \n",
                           "Kruskal-Wallis test for equal medians: ",round(test$kw_test$p.value,4),"\n",
                           "ANOVA test for equal means: ",round(summary(test$anova_test)[[1]]$'Pr(>F)'[1],4),"\n",
                           "Levene's test for equal variance:", round(test$levene_test$`Pr(>F)`[1],4)))
    
    ggplot2::ggsave(plot=plot,
                    path = "/Users/danielavolpatto/Documents/Università/dottorato/Connector/Plot connector/MCL0208/Quantitative/4 Cluster/",
                    filename = paste("plot_distribution_",variable,"_",tissue,".png",sep=""),
                    width = 18,
                    height = 10,
                    device="png",
                    units = "cm",
                    bg ='transparent')

  }
}

for(tissue in c("BM","PB")){
  for(variable in quant_var){
    test<-test_distr_unified(info_pz_MCL0208_connector,variable,tissue)
    plot<-test$plot+
      labs(caption = paste("Shapiro test for normality: * = pvalue<0.05, ** = p-value<0.01 \n",
                           "Kolmogorov-Smirnov test for equal distribution: ",round(test$kw_test$p.value,4),"\n",
                           "Wilcoxon-Mann-Whitney test:",round(test$wilcox_test$p.value,4),"\n",
                           "Levene's test for equal variances:", round(test$levene_test$`Pr(>F)`[1],4),"\n")
      )
    ggplot2::ggsave(plot=plot,
                    path = "/Users/danielavolpatto/Documents/Università/dottorato/Connector/Plot connector/MCL0208/Quantitative/Cluster uniti/",
                    filename = paste("plot_distribution_unified",variable,"_",tissue,".png",sep=""),
                    width = 18,
                    height = 10,
                    device="png",
                    units = "cm",
                    bg ='transparent')
    
  }
}
