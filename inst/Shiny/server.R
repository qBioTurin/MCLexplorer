WholeData <- readRDS( system.file("Data","WholeData.RDs", package = "ShinyMCL") )
SurvYOUNGERData <- readRDS( system.file("Data","YoungerMCLsurvival.RDs", package = "ShinyMCL") )
ColorCluster <- readRDS( system.file("Data","ColorCluster.RDs", package = "ShinyMCL") )
Condensed_mut <- readRDS( system.file("Data","Condensed_mut.RDs", package = "ShinyMCL") )
info_pz_MCL0208_connector <- readRDS( system.file("Data","info_pz_MCL0208_connector.RDs", package = "ShinyMCL") )
source( system.file("Shiny","AuxFunctions.R", package = "ShinyMCL") )
source( system.file("Rfunctions","test_funct_MCL0208.R", package = "ShinyMCL") )

axistheme = theme(
  plot.caption = element_text(size = 12, face = "italic", family = "Times"),
  axis.text.x = element_text(size = 14, family = "Times"),
  axis.text.y = element_text(size = 14, family = "Times"),
  axis.title.x = element_text(margin = margin(b = 1, unit = "cm"), size = 16, face = "bold", family = "Times"),
  axis.title.y = element_text(size = 16, face = "bold", family = "Times"),
  legend.text = element_text(size = 14, family = "Times"),
  legend.title = element_text(size = 16, face = "bold", family = "Times")
)

server <- function(input, output, session) {
  
  ### MCL208 exploration ####
  
  MCLvalues = reactiveValues(EntropyDF = NULL, ClusterDF = NULL , splinePlots = NULL, info_pz_MCL0208 = NULL)
  
  observe({
    show_modal_spinner()
    
    tissueMCL= input$tissueMCL
    #ArmMCL = input$ArmMCL
    ArmMCL = "NULL"
    
    isolate({
      
      if(ArmMCL == "NULL"){
        CONNECTORList.FCM = WholeData$MCL0208[[tissueMCL]]$CONNECTORList.FCM
      }else{
        CONNECTORList.FCM = WholeData$MCL0208[[tissueMCL]]$ClassSperimClusterOsserav[[paste0("Cluster",ArmMCL)]]
      }
      
      MCLvalues$splinePlots = connector::Spline.plots(FCM.plots = ClusterWithMeanCurve(CONNECTORList.FCM,feature = "Dynamic") )
      
      ## Entropy plot
      probs = CONNECTORList.FCM$FCM$cluster$cl.info$probs
      colnames(probs) =CONNECTORList.FCM$FCM$cluster$cluster.names
      Dataset = CONNECTORList.FCM$CONNECTORList$Dataset
      
      MatrixClass= as.data.frame(probs)
      MatrixClass$ClusterType <- colnames(MatrixClass)[apply(MatrixClass, MARGIN = 1, FUN = which.max)]
      MatrixClass= MatrixClass %>%
        mutate(MajorClusterValue = pmax(A,B,C,D))
      
      ProbMax = as.data.frame(t(sapply(1:dim(probs)[1], function(x) data.frame(ID = x, probMax = max(probs[x,])  ) )))
      
      df1 = merge(
        Dataset %>%
          group_by(ID) %>%
          dplyr::summarise(L = length(Time)),
        MatrixClass %>%
          mutate(ID = 1:length(ClusterType) ) %>%
          tidyr::gather(-ID,-MajorClusterValue,-ClusterType,key = "Cluster",value = "Prob") %>%
          group_by(ID) %>%
          mutate(Entropy = -sum(Prob*log2(Prob))) %>%
          ungroup() %>%
          tidyr::spread(key = "Cluster",value = "Prob") 
      )
      
      
      ## Curve Plot
      # means curve 
      
      meancurves = as.data.frame(CONNECTORList.FCM$FCM$prediction$meancurves)
      colnames(meancurves) = CONNECTORList.FCM$FCM$cluster$cluster.names
      meancurves$Time =CONNECTORList.FCM$CONNECTORList$TimeGrid
      meancurves = meancurves %>% tidyr::gather(-Time, key = "Cluster", value = "Observation")
      
      df = merge(
        merge(
          CONNECTORList.FCM$FCM$cluster$ClustCurve,
          MatrixClass %>%
            mutate(ID = 1:length(ClusterType) ) %>%
            tidyr::gather(-ID,-MajorClusterValue, -ClusterType, key = "Cluster", value = "Prob") %>%
            group_by(ID) %>%
            mutate(Entropy = -sum(Prob*log2(Prob))) %>%
            ungroup() %>%
            tidyr::spread(key = "Cluster",value = "Prob") ),
        CONNECTORList.FCM$CONNECTORList$Dataset
      ) %>%
        mutate(Cluster = CONNECTORList.FCM$FCM$cluster$cluster.names[Cluster] ) %>%
        group_by(ID) %>%
        dplyr::mutate(L = length(Time) )
      
      df = merge(df,CONNECTORList.FCM$CONNECTORList$LabCurv %>% select(ID,IDold),by="ID" )
      
      MCLvalues$EntropyDF = df1
      MCLvalues$ClusterDF = df
      MCLvalues$meancurves = meancurves
      
      maxL = max(df1$L)
      maxE = max(df1$Entropy)
      minL = min(df1$L)
      minE = min(df1$Entropy)
      
      updateSliderInput(session = session, inputId = "entropy", min = minE, max = maxE, value = c(minE,maxE) )
      updateSliderInput(session = session, inputId =  "length", min = minL, max = maxL, value = c(minL,maxL) )
      
    })
    
    remove_modal_spinner()
  })
  
  output$MCL_clusteringPlot <- renderPlot({
    req(MCLvalues$ClusterDF -> df)
    input$tissueMCL -> tissueMCL
    WholeData[["MCL0208"]][[tissueMCL]]$Dataset -> Dataset
    completeDF = merge(df,Dataset%>%select(-Cluster), by = "ID")
    Tissuecol = ColorCluster[grep(pattern = paste0(tissueMCL,"_"), names(ColorCluster))]
    names(Tissuecol) = gsub(pattern = paste0(tissueMCL,"_"),replacement = "",names(Tissuecol))
    PlotComplete = plot.genaration(completeDF,Tissuecol)
    grid.draw(PlotComplete)
  })
  
  output$MCL_survPlot <- renderPlot({
    req(MCLvalues$ClusterDF -> df)
    selectSurvMCL = input$selectSurvMCL
    input$tissueMCL -> tissueMCL
    # input$ArmMCL -> arm
    arm = "NULL"
    WholeData[["MCL0208"]][[tissueMCL]]$Dataset -> Dataset
    
    if(arm != "NULL"){
      Dataset = Dataset %>% filter(Arm == arm)
    }
    
    Dataset = Dataset %>% mutate(ArmLabel = ifelse(Arm == 1,"Lenalidomide Therapy", "No Therapy"))
    
    if(selectSurvMCL == "Arm - single cluster"){
      ### single Clusters w.r.t ARM
      ggsurvList = lapply( c(as.list(unique(Dataset$Cluster)),list(c("A","B"), c("C","D")) ), function(cl,Dataset){
        
        Dataset_cl = Dataset %>% filter(Cluster %in% cl)
        
        fit <- eval(parse(text = paste0("survfit(Surv(Dataset_cl$TimeTTPevent,Dataset_cl$TTP) ~ ArmLabel, data = Dataset_cl)")))
        list2env(list(Dataset_cl = Dataset_cl), envir = .GlobalEnv)
        
        ggsurv=ggsurvplot(fit = fit, data = Dataset_cl,
                          xlab = "Days",  ylab = "TTP",
                          size = 1, pval = TRUE, risk.table = TRUE,conf.int = T,
                          risk.table.col="strata",ggtheme = theme_bw() ,surv.median.line = "hv" )
        
        ggsurv$plot = ggsurv$plot + labs(title = paste0("Cluster ", paste0(cl,collapse = " "))  )
        
        surv_median_values = surv_median(fit = fit)
        if(! all(is.na(surv_median_values$median)))
          ggsurv$plot = ggsurv$plot +
          annotate("text", x = surv_median_values$median, y = 0,
                   label = round(surv_median_values$median/365,digits = 3), size = 3, hjust = 0)
        
      },Dataset)
      
      
      pl = cowplot::plot_grid(plotlist = ggsurvList, ncol = 2)
    }else{
      list2env(list(Dataset = Dataset), envir = .GlobalEnv)
      
      if(selectSurvMCL == "Arm - all clusters"){
        fit <- eval(parse(text = paste0("survfit(Surv(Dataset$TimeTTPevent,Dataset$TTP) ~ ArmLabel, data = Dataset)")))
      }else if(selectSurvMCL == "Clusters"){
        fit <- eval(parse(text = paste0("survfit(Surv(Dataset$TimeTTPevent,Dataset$TTP) ~ Cluster, data = Dataset)")))
      }else if(selectSurvMCL == "Merging Clusters"){
        Dataset = Dataset %>% mutate(Cluster2 = ifelse(Cluster %in% c("A","B"), "A-B", "C-D"))
        if( length(unique(Dataset$Cluster)) > 4 ) 
          Dataset = Dataset %>% 
            mutate(Cluster2 = ifelse(Cluster2 != "A-B", paste0(LETTERS[3:length(unique(Dataset$Cluster))],collapse = "-") ))
        
        fit <- eval(parse(text = paste0("survfit(Surv(Dataset$TimeTTPevent,Dataset$TTP) ~ Cluster2, data = Dataset)")))
      } 
      
      ggsurv=ggsurvplot(fit = fit, data = Dataset,
                        xlab = "Days",  ylab = "TTP",
                        size = 1, pval = TRUE, risk.table = TRUE,conf.int = T,
                        risk.table.col="strata",ggtheme = theme_bw() ,surv.median.line = "hv" )
      
      surv_median_values = surv_median(fit = fit)
      if(! all(is.na(surv_median_values$median)))
        ggsurv$plot = ggsurv$plot +
        annotate("text", x = surv_median_values$median, y = 0,
                 label = round(surv_median_values$median,digits = 3), size = 3, hjust = 0)
      
      ggsurv$plot -> pl
    }
    
    pl
  })
  
  output$MCL_survUIPlot <- renderUI({
    selectSurvMCL = input$selectSurvMCL
    if(selectSurvMCL == "Arm - single cluster") plot_height = "800px"
    else plot_height = "400px"
    
    plotOutput("MCL_survPlot", height = plot_height)
  })
  
  output$MCL_scatterPlot <- renderPlot({
    if(!is.null(MCLvalues$EntropyDF) ){
      df1_sub  = MCLvalues$EntropyDF %>%
        filter(
          Entropy <= max(input$entropy), Entropy >= min(input$entropy),
          L <= max(input$length), L >= min(input$length)
        )
      
      updateSelectInput("SplineID", session = session,
                        choices = unique(df1_sub$ID) )
      
      if(dim(df1_sub)[1] > 0 ){
        maxL = max(MCLvalues$EntropyDF$L)
        maxE = max(MCLvalues$EntropyDF$Entropy)
        minL = min(MCLvalues$EntropyDF$L)
        minE = min(MCLvalues$EntropyDF$Entropy)
        
        plot =  df1_sub %>%
          ggplot( ) +
          geom_text(aes(y = Entropy, x = L,col = Entropy, label = ID ))+
          ylim(0,maxE+0.1) +
          xlim(minL-1,maxL+1) + 
          geom_jitter(aes(y = Entropy, x = L,col = Entropy ),
                      width = 0.1,height = 0) +
          facet_wrap(~ClusterType) + 
          theme_bw() +
          labs(x = "Curve Length", y = "Entropy", col = "Entropy") +
          theme(legend.position = "right")+
          scale_color_continuous(low = "grey",
                                 high = "#53190b", 
                                 breaks = c(round(minE,digits = 2),
                                            round(maxE,digits = 2)),
                                 labels = c(round(minE,digits = 2),
                                            round(maxE,digits = 2))
          )
      }else{
        plot = ggplot()
      }
      
      plot
    }
  })
  
  output$MCL_linePlot  <- renderPlot({
    if(!is.null(MCLvalues$ClusterDF) ){
      df_sub = MCLvalues$ClusterDF %>% 
        filter(Entropy <= max(input$entropy),Entropy >= min(input$entropy),
               L <= max(input$length),L >= min(input$length))
      
      if(dim(df_sub)[1] > 0 ){ 
        
        plot = ggplot(df_sub) +
          geom_line(data = MCLvalues$meancurves, aes(x = Time, y = Observation), col = "red" )+
          geom_line(aes(x = Time, y = Observation,group = ID, col = as.numeric(Entropy)))+ 
          facet_wrap(~Cluster)+ 
          theme_bw() +
          labs(col = "Entropy") +
          theme(legend.position = "none")+
          scale_color_continuous(low = "grey",high = "#53190b")
        
      }else{
        plot = ggplot()
      }
      
      plot
    }
  })
  
  output$MCL_FittedCurve <- renderPlot({
    ID <- input$SplineID
    MCLvalues$splinePlots[[ID]]
  })
  
  
  #### Statistical part 
  output$DTtableMCL <- DT::renderDT({
    req(MCLvalues$ClusterDF -> df)
    input$tissueMCL -> tissue
    # input$ArmMCL -> arm
    arm = "NULL"
    input$qual_varsMCL -> variable
    
    if(arm == "NULL"){
      pz_data = info_pz_MCL0208_connector %>%
        mutate(Cluster = !!sym(paste0("Cluster_",tissue)) )%>%
        select(-Cluster_BM,-Cluster_PB)
    }else{
      dftmp = df %>% select(-ID) %>% mutate( ID = paste0(IDold) )  %>% select(ID,Cluster) %>% distinct()
      pz_data = merge(info_pz_MCL0208_connector %>% select(-Cluster_BM,-Cluster_PB) ,dftmp)
    }
    
    pz_data = pz_data  %>% 
      relocate(ID, Cluster,PFS, Arm) %>%
      mutate(ID = as.factor(ID),
             Cluster = as.factor(Cluster))
    
    DT::datatable(pz_data, filter = 'top',
                  options = list(
                    pageLength = 5, autoWidth = TRUE,scrollX = TRUE
                  )
    )
    
  })
  
  output$MCL_qualClinicalPlot <- renderPlot({
    req(MCLvalues$ClusterDF -> df)
    input$ClusterCheckMCL
    input$tissueMCL -> tissue
    # input$ArmMCL -> arm
    arm = "NULL"
    input$qual_varsMCL -> variable
    
    if(arm == "NULL"){
      pz_data = info_pz_MCL0208_connector
    }else{
      dftmp = df %>% select(-ID) %>% mutate( ID = paste0(IDold) )  %>% select(ID,Cluster) %>% distinct()
      tmp = merge(info_pz_MCL0208_connector,dftmp)
      pz_data = tmp %>% mutate(Cluster_BM = Cluster, Cluster_PB = Cluster) 
    }
    
    if(input$ClusterCheckMCL){
      pz_data[,paste0("Cluster_",tissue)]<- ifelse( pz_data[[paste0("Cluster_",tissue)]] %in% c("A","B"),"A/B","C/D")
      if(tissue=="BM"){
        palette<-c("#F1BB7B", "#FD6467", "#5B1A18", "#D67236")
      }
      else if(tissue=="PB"){
        palette<-c("#440154", "#31688E", "#35B779", "#8FD744")
      }
    }else{
      if(tissue=="BM"){
        palette<-c("#F1BB7B", "#FD6467", "#5B1A18", "#D67236")
      }
      else if(tissue=="PB"){
        palette<-c("#440154", "#31688E", "#35B779", "#8FD744")
      }
    }
    
    test<-test_indipendence(pz_data,variable,paste0("Cluster_",tissue),palette)
    
    test$plot & (
      theme(
        legend.position = "bottom",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent')
      ) +
        axistheme 
    )
    
  })
  
  output$MCL_quantClinicalPlot <- renderPlot({
    req(MCLvalues$ClusterDF -> df)
    input$tissueMCL -> tissue
    input$quant_varsMCL -> variable

    # input$ArmMCL -> arm
    arm = "NULL"
    if(arm == "NULL"){
      pz_data = info_pz_MCL0208_connector
    }else{
      dftmp = df %>% select(-ID) %>% mutate( ID = paste0(IDold) )  %>% select(ID,Cluster) %>% distinct()
      tmp = merge(info_pz_MCL0208_connector,dftmp)
      pz_data = tmp %>% mutate(Cluster_BM = Cluster, Cluster_PB = Cluster) 
    }
    
    
    if(input$ClusterCheckMCL){ 
      
      if(tissue=="BM"){
        palette<-c("#F79071","#994627")
      }
      else if(tissue=="PB"){
        palette<-c("#3B3571", "#62C75F")
      }
      
      pz_data[,paste0("Cluster_",tissue)]<- ifelse( pz_data[[paste0("Cluster_",tissue)]] %in% c("A","B"),"A/B","C/D")
      test<-test_distr_unified(pz_data,variable,paste0("Cluster_",tissue),palette)
      plot<-test$plot+
        labs(caption = paste("Shapiro test for normality: * = pvalue<0.05, ** = p-value<0.01 \n",
                             "Kolmogorov-Smirnov test for equal distribution: ",round(test$kw_test$p.value,4),"\n",
                             "T-test for equal means:",round(test$t_test$p.value,4),"\n",
                             "Wilcoxon-Mann-Whitney test:",round(test$wilcox_test$p.value,4),"\n",
                             "F-test for equal variances:",round(test$f_test$p.value,4),"\n",
                             "Levene's test for equal variances:", round(test$levene_test$`Pr(>F)`[1],4),"\n")
        )
    }else{ 
      if(tissue=="BM"){
        palette<-c("#F1BB7B", "#FD6467", "#5B1A18", "#D67236")
      }
      else if(tissue=="PB"){
        palette<-c("#440154", "#31688E", "#35B779", "#8FD744")
      }
      
      test<-test_distribution(pz_data,variable,paste0("Cluster_",tissue),palette)
      
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
    }
    
    plot & (theme(
      legend.position = "bottom")+axistheme)
    
  })
  
  observe({
    req(MCLvalues$ClusterDF -> df)
    input$tissueMCL -> tissue
    # input$cond_typeMCL -> type
    type = "MCL panel"
    # input$ArmMCL -> arm
    arm = "NULL"
    oldtypes = c("All_anyTimeAndTissue","Pavia_DIAG","Pavia_anyTimeAndTissue","Pavia_lost","Pavia_gained","Moia")
    names(oldtypes) = c("CHIP panel - anyTimeAndTissue","CHIP panel - DIAGNOSIS","CHIP panel - anyTimeAndTissue","CHIP panel - lost","CHIP panel - gained","MCL panel")
    oldtypes[type]-> type
    
    if(type=="All_anyTimeAndTissue"){
      mut<-Condensed_mut$All %>%
        rowwise()%>%
        mutate(across(contains("mutated"),~any(unlist(.),na.rm = TRUE)))%>%
        select(ID,contains("mutated"))
    }
    else if(type=="Pavia_anyTimeAndTissue"){
      mut<- Condensed_mut$Pavia %>%
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
    else if(type=="Moia"){
      mut<-Condensed_mut$Moia%>%
        group_by(ID)%>%
        mutate(across(contains("mutated"),~any(unlist(.),na.rm = TRUE)))%>%
        select(ID,contains("mutated"))
    }
    else if(type=="Pavia_DIAG"){
      mut<-Condensed_mut$Pavia%>%
        filter(time=="DIAG")%>%
        select(-time)%>%
        group_by(ID)%>%
        mutate(across(contains("mutated"),~any(.,na.rm = TRUE)))%>%
        distinct()%>%
        select(ID,contains("mutated"))
    }
    else if(type=="Pavia_lost"){
      mut<-Condensed_mut$Pavia%>%
        #order(time)
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
    else if(type=="Pavia_gained"){
      mut<-Condensed_mut$Pavia%>%
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
    
    if(arm == "NULL"){
      pz_data = info_pz_MCL0208_connector
    }else{
      dftmp = df %>% select(-ID) %>% mutate( ID = paste0(IDold) )  %>% select(ID,Cluster) %>% distinct()
      tmp = merge(info_pz_MCL0208_connector,dftmp)
      pz_data = tmp %>% mutate(Cluster_BM = Cluster, Cluster_PB = Cluster) 
    }
    
    info_tmp<-full_join(pz_data,mut)
    qual_vars<-colnames(mut[,2:ncol(mut)])
    
    num_TRUE<-info_tmp%>%
      filter(!is.na(get(paste("Cluster",tissue,sep = "_"))))%>%
      select(qual_vars)%>%
      colSums(na.rm = TRUE)%>%
      as.vector()
    
    updateSelectInput("var_typeMCL" , session = session, choices =  qual_vars[num_TRUE>0] )
    MCLvalues$info_pz_MCL0208 = info_tmp
  })
  
  output$MCL_MutationPlot <- renderPlot({
    req(MCLvalues$info_pz_MCL0208 -> info_tmp)
    req(input$tissueMCL -> tissue)
    req(input$var_typeMCL -> variable)
    
    if(input$ClusterCheckMCL){
      if(tissue=="BM"){
        palette<-c("#F79071","#994627")
      }
      else if(tissue=="PB"){
        palette<-c("#3B3571", "#62C75F")
      }
      info_tmp[,paste0("Cluster_",tissue)]<- ifelse( info_tmp[[paste0("Cluster_",tissue)]] %in% c("A","B"),"A/B","C/D")
    }
    else{
      if(tissue=="BM"){
        palette<-c("#F1BB7B", "#FD6467", "#5B1A18", "#D67236")
      }
      else if(tissue=="PB"){
        palette<-c("#440154", "#31688E", "#35B779", "#8FD744")
      }
    }
    
    test<-test_indipendence(info_tmp,variable,paste0("Cluster_",tissue),palette)
    test$plot  & (
      theme(
        legend.position = "bottom",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent')
      )+axistheme
    )
  })
  
  output$MCL_pharmaPlot <- renderPlot({
    req(MCLvalues$info_pz_MCL0208 -> info_tmp)
    req(input$tissueMCL -> tissue)
    req(input$pharma_varsMCL -> variable)
    
    if(input$ClusterCheckMCL){
      if(tissue=="BM"){
        palette<-c("#F79071","#994627")
      }
      else if(tissue=="PB"){
        palette<-c("#3B3571", "#62C75F")
      }
      info_tmp[,paste0("Cluster_",tissue)]<- ifelse( info_tmp[[paste0("Cluster_",tissue)]] %in% c("A","B"),"A/B","C/D")
    }
    else{
      if(tissue=="BM"){
        palette<-c("#F1BB7B", "#FD6467", "#5B1A18", "#D67236")
      }
      else if(tissue=="PB"){
        palette<-c("#440154", "#31688E", "#35B779", "#8FD744")
      }
    }
    
    test<-test_indipendence(info_tmp,variable,paste0("Cluster_",tissue),palette)
    test$plot  & (
      theme(
        legend.position = "bottom",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent')
      )+axistheme
    )
  })
  
  #### END: MCL208 exploration ####
  
  ### YoungMCL classification ####
  
  listV = reactiveValues(CONNECTORList.FCM = list(),
                         classificationFCM = NULL,
                         DataFCMtrunc = NULL,
                         SurvDataINIT=NULL,
                         SurvDataOUT=NULL,
                         ClassificationINIT = NULL,
                         ClassificationOUT = NULL,
                         DisTimeBX = NULL,
                         pl = NULL,
                         LPanalysis = NULL)
  
  ColorCluster <- readRDS( system.file("Data","ColorCluster.RDs", package = "ShinyMCL") )
  Steps <- c("BASELINE","RCHOP","ARAC","ASCT",paste0("M",seq(6,36,6)))
  
  observeEvent(input$tissueBox,{  
    #show_modal_spinner()# show the modal window
    output$classifiedCurve <- renderPlot({ggplot()})
    listV$DataFCMtrunc = NULL
    listV$LPanalysis = NULL
    listV$pl = ggplot()
    
    listV$CONNECTORList.FCM = WholeData$MCL0208[[input$tissueBox]]$CONNECTORList.FCM
    listV$classificationFCM = WholeData$YoungMCL[[input$tissueBox]]$Dataset
    
    # removing the lp analysis plots
    if(!is.null(listV$LPanalysis)){
      for( i in seq_along(listV$LPanalysis)){
        plotname <- paste("plot", i, sep="")
        removeUI(plotname)
      }
    }
    
    ## plot the boxplot regular time ####
    listV$CONNECTORList.FCM$FCM$cluster$ClustCurve$Cluster = listV$CONNECTORList.FCM$FCM$cluster$cluster.names[listV$CONNECTORList.FCM$FCM$cluster$ClustCurve$Cluster]
    df = merge(listV$CONNECTORList.FCM$FCM$cluster$ClustCurve,listV$CONNECTORList.FCM$CONNECTORList$LabCurv)
    df = merge(df,listV$CONNECTORList.FCM$CONNECTORList$Dataset)
    df$Step = factor(df$Step, levels = Steps)
    
    DisTimeBX = df %>%
      mutate(y = 1,Step = factor(Step,levels = Steps[-1])) %>%
      na.omit() %>%
      ggplot()+
      geom_boxplot(aes( x= Time/365, col = Step, fill = Step),alpha = 0.5,position="identity" ) +
      theme_bw() +
      labs(x = "")+
      theme(legend.position = "none",
            axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            plot.margin = unit(c(0, 0,0,0), "cm"))+
      geom_label( data = df %>% group_by(Step) %>% summarise(Time = median(Time / 365)),
                  aes(x = Time, label = Step, y=-0.4, col = Step) ,fill = "white", fontface = "bold",
                  #position = position_dodge(width = 0.75),  # Adjust the width as needed
                  #vjust = -0.5,  # Adjust the vertical justification as needed
                  size = 5  # Adjust the size of the text as needed
      ) 
    listV$DisTimeBX = DisTimeBX
    
    ####
    Info = WholeData$YoungMCL$Info
    
    updateSliderInput("Cut",
                      session = session,
                      min = min(listV$classificationFCM$Time),
                      max = max(listV$classificationFCM$Time),
                      value = max(listV$classificationFCM$Time) )
    
    listV$ClassificationINIT  = SurvYOUNGERData[[input$tissueBox]]
    #remove_modal_spinner() # remove it when done
  })
  
  observe({
    classification= req(listV$ClassificationINIT)
    Dataset = WholeData$YoungMCL[[input$tissueBox]]$Dataset
    ## change color depedning on input$youngMCLselectColor
    DatasetClusters= merge(Dataset,classification$Classification$ClassMatrix_entropy,by = "ID")
    
    colorVar = c("None", "ttpevent","rnd1","INIERG")
    names(colorVar) = c("None", "TTP","Random","INIERG")
    selectColor=unname(colorVar[input$youngMCLselectColor])
    
    output$classYoungPlot <- renderPlot({
      
      pl = ggplot(DatasetClusters,aes(x = Time, y = Observation, group=ID))+
        facet_wrap(~Cluster,ncol = 1) +
        labs(x = "Days", y = "") +
        scale_y_continuous(limits=c(-1, 8) ,
                           breaks = seq(0,8,2),
                           labels = c("NEG","POS",TeX("$10^{-3}$"),TeX("$10^{-2}$"),TeX("$10^{-1}$")) ) +
        theme_bw() +
        theme(legend.position = "top",
              plot.margin = unit(c(0, 0,0,0), "cm"))
      
      if("Unclassified" %in% unique(DatasetClusters$Cluster)  )
      {
        pl = pl +geom_rect(
          aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
          data.frame(Cluster = 'Unclassified', Time = 0,Observation = 0, ID= -1),
          fill = 'grey', alpha = 0.5
        )
      }
      
      if(selectColor == "None")
        pl = pl+ geom_line()
      else {
        DatasetClusters = DatasetClusters %>% mutate_at(.vars= selectColor, .funs = as.factor )
        pl = pl +
          geom_exec(geom_line, data =  DatasetClusters ,
                    x = "Time",y = "Observation", group = "ID", color = selectColor )+
          labs(col = input$youngMCLselectColor )
        
        if(selectColor %in% c("ttpevent","rnd1") )
          pl = pl + scale_color_manual(values = c("0" = "blue","1"="red"))
      }
      
      
      pl / listV$DisTimeBX + plot_layout(heights = c(2,0.5))
    })
    
    ggsurvData = classification$ggsurvData
    fit = classification$fit
    Info_tmp= classification$Info_tmp
    
    colors = ColorCluster[grep(x = names(ColorCluster),pattern =  input$tissueBox)]
    names(colors) = gsub(x = names(colors),replacement = "Cluster=",pattern =  paste0(input$tissueBox,"_"))
    list2env(list(Info_tmp = Info_tmp), envir = .GlobalEnv)
    ggsurv=ggsurvplot(fit = fit, data = Info_tmp,
                      xlab = "Year",  ylab = "TTP",
                      size = 1, pval = TRUE, risk.table = TRUE,conf.int = T,
                      risk.table.col="strata",ggtheme = theme_bw(),surv.median.line = "hv"  )
    pval = surv_pvalue(fit = fit, data = Info_tmp, get_coord = TRUE)
    
    ggsurvplot = ggplot() +
      geom_ribbon(
        data = ggsurv$plot$data,
        aes(
          x = time, 
          ymin = lower, 
          ymax = upper, 
          fill = strata
        ),
        alpha = 0.2  # Adjust transparency as needed
      )+
      geom_step( data = ggsurv$plot$data, aes(x = time, y = surv,  color = strata) ) +
      geom_point(data = ggsurv$plot$data[ggsurv$plot$data$n.censor > 0, , drop = FALSE],
                 aes(x = time, y = surv,  color = strata), size = 4, shape = "+")+
      theme_bw()+ labs(x = "Year", y = "TTP", fill = "Cluster", color = "Cluster")+
      scale_color_manual(values = colors)+scale_fill_manual(values = colors)+
      scale_x_continuous(breaks = c(0,3,6,9,12,15)*365.25, labels = c(0,3,6,9,12,15))+
      annotate("text", x = pval$pval.x, y = pval$pval.y,
               label = pval$pval.txt, size = 5, hjust = 0)
    
    surv_median_values = surv_median(fit = fit)
    if(! all(is.na(surv_median_values$median)) ){
      ggsurvplot = ggsurvplot +
        geom_segment(
          data = surv_median_values,
          aes(
            x = median, xend = median, 
            y = 0, yend = 0.5, 
            color = strata
          ),
          linetype = "dashed"
        )+
        geom_segment(
          data = surv_median_values,
          aes(
            x = 0, xend = median, 
            y = 0.5, yend = 0.5, 
            color = strata
          ),
          linetype = "dashed"
        )  +
        annotate("text", x = surv_median_values$median, y = 0,
                 label = round(surv_median_values$median/365,digits = 3), size = 3, hjust = 0)
    }
    
    output$survYoungPlot <- renderPlot({ggsurvplot})
    
    if(selectColor != "None")
    {
      counts =DatasetClusters[,c("ID","Cluster",selectColor)] %>% distinct()
      counts= table(counts$Cluster, counts[,selectColor] )
      row_names <- rownames(counts)
      col_names <- colnames(counts)
      
      # Initialize an empty vector to hold the results
      output_list <- c()
      
      # Loop through each cell in the table
      for (i in seq_along(row_names)) {
        for (j in seq_along(col_names)) {
          # Create a string for each combination of row and column
          cell <- paste(row_names[i], col_names[j], counts[i, j], sep=": ")
          # Add the string to the output list
          output_list <- c(output_list, cell)
        }
      }
      
      # Combine all cell strings into one final string
      outputCount <- paste(output_list, collapse=",\n ")
    }else{outputCount=""}
    
    YmclSummary <- paste(
      "### Summary of Young MCL data. ",
      "Number of curves: ",length(classification$Classification$ClassMatrix_entropy$ID),
      "Classification results:",
      paste(names(table(classification$Classification$ClassMatrix_entropy$Cluster)),": ",
            table(classification$Classification$ClassMatrix_entropy$Cluster),collapse = ",\n"),
      ifelse(selectColor != "None",
             paste("### Stratification\n", outputCount, "\n"), no = ""
      ),
      # Add more summary information if needed
      sep = "\n"
    )
    
    output$YmclSummary <- renderText({
      YmclSummary
    })
    
  })
  
  ##### Line Plot ####
  observe({
    input$youngMCLselectColor
    
    colorVar = c("None", "ttpevent","rnd1","INIERG")
    names(colorVar) = c("None", "TTP","Random","INIERG")
    selectColor=unname(colorVar[input$youngMCLselectColor])
    
    Data = WholeData$YoungMCL[[input$tissueBox]]$Dataset
    vline = input$Cut
    pl = Data %>% ggplot() +
      geom_vline(xintercept = vline,
                 col = "red",
                 linetype = "dashed")+
      scale_y_continuous(limits=c(-1, 8) ,
                         breaks = seq(0,8,2),
                         labels = c("NEG","POS",TeX("$10^{-3}$"),TeX("$10^{-2}$"),TeX("$10^{-1}$")) ) +
      theme_bw() +
      theme(legend.position = "top",
            plot.margin = unit(c(0, 0,0,0), "cm"))
    
    DataTrunc = listV$DataFCMtrunc
    
    isolate({
      output$linePlot <- renderPlot({
        if(!is.null(Data)){
          if(selectColor!="None"){
            Data = Data %>% mutate_at(.vars=selectColor, .funs = as.factor )
            
            pl = pl + geom_exec(geom_line, data =  Data ,
                                x = "Time",y = "Observation", group = "ID", color = selectColor )+
              labs(col = input$youngMCLselectColor )
            
            if(selectColor %in% c("ttpevent","rnd1") )
              pl = pl + scale_color_manual(values = c("0" = "blue","1"="red"))
            
          }else
            pl = pl + geom_line(data = Data, aes(x = Time,y = Observation, group = ID))
          
          pl = pl / listV$DisTimeBX + plot_layout(heights = c(2,0.5))
          return(pl)
        }
      })
      
      output$lineTruncPlot <- renderPlot({
        if(!is.null(listV$DataFCMtrunc)){
          
          pl = listV$DataFCMtrunc %>%
            as.data.frame() %>%
            ggplot() +
            theme_bw()+
            theme(legend.position = "none",
                  plot.margin = unit(c(0, 0,0,0), "cm"))+
            xlim(min(Data$Time),
                 max(Data$Time))
          
          if(selectColor!="None"){
            listV$DataFCMtrunc = listV$DataFCMtrunc %>% 
              mutate_at(.vars=selectColor,.funs = as.factor )
            
            listV$DataFCMtrunc[[selectColor]] = factor(x = as.character(listV$DataFCMtrunc[[selectColor]]), levels =  levels(as.factor(Data[[selectColor]])) )
            
            pl = pl + geom_exec(geom_line,
                                data =  listV$DataFCMtrunc %>%
                                  as.data.frame() ,
                                x = "Time",y = "Observation", group = "ID", color =selectColor )+
              labs(col = input$youngMCLselectColor )
            
            if(selectColor %in% c("ttpevent","rnd1") )
              pl = pl + scale_color_manual(values = c("0" = "blue","1"="red"))
            
          } else
            pl = pl + geom_line(aes(x = Time,y = Observation, group = ID)) 
          
          return(pl)
        }
      })
      
    })
  })
  
  observeEvent(listV$pl,{
    output$classifiedCurve <- renderPlot(
      listV$pl 
    )
  })
  
  #### LP analysis 
  
  observeEvent(input$Cut,{
    if(!is.null(listV$classificationFCM) && input$Cut > min(listV$classificationFCM$Time) ){
      DataFCMtrunc  =  listV$classificationFCM %>%
        filter(Time <= input$Cut )  %>% 
        group_by(ID) %>%
        mutate(L = length(ID)) 
      
      L =  length(unique(DataFCMtrunc$ID[DataFCMtrunc$L< 3]))
      
      DataFCMtrunc  =  DataFCMtrunc %>%
        filter(L>2) %>%
        select(-L) %>%
        mutate( ttpDiff = ttp*365.25 - input$Cut )  %>% 
        filter( ttpDiff >= 0 )
      
      if(nrow(DataFCMtrunc)==0 ){
        shinyalert("Warning", "No curves are selected", type = "error")
      }else{
        listV$DataFCMtrunc = DataFCMtrunc
        
        # Assuming 'data' is your dataset and 'curves' is the number of curves
        truncation_summary <- paste(
          "### Summary of Truncation: ",
          "Number of Curves:", length(unique(listV$classificationFCM$ID)),
          "Truncation Value:", input$Cut,
          "Number of curves <= 2 points:", L,
          "Number of curves > 2 points with TTP event before LP",
          length(unique(listV$classificationFCM$ID)) - length(unique(listV$DataFCMtrunc$ID)) - L,
          "Final number of truncated curves","for the analysis:", length(unique(listV$DataFCMtrunc$ID)),
          # Add more summary information if needed
          sep = "\n"
        )
        
        output$truncationSummary <- renderText({
          truncation_summary
        })
      }
    }
  })
  
  observeEvent(input$GoClass,{
    show_modal_spinner()# show the modal window
    
    if(dim(listV$DataFCMtrunc)[1] > 0 && input$Cut <=  max(listV$classificationFCM$Time) ){
      listV$ClassificationOUT= ClassificationNewCurves(listV$DataFCMtrunc,
                                                       listV$CONNECTORList.FCM, Cores = 2)
      
      
      ClassMatrix_entropy =  listV$ClassificationOUT$ClassMatrix_entropy %>%
        filter(Cluster != "Unclassified")
      
      Info_tmp= merge(listV$DataFCMtrunc %>% dplyr::select(ID,INIERG,rnd1, ttpevent,ttp, ttpDiff) %>% distinct(), 
                      as.data.frame(ClassMatrix_entropy), by = "ID")
      
      SmallNumberIDs = Info_tmp %>% group_by(Cluster) %>% dplyr::summarise(N = n()) %>% filter(N <3)
      if(dim(SmallNumberIDs)[1] > 0 )
      {
        shinyalert("Warning", paste0("The cluster(s) ", paste(SmallNumberIDs$Cluster,collapse = ", "), " can not be used in the LP analysis given the number of curves is less than 3"), "warning", 5000)
        Info_tmp = Info_tmp %>% filter(! Cluster %in% SmallNumberIDs$Cluster )
      }
      
      ############################ SURVIVAL ANALYSIS
      
      fit <- survfit(Surv(Info_tmp$ttpDiff, as.numeric(Info_tmp$ttpevent) )~ Cluster, data = Info_tmp)
      list2env(list(Info_tmp = Info_tmp, fit = fit), envir = .GlobalEnv)
      
      ggsurv = NULL
      
      tryCatch({
        ggsurv=ggsurvplot(fit = fit, data = Info_tmp,
                          xlab = "Year",  ylab = "TTP",
                          size = 1, pval = TRUE, risk.table = TRUE,conf.int = T,
                          risk.table.col="strata",ggtheme = theme_bw(),surv.median.line = "hv"  )
      }, error = function(e) {
        # Display the error message in a Shiny alert if an error occurs
        shinyalert(
          title = "Error",
          text = paste("There was an issue generating the survival plot:", e$message),
          type = "error"
        )
      })
      
      listV$ClassificationOUT$ggsurv = ggsurv
      listV$ClassificationOUT$Info_tmp = Info_tmp
      listV$ClassificationOUT$fit = fit
    }
    remove_modal_spinner() # remove it when done
    ### change tab panel 
    updateTabsetPanel(
      inputId="YoungPanels",
      selected = "panel_Classification"
    )
    
  })
  
  observe({
    input$GoClass
    input$youngMCLselectColor
    
    isolate({
      input$Cut -> CutTime
      #cluster = input$clusterBox
      if(!is.null(listV$ClassificationOUT) && !is.null(listV$ClassificationOUT$ggsurv) ){
        
        listV$ClassificationOUT$ggsurv -> ggsurv
        listV$ClassificationOUT$Info_tmp -> Info_tmp
        listV$ClassificationOUT$fit -> fit
        
        colors = ColorCluster[grep(x = names(ColorCluster),pattern =  input$tissueBox)]
        names(colors) = gsub(x = names(colors),replacement = "Cluster=",pattern =  paste0(input$tissueBox,"_"))
        
        list2env(list(Info_tmp = Info_tmp), envir = .GlobalEnv)
        pval = surv_pvalue(fit = fit, data = Info_tmp, get_coord = TRUE)
        
        # ggsurvplot = ggplot()  +
        #   geom_ribbon(
        #     data = ggsurv$plot$data,
        #     aes(
        #       x = CutTime+time, 
        #       ymin = lower, 
        #       ymax = upper, 
        #       fill = strata
        #     ),
        #     alpha = 0.2  
        #   ) 
        
        ggsurvplot =  ggplot(data = ggsurv$plot$data, aes( x = CutTime+time))  +
          #ggsurv$plot$layers[[3]] + 
          geom_step( data = ggsurv$plot$data, aes(x = CutTime+ time, y = surv,  color = strata),linewidth =1 ) +
          geom_point(data = ggsurv$plot$data[ggsurv$plot$data$n.censor > 0, , drop = FALSE],
                     aes(x = CutTime+ time, y = surv,  color = strata), size = 6, shape = "+")+
          theme_bw()+ labs(x = "Year", y = "TTP", color = "Cluster",fill = "Cluster")+
          scale_color_manual(values = colors)+
          scale_fill_manual(values = colors)+
          scale_x_continuous(breaks = c(0,3,6,9,12,15)*365.25, labels = c(0,3,6,9,12,15))+
          annotate("text", x = pval$pval.x, y = pval$pval.y,
                   label = pval$pval.txt, size = 5, hjust = 0)
        
        surv_median_values = surv_median(fit = fit)
        
        if(! all(is.na(surv_median_values$median)) ){
          ggsurvplot = ggsurvplot +
            geom_segment(
              data = surv_median_values,
              aes(
                x = CutTime+median, xend = CutTime+median, 
                y = 0, yend = 0.5, 
                color = strata
              ),
              linetype = "dashed"
            )+
            geom_segment(
              data = surv_median_values,
              aes(
                x = 0, xend = CutTime+median, 
                y = 0.5, yend = 0.5, 
                color = strata
              ),
              linetype = "dashed"
            )  +
            annotate("text", x = CutTime+surv_median_values$median, y = 0,
                     label = round(CutTime+surv_median_values$median/365,digits = 3), size = 3, hjust = 0)
        }
        
        Classification_new = listV$ClassificationOUT$ClassMatrix_entropy %>% select(ID,Cluster)
        meancurves = as.data.frame(listV$CONNECTORList.FCM$FCM$prediction$meancurves)
        colnames(meancurves) = listV$CONNECTORList.FCM$FCM$cluster$cluster.names
        meancurves$Time = listV$CONNECTORList.FCM$CONNECTORList$TimeGrid
        meancurves = meancurves %>% tidyr::gather(-Time, key = "Cluster", value = "Observation")
        
        DatasetClusters = merge(listV$DataFCMtrunc, Classification_new)
        
        pl1 = DatasetClusters %>%
          ggplot()+
          facet_grid(~Cluster) 
        
        if("Unclassified" %in% unique(DatasetClusters$Cluster)  )
        {
          pl1 = pl1 +geom_rect(
            aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
            data.frame(Cluster = 'Unclassified', Time = 0,Observation = 0, ID= -1),
            fill = 'grey', alpha = 0.5
          )
        }
        
        colorVar = c("None", "ttpevent","rnd1","INIERG")
        names(colorVar) = c("None", "TTP","Random","INIERG")
        selectColor=unname(colorVar[input$youngMCLselectColor])
        
        if(selectColor != "None"){
          DatasetClusters = DatasetClusters %>% mutate_at(.vars= selectColor, .funs = as.factor )
          pl1 = pl1 +
            geom_exec(geom_line, data =  DatasetClusters ,
                      x = "Time",y = "Observation", group = "ID", color = selectColor )+
            labs(col = input$youngMCLselectColor )
          
          if(selectColor %in% c("ttpevent","rnd1") )
            pl1 = pl1 + scale_color_manual(values = c("0" = "blue","1"="red", "Mean Curves" = "green"))
        }else{
          pl1 = pl1 + geom_line(aes(x=Time, y = Observation, group = ID))
        }
        
        pl1 = pl1 +
          geom_line(data = meancurves, aes(x=Time, y = Observation, col = "Mean Curves"),linewidth =1, linetype = "dashed")+
          labs(title = paste(input$tissueBox, " classification") )+
          theme_bw()
        listV$ggsurvLast = ggsurvplot
        listV$pl = pl1/ggsurvplot
        
      }else{
        listV$pl = ggplot()
      }
    })
  })
  
  observeEvent(list(input$saveLP,input$saveLP_fromClass),{
    isolate({
      if(!is.null(listV$ggsurvLast)){
        timeLP = input$Cut
        
        if(is.null(listV$LPanalysis)){
          listV$LPanalysis = list( pl = listV$ggsurvLast + labs(title = paste0("LP at ", timeLP) ) + theme(legend.position = "top")) 
        }else{
          listV$LPanalysis[[length(listV$LPanalysis)+1]] =  (listV$ggsurvLast + labs(title = paste0("LP at ", timeLP, " days"))  + theme(legend.position = "none"))
        }
        
        names(listV$LPanalysis) = paste0("plot",1:length(listV$LPanalysis))
        
        ### change tab panel 
        updateTabsetPanel(
          inputId="YoungPanels",
          selected = "panel_Summary"
        )
      }
    })
  })
  
  output$landmarkPoints <- renderUI({
    if (!is.null(listV$LPanalysis)) {
      for( i in seq_along(listV$LPanalysis) ){
        plotname <- paste("plot", i, sep="")
        removeUI(plotname)
      }
      
      plot_output_list <- lapply(seq_along(listV$LPanalysis), function(i) {
        plotname <- paste("plot", i, sep="")
        plotOutput(outputId = plotname)
      })
      
      # Convert the list to a tagList - this is necessary for the list of items
      # to display properly.
      do.call(tagList, plot_output_list)
      #div(plotOutput(outputId = paste("plot",length(listV$LPanalysis), sep="")))
    }
  })
  
  observeEvent(list(input$saveLP_fromClass,input$saveLP),{
    
    if (!is.null(listV$LPanalysis)) {
      lapply(seq_along(listV$LPanalysis), function(i) {
        p = listV$LPanalysis[[paste("plot", i, sep="")]]
        output[[paste("plot",i, sep="")]] <- renderPlot({ p   })
      })
    }
  })
  #### END: YoungMCL classification ####
  
  #### USER import #### 
  ListValuesUsers <- reactiveValues(
    excel_df = NULL,
    txt_df = NULL,
    merged_df = NULL,
    ClassificationINIT = NULL,
    DataFCMtrunc = NULL,
    ggsurvLast = NULL,
    LPanalysis = NULL,
    pl = NULL
  )
  
  # Process and render Excel file
  excel_data <- reactive({
    req(input$file1)
    read_excel(input$file1$datapath)
  })
  
  output$column_selectors <- renderUI({
    df <- excel_data()
    colnames <- names(df)
    tagList(
      selectInput("time_col", "Select Time Column:", choices = colnames),
      selectInput("obs_col", "Select Observation Column:", choices = colnames),
      selectInput("id_col", "Select ID Column:", choices = colnames)
    )
  })
  
  output$contents1 <- renderTable({
    req(ListValuesUsers$excel_df)
    head(ListValuesUsers$excel_df)
  })
  
  observe({
    req(input$time_col, input$obs_col, input$id_col)
    df <- excel_data()
    selected_df <- df %>%
      select(Time = !!sym(input$time_col), Observation = !!sym(input$obs_col), ID = !!sym(input$id_col))
    ListValuesUsers$excel_df <- selected_df
  })
  
  output$info1 <- renderPrint({
    req(ListValuesUsers$excel_df)
    df <- ListValuesUsers$excel_df
    num_ids <- n_distinct(df[[input$id_col]])
    num_obs <- nrow(df)
    cat(paste("Number of IDs:", num_ids, "\nNumber of Observations per ID:", num_obs))
  })
  
  # Process and render TXT file
  txt_data <- reactive({
    req(input$file2)
    read.table(input$file2$datapath, header = input$header2, sep = input$sep2, quote = input$quote2)
  })
  
  output$txt_id_selector <- renderUI({
    df <- txt_data()
    colnames <- names(df)
    selectInput("txt_id_col", "Select TXT ID Column:", choices = colnames)
  })
  
  observe({
    req(input$txt_id_col)
    df <- txt_data()
    
    tryCatch({
      selected_df <- df %>% rename(ID = !!sym(input$txt_id_col))
      ListValuesUsers$txt_df <- selected_df
    }, error = function(e) NULL)
    
  })
  
  output$contents2 <- renderTable({
    req(ListValuesUsers$txt_df)
    head(ListValuesUsers$txt_df)
  })
  
  output$info2 <- renderPrint({
    req(ListValuesUsers$txt_df)
    df <- ListValuesUsers$txt_df
    num_features <- ncol(df) - 1
    cat(paste("Number of Features:", num_features))
  })
  
  observe({
    req(ListValuesUsers$excel_df,ListValuesUsers$txt_df)
    excel_ids <- ListValuesUsers$excel_df$ID
    txt_ids <- ListValuesUsers$txt_df$ID
    
    if(length(colnames(ListValuesUsers$txt_df))<3){
      ListValuesUsers$merged_df <- NULL
      shinyalert("Error", "The features dataset in the TXT must have more than 2 columns, for the survival analysis!", type = "error")
      return()
    }
    
    if (all(excel_ids %in% txt_ids)) {
      merged_df <- merge(ListValuesUsers$excel_df, ListValuesUsers$txt_df, by.x = "ID", by.y = "ID", all.x = TRUE)
      ListValuesUsers$merged_df <- merged_df
    } else {
      ListValuesUsers$merged_df <- NULL
      shinyalert("Error", "No matching ID found between the datasets.", type = "error")
    }
  })
  
  observeEvent(ListValuesUsers$merged_df, {
    req(ListValuesUsers$merged_df)
    df <- ListValuesUsers$txt_df
    colnames <- names(df)
    updateSelectInput("color_col", session = session, choices = c("None",colnames), selected = "None" )
    updateSelectInput("selectColor_user", session = session, choices = c("None",colnames), selected = "None")
    updateSelectInput("featureEventSurv_user", session = session, choices = c("",colnames[-1]), selected = "" )
    updateSelectInput("featureTimeSurv_user", session = session, choices = c("",colnames[-1]), selected = "" )
    
  })
  
  # Generate Line Plot
  output$linePlotUser <- renderPlot({
    req(ListValuesUsers$merged_df)
    df <- ListValuesUsers$merged_df
    
    if(input$color_col == "None"){
      pl=  ggplot(df,
                  aes(x = Time, y = Observation, group = ID )
      ) +
        geom_line() +
        labs(title = "", x = "Time", y = "Observation")+
        theme_bw()
    }else{
      pl=  ggplot(df,
                  aes(x = Time, y = Observation, group = ID, color = as.factor(!!sym(input$color_col)) )
      ) +
        geom_line() +
        labs(title = "", x = "Time", y = "Observation")+
        theme_bw()
    }
    
    pl
    
  })
  
  observeEvent(input$GoClass_Tabuser,{
    updateTabsetPanel(session, "SideTabs",
                      selected = "user_classification")
  })
  
  
  ### Classification section
  
  observeEvent(input$tissueBox_user,{
    ListValuesUsers$CONNECTORList.FCM =  WholeData$MCL0208[[input$tissueBox_user]]$CONNECTORList.FCM
  })
  
  
  observe({  
    input$tissueBox_user
    req(ListValuesUsers$merged_df)
    
    isolate({
      show_modal_spinner()# show the modal window
      output$classifiedCurve_user <- renderPlot({ggplot()})
      ListValuesUsers$DataFCMtrunc = NULL
      ListValuesUsers$LPanalysis = NULL
      ListValuesUsers$pl = ggplot()
      ListValuesUsers$ClassificationINIT = NULL
      
      # removing the lp analysis plots
      if(!is.null(ListValuesUsers$LPanalysis)){
        for( i in seq_along(ListValuesUsers$LPanalysis)){
          plotname <- paste("plot", i, sep="")
          removeUI(plotname)
        }
      }
      
      ## plot the boxplot regular time ####
      ListValuesUsers$CONNECTORList.FCM$FCM$cluster$ClustCurve$Cluster = ListValuesUsers$CONNECTORList.FCM$FCM$cluster$cluster.names[ListValuesUsers$CONNECTORList.FCM$FCM$cluster$ClustCurve$Cluster]
      df = merge(ListValuesUsers$CONNECTORList.FCM$FCM$cluster$ClustCurve,ListValuesUsers$CONNECTORList.FCM$CONNECTORList$LabCurv)
      df = merge(df,ListValuesUsers$CONNECTORList.FCM$CONNECTORList$Dataset)
      df$Step = factor(df$Step, levels = Steps)
      
      DisTimeBX = df %>%
        mutate(y = 1,Step = factor(Step,levels = Steps[-1])) %>%
        na.omit() %>%
        ggplot()+
        geom_boxplot(aes( x= Time/365, col = Step, fill = Step),alpha = 0.5,position="identity" ) +
        theme_bw() +
        labs(x = "", y = "")+
        theme(legend.position = "none",
              axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.title.y=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              plot.margin = unit(c(0, 0,0,0), "cm"))+
        geom_label( data = df %>% group_by(Step) %>% summarise(Time = median(Time / 365)),
                    aes(x = Time, label = Step, y=-0.4, col = Step) ,fill = "white", fontface = "bold",
                    #position = position_dodge(width = 0.75),  # Adjust the width as needed
                    #vjust = -0.5,  # Adjust the vertical justification as needed
                    size = 5  # Adjust the size of the text as needed
        )  
      ListValuesUsers$DisTimeBX = DisTimeBX
      
      ####
      
      updateSliderInput("Cut_user",
                        session = session,
                        min = min(ListValuesUsers$merged_df$Time),
                        max = max(ListValuesUsers$merged_df$Time),
                        value = max(ListValuesUsers$merged_df$Time) )
      
      
      remove_modal_spinner() # remove it when done
      
    })
  })
  
  observe({
    req(ListValuesUsers$merged_df)
    
    EventSU = req(input$featureEventSurv_user)
    TimeSU = req(input$featureTimeSurv_user)
    
    if(EventSU == "" || TimeSU == "" ) return()
    
    if(!is.numeric(ListValuesUsers$merged_df[[EventSU]]) || !is.numeric(ListValuesUsers$merged_df[[TimeSU]])){ 
      shinyalert("Error", "The variables selected for the survival analysis are not numeric.", type = "error")
      return()
    }
    
    if( !all(unique(ListValuesUsers$merged_df[[EventSU]]) %in% c(0,1)) ){ 
      shinyalert("Error", "The event must be coded as 0 and 1.", type = "error")
      return()
    }
    
    show_modal_spinner()# show the modal window
    
    if(is.null(ListValuesUsers$ClassificationINIT)){
      Classification= ClassificationNewCurves(ListValuesUsers$merged_df[,c("ID","Observation","Time")], 
                                              ListValuesUsers$CONNECTORList.FCM,
                                              Cores = 2)
      
      ClassMatrix_entropy =  Classification$ClassMatrix_entropy %>%
        filter(Cluster != "Unclassified")
      
      DatasetClusters= merge(ListValuesUsers$merged_df, 
                             as.data.frame(ClassMatrix_entropy),
                             by = "ID")
      
      Info_tmp = DatasetClusters %>% dplyr::select(-Time,-Observation) %>% distinct()
      
      ############################ SURVIVAL ANALYSIS
      SmallNumberIDs = Info_tmp %>% group_by(Cluster) %>% dplyr::summarise(N = n()) %>% filter(N <3)
      if(dim(SmallNumberIDs)[1] > 0 )
      {
        shinyalert("Warning", paste0("The cluster(s) ", paste(SmallNumberIDs$Cluster,collapse = ", "), " can not be used in the LP analysis given the number of curves is less than 3"), "warning", 5000)
        Info_tmp = Info_tmp %>% filter(! Cluster %in% SmallNumberIDs$Cluster )
      }
      fit <- eval(parse(text = paste0("survfit(Surv(Info_tmp$",TimeSU,",as.numeric(Info_tmp$",EventSU,")) ~ Cluster, data = Info_tmp)")))
      
      list2env(list(Info_tmp = Info_tmp), envir = .GlobalEnv)
      ggsurv=ggsurvplot(fit = fit, data = Info_tmp,
                        xlab = "Time",  ylab = EventSU,
                        size = 1, pval = TRUE, risk.table = TRUE,conf.int = T,
                        risk.table.col="strata",ggtheme = theme_bw(),surv.median.line = "hv"  )
      
      ListValuesUsers$ClassificationINIT  = list(ggsurvData = ggsurv$data.survplot,
                                                 Classification = ClassMatrix_entropy)
      
      output$survUserPlot <- renderPlot({ggsurv})
      
      
    }else{
      DatasetClusters= merge(ListValuesUsers$merged_df, 
                             as.data.frame(ListValuesUsers$ClassificationINIT$Classification),
                             by = "ID")
    }
    
    ### Plot curves
    output$classUserPlot <- renderPlot({
      pl = ggplot(DatasetClusters,aes(x = Time, y = Observation, group=ID))+
        facet_wrap(~Cluster,ncol = 1) +
        labs(x = "Days", y = "") +
        theme_bw() +
        theme(legend.position = "top",
              plot.margin = unit(c(0, 0,0,0), "cm"))
      
      if("Unclassified" %in% unique(DatasetClusters$Cluster)  )
      {
        pl = pl +geom_rect(
          aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
          data.frame(Cluster = 'Unclassified', Time = 0,Observation = 0, ID= -1),
          fill = 'grey', alpha = 0.5
        )
      }
      
      if(input$selectColor_user != "None"){
        DatasetClusters = DatasetClusters %>% mutate_at(.vars=input$selectColor_user,.funs = as.factor )
        pl = pl + geom_exec(geom_line, data =  DatasetClusters ,
                            x = "Time",y = "Observation", group = "ID", color = input$selectColor_user )
      }else{
        pl = pl + geom_line()
      }
      
      pl 
    })
    
    remove_modal_spinner()
    
    if(input$selectColor_user != "None")
    {
      counts =DatasetClusters[,c("ID","Cluster",input$selectColor_user)] %>% distinct()
      counts= table(counts$Cluster, counts[,input$selectColor_user] )
      row_names <- rownames(counts)
      col_names <- colnames(counts)
      
      # Initialize an empty vector to hold the results
      output_list <- c()
      
      # Loop through each cell in the table
      for (i in seq_along(row_names)) {
        for (j in seq_along(col_names)) {
          # Create a string for each combination of row and column
          cell <- paste(row_names[i], col_names[j], counts[i, j], sep=": ")
          # Add the string to the output list
          output_list <- c(output_list, cell)
        }
      }
      
      # Combine all cell strings into one final string
      outputCount <- paste(output_list, collapse=",\n ")
    }else{outputCount=""}
    
    UserSummary <- paste(
      "### Summary of User Data. ",
      "Number of curves: ",length( ListValuesUsers$ClassificationINIT$Classification$ID),
      "Classification results:",
      paste(names(table( ListValuesUsers$ClassificationINIT$Classification$Cluster)),": ",
            table( ListValuesUsers$ClassificationINIT$Classification$Cluster),collapse = ",\n"),
      ifelse(input$selectColor_user != "None",
             paste("### Stratification\n", outputCount, "\n"), no = ""
      ),
      # Add more summary information if needed
      sep = "\n"
    )
    
    output$UserSummary <- renderText({
      UserSummary
    })
    
  })
  
  ##### Line PLot ####
  observe({
    req(ListValuesUsers$merged_df)
    
    input$selectColor_user
    
    Data = ListValuesUsers$merged_df
    vline = input$Cut_user
    
    pl = Data %>% ggplot() +
      theme_bw()+
      geom_vline(xintercept = vline,
                 col = "red",
                 linetype = "dashed")
    
    DataTrunc = ListValuesUsers$DataFCMtrunc
    
    isolate({
      output$linePlot_user <- renderPlot({
        if(!is.null(Data)){
          if(input$selectColor_user != "None"){
            Data = Data %>% mutate_at(.vars=input$selectColor_user,.funs = as.factor )
            
            pl = pl + geom_exec(geom_line, data =  Data ,
                                x = "Time",y = "Observation", group = "ID", color = input$selectColor_user )
          }else
            pl = pl + geom_line(data = Data, aes(x = Time,y = Observation, group = ID))
          
          pl = pl / ListValuesUsers$DisTimeBX + plot_layout(heights = c(2,0.5))
          return(pl)
        }
      })
      
      output$lineTruncPlot_user <- renderPlot({
        if(!is.null(ListValuesUsers$DataFCMtrunc)){
          
          pl = ListValuesUsers$DataFCMtrunc %>%
            as.data.frame() %>%
            ggplot() +
            theme_bw()+
            xlim(min(Data$Time),
                 max(Data$Time))
          
          if(input$selectColor_user!="None"){
            ListValuesUsers$DataFCMtrunc = ListValuesUsers$DataFCMtrunc %>% 
              mutate_at(.vars=input$selectColor_user,.funs = as.factor )
            
            pl = pl + geom_exec(geom_line,
                                data =  ListValuesUsers$DataFCMtrunc %>%
                                  as.data.frame() ,
                                x = "Time",y = "Observation", group = "ID", color =input$selectColor_user )
          } else{
            pl = pl + geom_line(aes(x = Time,y = Observation, group = ID))
          }
          
          return(pl)
        }
      })
      
    })
  })
  
  observeEvent(ListValuesUsers$pl,{
    output$classifiedCurve <- renderPlot(
      ListValuesUsers$pl
    )
  })
  
  
  #### LP analysis ####
  
  observeEvent(input$Cut_user,{
    
    if(!is.null(ListValuesUsers$merged_df) && input$Cut_user > min(ListValuesUsers$merged_df$Time) ){
      
      EventSU = req(input$featureEventSurv_user)
      TimeSU = req(input$featureTimeSurv_user)
      
      if(EventSU == "" || TimeSU == "" ) return()
      
      if( !all(unique(ListValuesUsers$merged_df[[EventSU]]) %in% c(0,1)) ){ 
        shinyalert("Error", "The event must be coded as 0 and 1.", type = "error")
        return()
      }
      
      DataFCMtrunc  =  ListValuesUsers$merged_df %>%
        filter(Time <= input$Cut_user )  %>% 
        group_by(ID) %>%
        mutate(L = length(ID)) 
      
      L =  length(unique(DataFCMtrunc$ID[DataFCMtrunc$L< 3]))
      
      DataFCMtrunc  =  DataFCMtrunc %>%
        filter(L>2) %>%
        select(-L) %>%
        mutate( TimeEventDiff = !!sym(TimeSU) - input$Cut_user )  %>% 
        filter( TimeEventDiff >= 0 )
      
      if(nrow(DataFCMtrunc)==0 ){
        shinyalert("Warning", "No curves are selected", type = "error")
      }else{
        ListValuesUsers$DataFCMtrunc = DataFCMtrunc
        
        # Assuming 'data' is your dataset and 'curves' is the number of curves
        truncation_summary <- paste(
          "Summary of Truncation",
          "Number of Curves:", length(unique(ListValuesUsers$merged_df$ID)),
          "Truncation Value:", input$Cut_user,
          "Number of curves <= 2 points:", L,
          "Number of curves > 2 points with TTP event before LP",
          length(unique(ListValuesUsers$merged_df$ID)) - length(unique(ListValuesUsers$DataFCMtrunc$ID)) - L,
          "Number of truncated Curves:", length(unique(ListValuesUsers$DataFCMtrunc$ID)),
          # Add more summary information if needed
          sep = "\n"
        )
        
        output$UserTruncationSummary <- renderText({
          truncation_summary
        })
      }
    }
  })
  
  observe({
    req(ListValuesUsers$merged_df)
    input$selectColor_user
    input$GoClass_user
    
    isolate({
      EventSU = req(input$featureEventSurv_user)
      TimeSU = req(input$featureTimeSurv_user)
      
      if(EventSU == "" || TimeSU == "" ) return()
      
      if(!is.numeric(ListValuesUsers$merged_df[[EventSU]]) || !is.numeric(ListValuesUsers$merged_df[[TimeSU]])){ 
        shinyalert("Error", "The variables selected for the survival analysis are not numeric.", type = "error")
        return()
      }
      
      show_modal_spinner()# show the modal window
      
      if(!is.null(ListValuesUsers$DataFCMtrunc) && dim(ListValuesUsers$DataFCMtrunc)[1] > 0 ){
        CutTime = input$Cut_user
        
        ListValuesUsers$ClassificationOUT = ClassificationNewCurves(ListValuesUsers$DataFCMtrunc[,c("ID","Observation","Time")], 
                                                                    ListValuesUsers$CONNECTORList.FCM,
                                                                    Cores = 2)
        
        ClassMatrix_entropy =  ListValuesUsers$ClassificationOUT$ClassMatrix_entropy %>%
          filter(Cluster != "Unclassified")
        
        Info_tmp = merge(ListValuesUsers$DataFCMtrunc, ClassMatrix_entropy) %>% dplyr::select(-Time,-Observation) %>% distinct()
        
        SmallNumberIDs = Info_tmp %>% group_by(Cluster) %>% dplyr::summarise(N = n()) %>% filter(N <3)
        if(dim(SmallNumberIDs)[1] > 0 )
        {
          shinyalert("Warning", paste0("The cluster(s) ", paste(SmallNumberIDs$Cluster,collapse = ", "), " can not be used in the LP analysis given the number of curves is less than 3"), "warning", 5000)
          Info_tmp = Info_tmp %>% filter(! Cluster %in% SmallNumberIDs$Cluster )
        }
        
        ############################ SURVIVAL ANALYSIS
        fit <- eval(parse(text = paste0("survfit(Surv(Info_tmp$",TimeSU,",as.numeric(Info_tmp$",EventSU,")) ~ Cluster, data = Info_tmp)")))
        
        list2env(list(Info_tmp = Info_tmp), envir = .GlobalEnv)
        ggsurv=ggsurvplot(fit = fit, data = Info_tmp,
                          xlab = "Time",  ylab = EventSU,
                          size = 1, pval = TRUE, risk.table = TRUE,conf.int = T,
                          risk.table.col="strata",ggtheme = theme_bw(),surv.median.line = "hv"  )
        
        ggsurvDataINIT = ListValuesUsers$ClassificationINIT$ggsurvData
        
        colors = ColorCluster[grep(x = names(ColorCluster),pattern =  input$tissueBox_user)]
        names(colors) = gsub(x = names(colors),replacement = "Cluster=",pattern =  paste0(input$tissueBox_user,"_"))
        
        pval = surv_pvalue(fit = fit, data = Info_tmp, get_coord = TRUE)
        
        ggsurvplot = ggplot() +    
          geom_ribbon(
            data = ggsurv$plot$data,
            aes(
              x = CutTime+time, 
              ymin = lower, 
              ymax = upper, 
              fill = strata
            ),
            alpha = 0.2  # Adjust transparency as needed
          ) +
          geom_step( data = ggsurv$plot$data, aes(x = CutTime + time, y = surv,  color = strata) ) +
          geom_point(data = ggsurv$plot$data[ggsurv$plot$data$n.censor > 0, , drop = FALSE],
                     aes(x = CutTime + time, y = surv,  color = strata), size = 4, shape = "+")+
          theme_bw()+ labs(x = "Year", fill = "Cluster", color = "Cluster")+
          scale_color_manual(values = colors)+
          scale_fill_manual(values = colors)+
          annotate("text", x = pval$pval.x, y = pval$pval.y,
                   label = pval$pval.txt, size = 5, hjust = 0)
        
        
        surv_median_values = surv_median(fit = fit)
        
        if(! all(is.na(surv_median_values$median)) ){
          ggsurvplot = ggsurvplot +
            geom_segment(
              data = surv_median_values,
              aes(
                x = CutTime+median, xend = CutTime+median, 
                y = 0, yend = 0.5, 
                color = strata
              ),
              linetype = "dashed"
            )+
            geom_segment(
              data = surv_median_values,
              aes(
                x = 0, xend = CutTime+median, 
                y = 0.5, yend = 0.5, 
                color = strata
              ),
              linetype = "dashed"
            )  +
            annotate("text", x = CutTime+surv_median_values$median, y = 0,
                     label = round(CutTime+surv_median_values$median/365,digits = 3), size = 3, hjust = 0)
        }
        
        Classification_new = ListValuesUsers$ClassificationOUT$ClassMatrix_entropy %>% select(ID,Cluster)
        meancurves = as.data.frame(ListValuesUsers$CONNECTORList.FCM$FCM$prediction$meancurves)
        colnames(meancurves) = ListValuesUsers$CONNECTORList.FCM$FCM$cluster$cluster.names
        meancurves$Time = ListValuesUsers$CONNECTORList.FCM$CONNECTORList$TimeGrid
        meancurves = meancurves %>% tidyr::gather(-Time, key = "Cluster", value = "Observation")
        
        DatasetClusters= merge(Classification_new,ListValuesUsers$DataFCMtrunc)
        
        pl1 =  DatasetClusters %>%
          ggplot()+ facet_grid(~Cluster) 
        
        if("Unclassified" %in% unique(Classification_new$Cluster)  )
        {
          pl1 = pl1 +geom_rect(
            aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
            data.frame(Cluster = 'Unclassified', Time = 0,Observation = 0, ID= -1),
            fill = 'grey', alpha = 0.5
          )
        }
        
        if(input$selectColor_user!="None"){
          DatasetClusters = DatasetClusters %>% 
            mutate_at(.vars=input$selectColor_user,.funs = as.factor )
          
          pl1 = pl1 + geom_exec(geom_line,
                                data =  DatasetClusters ,
                                x = "Time",y = "Observation", group = "ID", color = input$selectColor_user )
        } else{
          pl1 = pl1 + geom_line(aes(x = Time,y = Observation, group = ID))
        }
        
        pl1 = pl1 +
          geom_line(data = meancurves, aes(x=Time, y = Observation, col = "Mean Curves"), linewidth = 1, linetype = "dashed")+
          labs(title = paste(input$tissueBox, " classification") )+
          theme_bw()
        
        ListValuesUsers$ggsurvLast = ggsurvplot
        ListValuesUsers$pl = pl1/ggsurvplot
        
        ### change tab panel 
        updateTabsetPanel(
          inputId="UserPanels",
          selected = "panel_UserClassification"
        )
      }
      
      remove_modal_spinner() # remove it when done
    })
  })
  
  observeEvent(ListValuesUsers$pl,{
    output$classifiedCurve_user <- renderPlot(
      ListValuesUsers$pl
    )
  })
  
  observeEvent(list(input$saveLP_user,input$saveLP_fromClass_user),{
    if(!is.null(ListValuesUsers$ggsurvLast)){
      timeLP = input$Cut_user
      
      if(is.null(ListValuesUsers$LPanalysis)){
        ListValuesUsers$LPanalysis = list( pl = ListValuesUsers$ggsurvLast + labs(title = paste0("LP at ", timeLP) ) + theme(legend.position = "top")) 
      }else{
        ListValuesUsers$LPanalysis[[length(ListValuesUsers$LPanalysis)+1]] =  (ListValuesUsers$ggsurvLast + labs(title = paste0("LP at ", timeLP, " days"))  + theme(legend.position = "none"))
      }
      
      names(ListValuesUsers$LPanalysis) = paste0("plot_user",1:length(ListValuesUsers$LPanalysis))
      
      ### change tab panel 
      updateTabsetPanel(
        inputId="UserPanels",
        selected = "panel_UserSummary"
      )
    }
  })
  
  output$landmarkPoints_user <- renderUI({
    if (!is.null(ListValuesUsers$LPanalysis)) {
      for( i in seq_along(ListValuesUsers$LPanalysis) ){
        plotname <- paste("plot_user", i, sep="")
        removeUI(plotname)
      }
      
      plot_output_list <- lapply(seq_along(ListValuesUsers$LPanalysis), function(i) {
        plotname <- paste("plot_user", i, sep="")
        plotOutput(outputId = plotname)
      })
      
      # Convert the list to a tagList - this is necessary for the list of items
      # to display properly.
      do.call(tagList, plot_output_list)
      #div(plotOutput(outputId = paste("plot",length(ListValuesUsers$LPanalysis), sep="")))
    }
  })
  
  observeEvent(list(input$saveLP_fromClass_user,input$saveLP_user),{
    
    if (!is.null(ListValuesUsers$LPanalysis)) {
      lapply(seq_along(ListValuesUsers$LPanalysis), function(i) {
        p = ListValuesUsers$LPanalysis[[paste("plot_user", i, sep="")]]
        output[[paste("plot_user",i, sep="")]] <- renderPlot({ p })
      })
    }
  })
  ####
  
}