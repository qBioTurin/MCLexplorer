test_indipendence<-function(info_tibble,variable_name,tissue){
  if(tissue=="BM"){
    Cluster_tissue<-"Cluster_BM"
    palette<-c("#F1BB7B", "#FD6467", "#5B1A18", "#D67236")
    }
  else if(tissue=="PB"){
    Cluster_tissue<-"Cluster_PB"
    palette<-c("#440154", "#31688E", "#35B779", "#8FD744")
  }
  table<-info_tibble%>%select(all_of(c(Cluster_tissue,variable_name)))%>%table()
  table<-table[rowSums(table)!=0,colSums(table)!=0]
  df<-as.data.frame(table)
  colnames(df)<-c("Var1","Var2","Count")
  mosaicplot<-ggplot(data = df) +
    geom_mosaic(aes(weight = Count, x = product(Var2), fill = Var1),alpha=1)+
    xlab(variable_name)+
    ylab("Clusters")+
    scale_fill_manual(values=palette)+
    scale_color_manual(values=palette)+
    guides(fill=guide_legend(title=Cluster_tissue),
           color=guide_legend(title=Cluster_tissue))+
    theme(axis.title.y = element_blank(),
          axis.line = element_blank(),
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          plot.background = element_rect(fill = "transparent",linewidth = 0),
          panel.background= element_rect(fill = "transparent",linewidth = 0),
          legend.box.background =  element_rect(fill = "transparent",linewidth = 0),
          legend.background = element_rect(fill = "transparent",linewidth = 0),
          strip.background=element_rect(fill = "transparent",linewidth = 1))
  chisq_test<-chisq.test(table)
  df_residuals<-as.data.frame.table(chisq_test$residuals)
  colnames(df_residuals)<-c("Var1","Var2","Residual")
  plot<-mosaicplot+
    ggplot() +
    geom_tile(data =df_residuals ,
              aes(x = Var2, y = Var1, fill = Residual,width=0.9, height=0.9),color = "white") +
    geom_text(data =df,aes(x=Var2,y=Var1,label=Count)) +
    scale_fill_gradient2(low = "#583E60", high = "#B6A21B", mid = "white",
                         midpoint = 0, limit = c(min(df_residuals$Residual), max(df_residuals$Residual)),
                         space = "Lab", name = "Residuals")+
    xlab(variable_name)+
    ylab("Clusters")+
    theme(axis.title = element_blank(),
          axis.line = element_blank(),
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          plot.background = element_rect(fill = "transparent",linewidth = 0),
          panel.background= element_rect(fill = "transparent",linewidth = 0),
          legend.box.background =  element_rect(fill = "transparent",linewidth = 0),
          legend.background = element_rect(fill = "transparent",linewidth = 0),
          strip.background=element_rect(fill = "transparent",linewidth = 1))
  if(any(chisq_test$expected<5)|all(dim(table)==c(2,2))){
    test_type<-"Fisher"
    test<-fisher.test(table)
    plot<-plot+
      plot_annotation(
        caption = paste("p-value of the Fisher test: ",round(test$p.value,4)),
        theme = theme(plot.background = element_rect(fill = "transparent",linewidth = 0),
                      panel.background= element_rect(fill = "transparent",linewidth = 0),
                      legend.box.background =  element_rect(fill = "transparent",linewidth = 0),
                      legend.background = element_rect(fill = "transparent",linewidth = 0),
                      strip.background=element_rect(fill = "transparent",linewidth = 1))
      )
  }
  else{
    test_type<-"Chi Squared"
    test<-chisq_test
    plot<-plot+
      plot_annotation(
        caption = paste("p-value of the Chi Squared test: ",round(test$p.value,4)),
        theme = theme(plot.background = element_rect(fill = "transparent",linewidth = 0),
                      panel.background= element_rect(fill = "transparent",linewidth = 0),
                      legend.box.background =  element_rect(fill = "transparent",linewidth = 0),
                      legend.background = element_rect(fill = "transparent",linewidth = 0),
                      strip.background=element_rect(fill = "transparent",linewidth = 1))
      )
  }
  return(list(test_type=test_type,test=test,plot=plot,table=table))
}

test_unified<-function(info_tibble,variable_name,tissue){
  if(tissue=="BM"){
    Cluster_tissue<-"Cluster_BM"
    palette<-c("#F79071","#994627")
  }
  else if(tissue=="PB"){
    Cluster_tissue<-"Cluster_PB"
    palette<-c("#3B3571", "#62C75F")
  }
  table<-info_tibble%>%
    mutate(Cluster_unified=ifelse(!!sym(Cluster_tissue)=="A"|!!sym(Cluster_tissue)=="B","A/B","C/D"))%>%
    select(Cluster_unified,all_of(variable_name))%>%
    table()
  table<-table[rowSums(table)!=0,colSums(table)!=0]
  df<-as.data.frame(table)
  colnames(df)<-c("Var1","Var2","Count")
  mosaicplot<-ggplot(data = df) +
    geom_mosaic(aes(weight = Count, x = product(Var2), fill = Var1),alpha=1)+
    xlab(variable_name)+ylab("Clusters")+
    scale_fill_manual(values=palette)+
    scale_color_manual(values=palette)+
    guides(fill=guide_legend(title=Cluster_tissue),
           color=guide_legend(title=Cluster_tissue))+
    theme(axis.title.y = element_blank(),
          axis.line = element_blank(),
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          plot.background = element_rect(fill = "transparent",linewidth = 0),
          panel.background= element_rect(fill = "transparent",linewidth = 0),
          legend.box.background =  element_rect(fill = "transparent",linewidth = 0),
          legend.background = element_rect(fill = "transparent",linewidth = 0),
          strip.background=element_rect(fill = "transparent",linewidth = 1))
  chisq_test<-chisq.test(table)
  df_residuals<-as.data.frame.table(chisq_test$residuals)
  colnames(df_residuals)<-c("Var1","Var2","Residual")
  plot<-mosaicplot+
    ggplot() +
    geom_tile(data =df_residuals ,
              aes(x = Var2, y = Var1, fill = Residual,width=0.9, height=0.9),color = "white") +
    geom_text(data =df,aes(x=Var2,y=Var1,label=Count)) +
    scale_fill_gradient2(low = "#583E60", high = "#B6A21B", mid = "white",
                         midpoint = 0, limit = c(min(df_residuals$Residual), max(df_residuals$Residual)),
                         space = "Lab", name = "Residuals")+
    xlab(variable_name)+ylab("Clusters")+
    theme(axis.title = element_blank(),
          axis.line = element_blank(),
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          plot.background = element_rect(fill = "transparent",linewidth = 0),
          panel.background= element_rect(fill = "transparent",linewidth = 0),
          legend.box.background =  element_rect(fill = "transparent",linewidth = 0),
          legend.background = element_rect(fill = "transparent",linewidth = 0),
          strip.background=element_rect(fill = "transparent",linewidth = 1))
  if(any(chisq_test$expected<5)|all(dim(table)==c(2,2))){
    test_type<-"Fisher"
    test<-fisher.test(table)
    plot<-plot+
      labs(
        caption = paste("p-value of the Fisher test: ",round(test$p.value,4)),
        theme = theme(plot.background = element_rect(fill = "transparent",linewidth = 0),
                      panel.background= element_rect(fill = "transparent",linewidth = 0),
                      legend.box.background =  element_rect(fill = "transparent",linewidth = 0),
                      legend.background = element_rect(fill = "transparent",linewidth = 0),
                      strip.background=element_rect(fill = "transparent",linewidth = 1))
      )
  }
  else{
    test_type<-"Chi Squared"
    test<-chisq_test
    plot<-plot+
      plot_annotation(
              caption = paste("p-value of the Chi Squared test: ",round(test$p.value,4)),
              theme = theme(plot.background = element_rect(fill = "transparent",linewidth = 0),
                            panel.background= element_rect(fill = "transparent",linewidth = 0),
                            legend.box.background =  element_rect(fill = "transparent",linewidth = 0),
                            legend.background = element_rect(fill = "transparent",linewidth = 0),
                            strip.background=element_rect(fill = "transparent",linewidth = 1))
            )
  }
  return(list(test_type=test_type,test=test,plot=plot,table=table))
}

test_distribution<-function(info_tibble,variable_name,tissue){
  if(tissue=="BM"){
    Cluster_tissue<-"Cluster_BM"
    palette<-c("#F1BB7B", "#FD6467", "#5B1A18", "#D67236")
  }
  else if(tissue=="PB"){
    Cluster_tissue<-"Cluster_PB"
    palette<-c("#440154", "#31688E", "#35B779", "#8FD744")
  }
  df<-info_tibble%>%
    filter(!is.na(!!sym(Cluster_tissue))&!is.na(!!sym(variable_name)))%>%
    select(all_of(c(Cluster_tissue,variable_name)))
  colnames(df)<-c("Cluster_tissue","variable_name")
  df$Cluster_tissue<-as.factor(df$Cluster_tissue)
  df$variable_name<-as.numeric(df$variable_name)
  df<-df%>%
    filter(!is.na(variable_name))
  df_count<-count(df,Cluster_tissue)
  
  kw_test<-kruskal.test(variable_name ~ Cluster_tissue, data = df)
  levene_test_result <- leveneTest(variable_name ~ Cluster_tissue, data = df)
  shapiro_results <- df %>%
    group_by(Cluster_tissue) %>%
    summarise(shapiro_p = shapiro.test(variable_name)$p.value)
  bartlett_test_result<-bartlett.test(variable_name ~ Cluster_tissue, data = df)
  anova_test_results <- aov(variable_name ~ Cluster_tissue, data = df)
  
  boxplot<-ggplot() +
    geom_boxplot(data = df,
                 aes(x=Cluster_tissue,
                     y=variable_name,
                     fill=Cluster_tissue,
                     color=Cluster_tissue,
                     group=Cluster_tissue),
                 alpha=0.5)+
    ylab(variable_name)+
    scale_fill_manual(values=palette)+
    scale_color_manual(values=palette)+
    guides(fill=guide_legend(title=Cluster_tissue),
           color=guide_legend(title=Cluster_tissue))+
    theme(plot.title = element_text(v = 8),
          plot.margin = margin(30, 10, 10, 10),
          plot.background = element_rect(fill = "transparent",linewidth = 0),
          panel.background= element_rect(fill = "transparent",linewidth = 0),
          legend.box.background =  element_rect(fill = "transparent",linewidth = 0),
          legend.background = element_rect(fill = "transparent",linewidth = 0),
          strip.background=element_rect(fill = "transparent",linewidth = 1),
          axis.title.x = element_blank())
  
  for(i in 1:nrow(df_count)) {
    if(shapiro_results$shapiro_p[i]>0.05){pvalue_shapiro<-""}
    else if(shapiro_results$shapiro_p[i]>0.01){pvalue_shapiro<-"*"}
    else {pvalue_shapiro<-"**"}
    boxplot <- boxplot + annotation_custom(
      grob = textGrob(label = paste(df_count$n[i],pvalue_shapiro), hjust = 0.5,vjust = -2, gp = gpar(cex = 1.2)),
      xmin = i - 0.3, xmax = i + 0.3,
      ymin = max(df$variable_name), ymax = max(df$variable_name)
    )
  }
  
  boxplot<-boxplot+
    coord_cartesian(clip = "off")+
    theme(plot.title = element_text(v = 8),
          plot.margin = margin(30, 10, 10, 10),
          plot.background = element_rect(fill = "transparent",linewidth = 0),
          panel.background= element_rect(fill = "transparent",linewidth = 0),
          legend.box.background =  element_rect(fill = "transparent",linewidth = 0),
          legend.background = element_rect(fill = "transparent",linewidth = 0),
          strip.background=element_rect(fill = "transparent",linewidth = 1))
  #scale_fill_manual(values=c("#875E6E","#BC8198"))+
  
  return(list(kw_test=kw_test,
              levene_test=levene_test_result,
              shapiro_test=shapiro_results,
              bartlett_test=bartlett_test_result,
              anova_test=anova_test_results,
              plot=boxplot))
}

test_distr_unified<-function(info_tibble,variable_name,tissue){
  if(tissue=="BM"){
    Cluster_tissue<-"Cluster_BM"
    palette<-c("#F79071","#994627")
  }
  else if(tissue=="PB"){
    Cluster_tissue<-"Cluster_PB"
    palette<-c("#3B3571", "#62C75F")
  }
  df<-info_tibble%>%
    filter(!is.na(!!sym(Cluster_tissue))&!is.na(!!sym(variable_name)))%>%
    mutate(Cluster_unified=ifelse(!!sym(Cluster_tissue)=="A"|!!sym(Cluster_tissue)=="B","A/B","C/D"))%>%
    select(Cluster_unified,all_of(variable_name))
  colnames(df)<-c("Cluster_tissue","variable_name")
  df$Cluster_tissue<-as.factor(df$Cluster_tissue)
  df$variable_name<-as.numeric(df$variable_name)
  df<-df%>%
    filter(!is.na(variable_name))
  df_count<-count(df,Cluster_tissue)
  
  kw_test<-ks.test(variable_name ~ Cluster_tissue, data = df)
  levene_test_result <- leveneTest(variable_name ~ Cluster_tissue, data = df)
  shapiro_results <- df %>%
    group_by(Cluster_tissue) %>%
    summarise(shapiro_p = shapiro.test(variable_name)$p.value)
  t_test_result <- t.test(variable_name ~ Cluster_tissue, data = df)
  wilcox_test_result <- wilcox.test(variable_name ~ Cluster_tissue, data = df)
  f_test<-var.test(variable_name ~ Cluster_tissue, data = df)
  
  boxplot<-ggplot() +
    geom_boxplot(data = df,
                 aes(x=Cluster_tissue,
                     y=variable_name,
                     fill=Cluster_tissue,
                     color=Cluster_tissue,
                     group=Cluster_tissue),
                 alpha=0.5)+
    ylab(variable_name)+
    scale_fill_manual(values=palette)+
    scale_color_manual(values=palette)
    theme(axis.title.x = element_blank())
  
  for(i in 1:nrow(df_count)) {
    if(shapiro_results$shapiro_p[i]>0.05){pvalue_shapiro<-""}
    else if(shapiro_results$shapiro_p[i]>0.01){pvalue_shapiro<-"*"}
    else {pvalue_shapiro<-"**"}
    boxplot <- boxplot + annotation_custom(
      grob = textGrob(label = paste(df_count$n[i],pvalue_shapiro),  hjust = 0.5,vjust = -2, gp = gpar(cex = 1.2)),
      xmin = i - 0.3, xmax = i + 0.3,
      ymin = max(df$variable_name), ymax = max(df$variable_name)
    )
  }
  
  boxplot<-boxplot+
    coord_cartesian(clip = "off") +
    guides(fill=guide_legend(title=Cluster_tissue),
           color=guide_legend(title=Cluster_tissue))+
    theme(plot.title = element_text(v = 8),
          plot.margin = margin(50, 10, 10, 10),
          plot.background = element_rect(fill = "transparent",linewidth = 0),
          panel.background= element_rect(fill = "transparent",linewidth = 0),
          legend.box.background =  element_rect(fill = "transparent",linewidth = 0),
          legend.background = element_rect(fill = "transparent",linewidth = 0),
          strip.background=element_rect(fill = "transparent",linewidth = 1))
  #scale_fill_manual(values=c("#875E6E","#BC8198"))+

  return(list(kw_test=kw_test,
              levene_test=levene_test_result,
              f_test=f_test,
              shapiro_test=shapiro_results,
              kw_test=kw_test,
              t_test=t_test_result,
              wilcox_test=wilcox_test_result,
              plot=boxplot))
}
