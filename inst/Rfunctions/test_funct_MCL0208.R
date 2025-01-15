test_indipendence<-function(info_tibble,variable_1,variable_2,palette){
  table<-info_tibble%>%select(all_of(c(variable_2,variable_1)))%>%table()
  table<-table[rowSums(table)!=0,colSums(table)!=0]
  if(!is.table(table)){return(NULL)}
  else{
    df<-as.data.frame(table)
    colnames(df)<-c("Var1","Var2","Count")
    mosaicplot<-ggplot(data = df) +
      geom_mosaic(aes(weight = Count, x = product(Var2), fill = Var1),alpha=1)+
      xlab(variable_1)+
      ylab(variable_2)+
      scale_fill_manual(values=palette)+
      scale_color_manual(values=palette)+
      guides(fill=guide_legend(title=variable_2),
             color=guide_legend(title=variable_2))+
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
      xlab(variable_1)+
      ylab(variable_2)+
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
}

test_distribution<-function(info_tibble,variable_1,variable_2,palette){
  df<-info_tibble%>%
    filter(!is.na(!!sym(variable_2))&!is.na(!!sym(variable_1)))%>%
    select(all_of(c(variable_2,variable_1)))
  colnames(df)<-c("variable_2","variable_1")
  df$variable_2<-as.factor(df$variable_2)
  df$variable_1<-as.numeric(df$variable_1)
  df<-df%>%
    filter(!is.na(variable_1))
  df_count<-count(df,variable_2)
  
  kw_test<-kruskal.test(variable_1 ~ variable_2, data = df)
  levene_test_result <- leveneTest(variable_1 ~ variable_2, data = df)
  shapiro_results <- df %>%
    group_by(variable_2) %>%
    mutate(count=n())%>%
    #filter(count>2)%>%
    mutate(shapiro_p = ifelse(count>2,shapiro.test(variable_1)$p.value,1))%>%
    ungroup()%>%
    select(variable_2,shapiro_p)%>%
    distinct()
  bartlett_test_result<-tryCatch(
    {bartlett.test(variable_1 ~ variable_2, data = df)},
    error = function(e) {NULL})
  anova_test_results <- aov(variable_1 ~ variable_2, data = df)
  
  boxplot<-ggplot() +
    geom_boxplot(data = df,
                 aes(x=variable_2,
                     y=variable_1,
                     fill=variable_2,
                     color=variable_2,
                     group=variable_2),
                 alpha=0.5)+
    ylab(variable_1)+
    scale_fill_manual(values=palette)+
    scale_color_manual(values=palette)+
    guides(fill=guide_legend(title=variable_2),
           color=guide_legend(title=variable_2))+
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
      ymin = max(df$variable_1), ymax = max(df$variable_1)
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
  
  return(list(#shapiro_test=shapiro_results,
    kw_test=kw_test,
    anova_test=anova_test_results,
    bartlett_test=bartlett_test_result,
    levene_test=levene_test_result,
    plot=boxplot))
}

test_distr_unified<-function(info_tibble,variable_1,variable_2,palette){
  df<-info_tibble%>%
    select(all_of(c(variable_2,variable_1)))
  colnames(df)<-c("variable_2","variable_1")
  df$variable_2<-as.factor(df$variable_2)
  df$variable_1<-as.numeric(df$variable_1)
  df_count<-count(df,variable_2)
  
  boxplot<-ggplot() +
    geom_boxplot(data = df,
                 aes(x=variable_2,
                     y=variable_1,
                     fill=variable_2,
                     color=variable_2,
                     group=variable_2),
                 alpha=0.5)+
    ylab(variable_1)+
    scale_fill_manual(values=palette)+
    scale_color_manual(values=palette)
  
  if(any(df_count$n<3)){
    kw_test<-NULL
    levene_test_result<-NULL
    shapiro_results<-NULL
    t_test_result<-NULL
    wilcox_test_result<-NULL
    f_test<-NULL
    
    for(i in 1:nrow(df_count)) {
      boxplot <- boxplot + annotation_custom(
        grob = textGrob(label = df_count$n[i],  hjust = 0.5,vjust = -2, gp = gpar(cex = 1.2)),
        xmin = i - 0.3, xmax = i + 0.3,
        ymin = max(df$variable_1,na.rm = TRUE), ymax = max(df$variable_1,na.rm = TRUE)
      )
    }
  }
  else{
    kw_test<-ks.test(variable_1 ~ variable_2, data = df)
    levene_test_result <- leveneTest(variable_1 ~ variable_2, data = df)
    shapiro_results <- df %>%
      mutate(count=n())%>%
      mutate(shapiro_p = ifelse(count>2,shapiro.test(variable_1)$p.value,1))%>%
      ungroup()%>%
      select(variable_2,shapiro_p)%>%
      distinct()
    t_test_result <- t.test(variable_1 ~ variable_2, data = df)
    wilcox_test_result <- wilcox.test(variable_1 ~ variable_2, data = df)
    f_test<-var.test(variable_1 ~ variable_2, data = df)
    
    for(i in 1:nrow(df_count)) {
      if(shapiro_results$shapiro_p[i]>0.05){pvalue_shapiro<-""}
      else if(shapiro_results$shapiro_p[i]>0.01){pvalue_shapiro<-"*"}
      else {pvalue_shapiro<-"**"}
      boxplot <- boxplot + annotation_custom(
        grob = textGrob(label = paste(df_count$n[i],pvalue_shapiro),  hjust = 0.5,vjust = -2, gp = gpar(cex = 1.2)),
        xmin = i - 0.3, xmax = i + 0.3,
        ymin = max(df$variable_1,na.rm = TRUE), ymax = max(df$variable_1,na.rm = TRUE)
      )
    }
  }
  boxplot<-boxplot+
    coord_cartesian(clip = "off") +
    guides(fill=guide_legend(title=variable_2),
           color=guide_legend(title=variable_2))+
    xlab(variable_2)+
    theme(plot.title = element_text(v = 8),
          plot.margin = margin(50, 10, 10, 10),
          plot.background = element_rect(fill = "transparent",linewidth = 0),
          panel.background= element_rect(fill = "transparent",linewidth = 0),
          legend.box.background =  element_rect(fill = "transparent",linewidth = 0),
          legend.background = element_rect(fill = "transparent",linewidth = 0),
          strip.background=element_rect(fill = "transparent",linewidth = 1))
  
  return(list(kw_test=kw_test,
              levene_test=levene_test_result,
              f_test=f_test,
              shapiro_test=shapiro_results,
              kw_test=kw_test,
              t_test=t_test_result,
              wilcox_test=wilcox_test_result,
              plot=boxplot))
}

