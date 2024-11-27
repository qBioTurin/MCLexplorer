library(connector)
library(ggplot2)
library(viridis)
library(labeling)
library(latex2exp)
library(grid)
library(patchwork)
library(dplyr)

Steps <- c("BASELINE","RCHOP","ARAC","ASCT",paste0("M",seq(6,36,6)))
plot.genaration = function(df,col){
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
  ### One column with white bk and color depending on the circos
  pl1 = ggplot(df,aes(x = Time, y = Observation, group=ID))+
    geom_line(aes(col=as.factor(TTP))) +
    facet_wrap(~Cluster,ncol = 1) +
    scale_color_manual(values = c("0" = "blue","1"="red"))+
    labs(x = "", y = "") +
    scale_y_continuous(limits=c(-1, 8) ,
                       breaks = seq(0,8,2),
                       labels = c("NEG","POS",TeX("$10^{-3}$"),TeX("$10^{-2}$"),TeX("$10^{-1}$")) ) +
    theme_bw() +
    theme(legend.position = "none",
          strip.text = element_text(face="bold",color = "white"),
          plot.margin = unit(c(0, 0,0,0), "cm"))
  pl = pl1/DisTimeBX + plot_layout(heights = c(4,0.5))
  g = patchworkGrob(pl)
  ### tables
  df$TTP=as.factor(df$TTP)
  levels(df$TTP)=c("NR", "R")
  CountsTTP = df %>% mutate(Cluster = as.factor(Cluster), TTP = as.factor(TTP)) %>% select(ID, Cluster,TTP) %>%
    distinct() %>% group_by(TTP,Cluster,.drop = FALSE) %>% count() %>% ungroup()
  ttheme = gridExtra::ttheme_minimal(
    core = list(bg_params = list(fill = c("#779ECB", "#F9665E"), col = "black", lwd = 3)),
    colhead = list(bg_params = list(fill = "grey", col = "black", lwd = 2)),
    rowhead = list(bg_params = list(fill = "white", col = "black", lwd = 2))
  )
  lay <- rbind(c(1,1,1,2),
               c(1,1,1,3),
               c(1,1,1,4),
               c(1,1,1,5),
               c(1,1,1,NA))
  g = gridExtra::grid.arrange(
    g ,
    gridExtra::tableGrob(CountsTTP %>% filter(Cluster=="A") %>% select(-Cluster),
                         rows=NULL,theme = ttheme),
    gridExtra::tableGrob(CountsTTP %>% filter(Cluster=="B") %>% select(-Cluster),
                         rows=NULL,theme = ttheme),
    gridExtra::tableGrob(CountsTTP %>% filter(Cluster=="C") %>% select(-Cluster),
                         rows=NULL,theme = ttheme),
    gridExtra::tableGrob(CountsTTP %>% filter(Cluster=="D") %>% select(-Cluster),
                         rows=NULL,theme = ttheme), layout_matrix = lay
  )
  ##### Change color
  #g <- ggplot_gtable(ggplot_build(pl))
  #g$layout$name
  #stript <- which(grepl('strip-t',  g$grobs[[1]]$layout$name))
  i=3
  j <- which(grepl('rect',  g$grobs[[1]]$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[1]]$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- col[ g$grobs[[1]]$grobs[[i]]$grobs[[1]]$children[[2]]$children[[1]]$label]
  stript <- which(grepl('strip-t', g$grobs[[1]]$grobs[[18]]$layout$name))
  k <- 1
  for (i in stript) {
    j <- which(grepl('rect', g$grobs[[1]]$grobs[[18]]$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[1]]$grobs[[18]]$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- col[g$grobs[[1]]$grobs[[18]]$grobs[[i]]$grobs[[1]]$children[[2]]$children[[1]]$label]
    k <- k+1
  }
  return(g)
}

plotinfo = function(tissue,Dataset){
  
  list2env(list(Dataset = Dataset), envir = .GlobalEnv)
  
  
  Dataset = Dataset %>% mutate(ArmLabel = ifelse(Arm == 1,"Lenalidomide \nTherapy", "No Therapy"))
  print(table(Dataset$Cluster,Dataset$ArmLabel))
  
  fit <- eval(parse(text = paste0("survfit(Surv(Dataset$TimeTTPevent,Dataset$TTP) ~ ArmLabel, data = Dataset)")))
  surv_median <- as.vector(summary(fit)$table[, "median"])
  df <- data.frame(x1 = surv_median, x2 = surv_median,
                   y1 = rep(0, length(surv_median)), y2 = rep(0.5, length(surv_median)))
  ggsurv=ggsurvplot(fit = fit, data = Dataset,
                    xlab = "Year",  ylab = "TTP",
                    size = 1, pval = TRUE, risk.table = TRUE,
                    risk.table.col="strata",ggtheme = theme_bw() )
  
  ggsurvList = lapply( unique(Dataset$Cluster), function(cl,Dataset){
    print(cl)
    Dataset_cl = Dataset %>% filter(Cluster == cl)
    
    fit <- eval(parse(text = paste0("survfit(Surv(Dataset_cl$TimeTTPevent,Dataset_cl$TTP) ~ ArmLabel, data = Dataset_cl)")))
    surv_median <- as.vector(summary(fit)$table[, "median"])
    df <- data.frame(x1 = surv_median, x2 = surv_median,
                     y1 = rep(0, length(surv_median)),
                     y2 = rep(0.5, length(surv_median)) )
    list2env(list(Dataset_cl = Dataset_cl), envir = .GlobalEnv)
    
    ggsurv=ggsurvplot(fit = fit, data = Dataset_cl,
                      xlab = "Year",  ylab = "TTP",
                      size = 1, pval = TRUE, risk.table = TRUE,
                      risk.table.col="strata",ggtheme = theme_bw() )
    
    ggsurv$plot + labs(title = paste0("Cluster ", cl) )
  },Dataset)
  
  return(list(SUall = ggsurv, SUclusters = ggsurvList))
  
}
