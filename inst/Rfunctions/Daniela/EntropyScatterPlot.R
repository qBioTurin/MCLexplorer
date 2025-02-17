
tissueMCL= "BM"

CONNECTORList.FCM = WholeData$MCL0208[[tissueMCL]]$CONNECTORList.FCM

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

  maxL = max(df1$L)
  maxE = max(df1Entropy)
  minL = min(df1$L)
  minE = min(df1$Entropy)
  
  plot =  df1 %>%
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
}