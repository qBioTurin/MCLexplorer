### ### ### ### ### ### 
### Data generation ### 
### ### ### ### ### ### 
library(dplyr)
library(readxl)

WholeData = list(MCL0208 = list(BM = list(), PB = list()),
                 EVA = list(BM = list(), PB = list()))
###  MCL208

load("MCL0208/PB/Data/Cluster/PB_mostp_clustering_4.RData")
CONNECTORList.FCM.PB -> WholeData$MCL0208$PB$CONNECTORList.FCM
df = readRDS("FeatureInfoPB.RDs")
WholeData$MCL0208$PB$Dataset <- df

load("MCL0208/BM/Data/Cluster/BM_mostp_clustering_4.RData")
CONNECTORList.FCM.BM ->  WholeData$MCL0208$BM$CONNECTORList.FCM
df = readRDS("FeatureInfoBM.RDs")
WholeData$MCL0208$BM$Dataset <- df

###  EVA

df <- haven::read_sas("./EVA/Data/connectormrdforexport20230425.sas7bdat")
legend <- read_excel("./EVA/Data/connectorMRDdata_legend.xlsx")
subdf = df %>% dplyr::select(-R1MELDUNG, -INI, -INIERG, -INIERGDAT, -PBSCT,
                             -RETRANDAT, -LREMDATE, -REZIDIV, -DEATH, -rnd1,
                             -survr1, -event, -ttp, -ttpevent,-sensitivity,-QR)
InfoComparison = df %>% dplyr::select(patid, R1MELDUNG, INI, INIERG, INIERGDAT, PBSCT,RETRANDAT, LREMDATE, REZIDIV, DEATH, rnd1, survr1, event, ttp, ttpevent,sensitivity,QR)
Info=as.data.frame(InfoComparison %>% select(patid,INIERG,rnd1, ttpevent, ttp))

colnames(Info)=c("ID","INIERG","rnd1", "ttpevent", "ttp")
WholeData$EVA$Info = Info

classificationFCM = readRDS("./EVA/RData/dataNewBMForClassification.RDs")
WholeData$EVA$BM$Dataset= left_join(classificationFCM, Info, "ID")

classificationFCM = readRDS("./EVA/RData/dataNewPBForClassification.RDs")
WholeData$EVA$PB$Dataset= left_join(classificationFCM, Info, "ID")

saveRDS(WholeData,file = "ShinyApp/ins/Data/WholeData.RDs")
