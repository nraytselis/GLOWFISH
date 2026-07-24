#Great Lakes GLNPO Data 
#Downloaded from EPA GLENDA on 7/23/26
setwd("~/Desktop/Whitefish/GreatLakes")
GLNPO_zoop = read_csv("GLNPO_GLENDA_Zoop_7-24-26.csv")

#Taxanomic Codes 
#CAL	Calanoid copepod adults (excluding Limnocalanus macrurus)
#CALIM	Calanoid copepod copepodites (excluding Limnocalanus macrurus)
#LIMNO	Limnocalanus macrurus adults
#LIMNOIM	Limnocalanus macrurus copepodites
#CYC	Cyclopoid copepod adults
#CYCIM	Cyclopoid copepod copepodites
#DAP	Daphnia
#BOS	Bosminid cladocerans
#HOLGIBB	Holopedium gibberum
#CLAOTH	Other herbivorous cladocerans
#BYTLONG	Bythotrephes longimanus
#CECPENG	Cercopagis pengoi
#NPRED	Native predatory cladocerans

#some columns are characters so convert to numeric
GLNPO_zoop$"BYTLONG_Num/m3" = as.numeric(GLNPO_zoop$"BYTLONG_Num/m3")
GLNPO_zoop$"BYTLONG_ugDW/m3" = as.numeric(GLNPO_zoop$"BYTLONG_ugDW/m3")
#Replace rows with NA with zero since no Bythotrephes longimanus observed 
GLNPO_zoop$"BYTLONG_Num/m3"[is.na(GLNPO_zoop$"BYTLONG_Num/m3")] <- 0
GLNPO_zoop$"BYTLONG_ugDW/m3"[is.na(GLNPO_zoop$"BYTLONG_ugDW/m3")] <- 0

##Lump into three groups: calanoids, cyclopoids, cladocerans 

#Calanoids:
#CAL, CALIM, LIMNO, LIMNOIM

#Cyclopoids
#CYC, CYCIM

#Cladocerans
#DAP, BOS, HOLGIBB, CLAOTH, BYTLONG, CECPENG, NPRED 

GLNPO_zoop$TotalCals <- rowSums(GLNPO_zoop[, 9:16])
GLNPO_zoop$TotalCyclos <- rowSums(GLNPO_zoop[, 17:20])
GLNPO_zoop$TotalClads <- rowSums(GLNPO_zoop[, 21:34])

#Break up by Lake
list_of_dfs_GLNPO_zoop <- split(GLNPO_zoop, GLNPO_zoop$LAKE)
Erie <- list_of_dfs_GLNPO_zoop$ER
Huron <- list_of_dfs_GLNPO_zoop$HU
Michigan <- list_of_dfs_GLNPO_zoop$MI
Ontario <- list_of_dfs_GLNPO_zoop$ON
Superior <- list_of_dfs_GLNPO_zoop$SU

#Bring in Taylor's other Data






