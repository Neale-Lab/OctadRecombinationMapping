
library(plyr)

#setwd("/mnt/nfs2/gdsc/mc482/Event_Table_combining_and_imaging/Output_Files")
master1<-read.delim("masterAEventTable_woflank_29-11-17")
master2<-read.delim("masterAEventTable_TCM")
master3<-read.delim("masterAEventTable_TS")
master<-rbind.fill(master1,master2,master3)
#master <- masterB

master<-read.delim("MasterAEventTable")

master$LCO <- as.numeric(master$LCO)
master$LCO[is.na(master$LCO)] <- 0
master$LNCO <- as.numeric(master$LNCO)
master$LNCO[is.na(master$LNCO)] <- 0

Meiosis <- unique(master$Meiosis)
CCT<- data.frame(Meiosis)
CCT$Meiosis <- as.character(CCT$Meiosis)
CCT$COs <- 0
CCT$NCOs <- 0
CCT$dCOs <- 0
CCT$dNCOs <- 0
CCT$CO_NCO <- 0
CCT$COlessChroms <- 0
CCT$WhichC <- 'None'
CCT$NCOlessChroms <- 0
CCT$WhichN <- 'None'
CCT$EZs <- 0
CCT$SOs <- 0
CCT$COs_without_tracts <- 0
CCT$NCOs_hDNA_only <- 0
CCT$CO2 <-0
CCT$CO3 <-0
CCT$CO4 <-0
CCT$NCO1 <-0
CCT$NCO2s <-0
CCT$NCO2ns <-0
CCT$NCO3 <-0
CCT$NCO4 <-0


for (j in 1:nrow(CCT)){
  subset1x <- master[which(master$Meiosis == CCT[j,1] & master$LCO == 1 & master$LNCO == 0),]
  if(nrow(subset1x)>0){CCT[j,2] <- nrow(subset1x)}
  subset2x <- master[which(master$Meiosis == CCT[j,1] & master$LNCO == 1  & master$LCO == 0),]
  if(nrow(subset2x)>0){CCT[j,3] <- nrow(subset2x)}
  
  #count number of dCOs
  subset2 <- master[which(master$Meiosis == CCT[j,1] & master$LCO == 2),]
  if(nrow(subset2)>0){CCT[j,4] <- nrow(subset2)}
  #count number of dNCOs
  subset2.1 <- master[which(master$Meiosis == CCT[j,1] & master$LNCO == 2),]
  if(nrow(subset2.1)>0){CCT[j,5] <- nrow(subset2.1)}
  #count number of CO+NCOs
  subset2.2 <- master[which(master$Meiosis == CCT[j,1] & master$LCO > 0 & master$LNCO > 0),]
  if(nrow(subset2.2)>0){CCT[j,6] <- nrow(subset2.2)}
  
  count=0
  string='Chr'
  count2=0
  string2='Chr'
for (i in 1:16){ 
  subset <- master[which(master$Meiosis == CCT[j,1] & master$chr == i),]
  #############CO-less and NCO-less chromos#############
  if(nrow(subset)>0){
  total <- sum(subset$LCO)
  total2 <- sum(subset$LNCO)
  if(total==0){count <- count +1}
  if(total==0){string <- paste(string, i, sep = " ")}
  if(total2==0){count2 <- count2 +1}
  if(total2==0){string2 <- paste(string2, i, sep = " ")}
  }
  if(nrow(subset)==0){count<- count+1}
  if(nrow(subset)==0){string <- paste(string, i, sep = " ")}
  if(nrow(subset)==0){count2<- count2+1}
  if(nrow(subset)==0){string2 <- paste(string2, i, sep = " ")}
  
}
CCT[j,7] <- count
if(count >0) {CCT[j,8] <- string}
CCT[j,9] <- count2
if(count2 >0) {CCT[j,10] <- string2}

  #events with 8-0
  subset3 <- master[grep("8", master$type), ] # Type contains 8 (i.e. 8:0 and 0:8 events)
  subset3 <- subset3[which(subset3$Meiosis == CCT[j,1]),]
  if(nrow(subset3)>0){CCT[j,11] <- nrow(subset3)}
  #events with 7-1
  subset6 <- master[grep("7", master$type), ] # Type contains 7 (i.e. 7:1 and 1:7 events)
  subset6 <- subset6[which(subset6$Meiosis == CCT[j,1]),]
  if(nrow(subset6)>0){CCT[j,12] <- nrow(subset6)}
  #cos with length 0
  subset4 <- master[which(master$Meiosis == CCT[j,1] & master$LCO > 0 & master$len_min == 0),]
  if(nrow(subset4)>0){CCT[j,13] <- nrow(subset4)}
  #no of NCOs composed only of hdNA
  subset5 <- master[intersect(grep("3", master$type),grep("1|2|4",master$type,invert=TRUE)),] # Type contains 3 (i.e. 3:5 and 5:3 events), but not any other segments
  subset5 <- subset5[which(subset5$Meiosis == CCT[j,1]),]
  if(nrow(subset5)>0){CCT[j,14] <- nrow(subset5)}
  
  #count COs affecting 2 chromatids
  subset6 <- master[which(master$Meiosis == CCT[j,1] & master$LCO > 0 & master$chromatids == "2_nonsis"),]
  if(nrow(subset6)>0){CCT[j,15] <- nrow(subset6)}
  #count COs affecting 3 chromatids
  subset7 <- master[which(master$Meiosis == CCT[j,1] & master$LCO > 0 & master$chromatids == 3),]
  if(nrow(subset7)>0){CCT[j,16] <- nrow(subset7)}
  #count COs affecting 4 chromatids
  subset8 <- master[which(master$Meiosis == CCT[j,1] & master$LCO > 0 & master$chromatids == 4),]
  if(nrow(subset8)>0){CCT[j,17] <- nrow(subset8)}
  #count NCOs affecting 1 chromatid
  subset9 <- master[which(master$Meiosis == CCT[j,1] & master$LNCO > 0 & master$LCO == 0 & master$chromatids == 1),]
  if(nrow(subset9)>0){CCT[j,18] <- nrow(subset9)}
  #count NCOs affecting 2_nonsis chromatids
  subset10 <- master[which(master$Meiosis == CCT[j,1] & master$LNCO > 0 & master$LCO == 0 & master$chromatids == "2_nonsis"),]
  if(nrow(subset10)>0){CCT[j,19] <- nrow(subset10)}
  #count NCOs affecting 2_sis chromatids
  subset11 <- master[which(master$Meiosis == CCT[j,1] & master$LNCO > 0 & master$LCO == 0 & master$chromatids == "2_sis"),]
  if(nrow(subset11)>0){CCT[j,20] <- nrow(subset11)}
  #count NCOs affecting 3 chromatids
  subset12 <- master[which(master$Meiosis == CCT[j,1] & master$LNCO > 0 & master$LCO == 0 & master$chromatids == 3),]
  if(nrow(subset12)>0){CCT[j,21] <- nrow(subset12)}
  #count NCOs affecting 4 chromatids
  subset13 <- master[which(master$Meiosis == CCT[j,1] & master$LNCO > 0 & master$LCO == 0 & master$chromatids == 4),]
  if(nrow(subset13)>0){CCT[j,22] <- nrow(subset13)}



}

write.table((CCT),file="COlessChromosomesTable.txt", col.names =TRUE, row.names = FALSE, quote = FALSE, sep="\t")

##events per chromo###

Meiosis <- unique(master$Meiosis)
CCT1<- data.frame(Meiosis)
CCT1$Meiosis <- as.character(CCT1$Meiosis)
CCT2<- data.frame(Meiosis)
CCT2$Meiosis <- as.character(CCT2$Meiosis)

for(x in 1:length(Meiosis)){
  subset1 <- master[which(master$Meiosis == Meiosis[x] & master$LCO > 0),]
  subset2 <- master[which(master$Meiosis == Meiosis[x] & master$LNCO > 0),]

for(i in 1:16){
  subset1.1 <- subset1[which(subset1$chr == i),]
  subset2.1 <- subset2[which(subset2$chr == i),]
  if(nrow(subset1.1)>0){CCT1[x,i+1] <- sum(subset1.1$LCO)}
  if(nrow(subset2.1)>0){CCT2[x,i+1] <- sum(subset2.1$LNCO)}
}}


write.table((CCT1),file="COperChrom2", col.names =TRUE, row.names = FALSE, quote = FALSE, sep="\t")
write.table((CCT2),file="NCOperChrom2", col.names =TRUE, row.names = FALSE, quote = FALSE, sep="\t")
