event_sorter <- function(b, u, Pan, filename){
  
  b$nb_seg <- 0
  b$groupe <- 0
  
  #split file into three groups: sortable NCOs, sortable COs, anything else
  
  #COs to sort - consider only the CO events with 2 non-sister chromatids.
  b1<-b[which(b$CO==1 & b$chromatids=="2_nonsis"),]
  #NCOs to sort - consider only the non-CO events with 1 chromatid or 2 non-sister chromatids
  b2<-b[which(b$CO==0 & (b$chromatids=="1" | b$chromatids=="2_nonsis")),]
  
  #remove sortable entries from original list
  b3 <- b[which(b$CO>1),] #more than 1 CO
  if(nrow(b3)>0){b3$classe <-3}
  b4 <- b[which(b$CO==1 & b$chromatids!="2_nonsis"),] 
  if(nrow(b4)>0){b4$classe <-2}
  b5 <- b[which(b$CO==0 & (b$chromatids=="2_sis" | b$chromatids=="3"| b$chromatids=="4")),] 
  if(nrow(b5)>0){b5$classe <-1}
  
  #u <- rename(u, c("CO_NCO"="CO"))
 
  b6 <- u[which(u$CO=="U"),] 
  if(nrow(b6)>0){
  b6$groupe <-0
  b6$seg_len_mid <-0
  b6$seg_len_min <-0
  b6$seg_len_max <-0
  b6$seg_start5 <-0
  b6$seg_stop5 <-0
  b6$seg_start3 <-0
  b6$seg_stop3 <-0
  b6$nb_seg <-0
  b6$classe <-4
  b6$segnb_snp<-0
  b6$segspaff <-0
  b6$chr_aff <- 0
  b6 <-b6[c("id","chr", "start5", "start3","stop5", "stop3", "type","CO", "chromatids", "len_min","len_mid", "len_max", 
            "commande", "groupe", "seg_len_mid","seg_len_min","seg_len_max","seg_start5","seg_start3", "seg_stop5","seg_stop3", "nb_seg","classe", "nb_snp", "segnb_snp", "segspaff","chr_aff")]
  #To keep 'U' events from at the end of the chromosomes, use the information from the old event caller
  #and match the format of the file from the new event caller
  }
  
  #b2: sortable NCOs
  #b5: unsortable NCOs (2_sis)
  #b1: sortable COs
  #b4: unsortable COs (3 and 4 chromatids)
  #b3: Unsortable COs (multiple COs)
  #b6: Class U events - at end of chromosome, never return to 4:4 so impossible to know if they are CO or NCO
  
  ###(Run CO and NCO sorters)###
  # Classification_CO_2chr_Marg
  
  #file_event<-"Events_OMT1TN_1500_5.5.CK.HP.EZIn" # the name of the file describing the events
  #b<-read.delim(file_event)
  # We consider only the CO events with 2 non-sister chromatids.
  #b<-subset(b,CO==1 & chromatids=="2_nonsis")
  
  # We call trans hDNA the segments 5:3 and 5:3a or the segments 3:5 and 3:5a belonging to a given event.
  
  # We call "incompatible events" the events containing both (5:3 or 6:2 segments) and (3:5 or 2:6 segments).
  
  # We make 3 groups :
  # - group 1 = the events without trans hDNA nor incompatibility, 
  # - group 2 = the events with trans hDNA but without incompatibility,
  # - group 3 = the incompatible events.
  if(nrow(b1)>0){
  for (i in 1:length(b1$type)){
    s<-as.character(b1$type[i])
    b1$groupe[i]<-ifelse(length(grep("[5 6]:",s))*length(grep("[2 3]:",s)) > 0,3,
                         ifelse(length(grep("[5:3 3:5]a",s))==1,2,1))
  }}
  
  # Inside each group, we classify the events into different classes :
  if(nrow(b1)>0){
  for (i in 1:length(b1$type)){
    s<-as.character(b1$type[i])
    len<-strsplit(s,"_")
    b1$nb_seg[i]<-length(len[[1]])
    
    if (b1$groupe[i]==1){
      b1$classe[i]<-ifelse(s=="(4:4aCO)",100,
                           ifelse(s=="(5:3)_(4:4aCO)",1,
                                  ifelse(s=="(3:5)_(4:4aCO)",2,
                                         ifelse(s=="(5:3)_(4:4)_(4:4aCO)",3,
                                                ifelse(s=="(4:4aCO)_(5:3)_(4:4a)",4,
                                                       ifelse(s=="(3:5)_(4:4)_(4:4aCO)",5,
                                                              ifelse(s=="(4:4aCO)_(3:5)_(4:4a)",6,
                                                                     
                                                                     ifelse(s=="(5:3)_(6:2)_(4:4aCO)",7,
                                                                            ifelse(s=="(6:2)_(5:3)_(4:4aCO)",8,
                                                                                   ifelse(s=="(3:5)_(2:6)_(4:4aCO)",9,
                                                                                          ifelse(s=="(2:6)_(3:5)_(4:4aCO)",10,
                                                                                                 ifelse(length(grep("4:4.i",s))==1,11,
                                                                                                        
                                                                                                        ifelse(s=="(6:2)_(4:4aCO)",15,
                                                                                                               ifelse(s=="(2:6)_(4:4aCO)",16,12))))))))))))))}
    
    if (b1$groupe[i]==3){
      b1$classe[i]<-ifelse(length(grep("4:4.i",s))==1,21,20)}
    
    if (b1$groupe[i]==2){
      b1$classe[i]<-ifelse(length(grep("5:3\\)_\\(5:3a",s))==1,30,
                           ifelse(length(grep("3:5\\)_\\(3:5a",s))==1,31,32))}
  }
  
  # The b file is then ordered by 1) class, 2) nb of segments and 3) len_mid.
  
  b1<-b1[order(b1$groupe, b1$classe,b1$nb_seg,b1$len_mid),]
  }
  ####NCO sorter####
  # Classification_NCO
  
  #file_event<-"Events_OMT1TN_1500_5.5.CK.HP.EZIn" # the name of the file describing the events
  #b2<-read.delim(file_event)
  # We consider only the non-CO events with 1 chromatid or 2 non-sister chromatids
  #b2<-subset(b2,CO==0 & (chromatids=="1" | chromatids=="2_nonsis"))
  
  # Élimination du (4:4) final et des parenthèses.
  # Elimination of ( 4: 4) final and parentheses.
  
  b2$type<-gsub("_\\(4:4\\)$","",b2$type)
  b2$type<-gsub("\\(","",b2$type)
  b2$type<-gsub("\\)","",b2$type)
  
  # We call trans hDNA the segments 5:3 and 5:3a or the segments 3:5 and 3:5a belonging to a given event.
  
  # We make 3 groups :
  # - group 1 with chromatids=="1" but without TRANS hDNA,
  # - group 2 with chromatids=="1" and with TRANS hDNA,
  # - group 3 with chromatids = "2_nonsis".
  
  for (i in 1:length(b2$type)){
    s<-as.character(b2$type[i])
    len<-strsplit(s,"_")
    b2$nb_seg[i]<-length(len[[1]])
    b2$groupe[i]<-ifelse(b2$chromatids[i]=="2_nonsis",3,
                         ifelse(length(grep("[5:3 3:5]a",s))==1,2,1))}
  
  # Inside each group, we classify the events into different classes :
  
  for (i in 1:length(b2$type)){
    s<-as.character(b2$type[i])
    
    if (b2$groupe[i]==1){
      b2$classe[i]<-ifelse(s=="5:3",1,
                           
                           ifelse(s=="3:5",2,
                                  ifelse(length(grep("^5:3",s))*length(grep("5:3$",s))==1,3,
                                         ifelse(length(grep("^3:5",s))*length(grep("3:5$",s))==1,4,
                                                ifelse(s=="6:2",10,
                                                       ifelse(s=="2:6",10.1,
                                                              ifelse(length(grep("6:2",s))==1,7,
                                                                     ifelse(length(grep("2:6",s))==1,8,9))))))))} 
    
    if (b2$groupe[i]==2) {
      b2$classe[i]<-ifelse(s=="5:3_5:3a",11,
                           ifelse(s=="3:5_3:5a",12,
                                  ifelse(length(grep("5:3_5:3a",s))==1,13,
                                         ifelse(length(grep("3:5_3:5a",s))==1,14,
                                                ifelse(s=="5:3_4:4_5:3a",15,
                                                       ifelse(s=="3:5_4:4_3:5a",16,
                                                              ifelse(length(grep("5:3_4:4_5:3a",s))==1,17,
                                                                     ifelse(length(grep("3:5_4:4_3:5a",s))==1,18,
                                                                            ifelse(s=="5:3_6:2_5:3a",19,
                                                                                   ifelse(s=="3:5_2:6_3:5a",20,
                                                                                          ifelse(length(grep("5:3_6:2_5:3a",s))==1,21,
                                                                                                 ifelse(length(grep("3:5_2:6_3:5a",s))==1,22,23))))))))))))}	
    
    if (b2$groupe[i]==3) {
      b2$classe[i]<-ifelse(length(grep("4:4.i",s))==1,31,30)}
    
  }
  
  # The b file is then ordered by 1) class, 2) nb of segments and 3) len_mid.
  
  b2<-b2[order(b2$groupe, b2$classe,b2$nb_seg,b2$len_mid),]
  
  ####Continue####
  
  b7 <- rbind(b2, b1, b5, b4, b3, b6) #adding the sets back together
  #b2: sortable NCOs
  #b1: sortable COs
  #b5: unsortable NCOs (2_sis)
  #b4: unsortable COs (3 and 4 chromatids)
  #b3: Unsortable COs (multiple COs)
  #b6: Unclassified events
  
  #############################################################
  # Calculate Spo11 hits per event here and add column
  #Loop through events and calculate Pan HpM for each event, then add to table b6$PanHpM. Looping is slow, but I don't have a better solution.
  for (i in 1:nrow(b7)){
    Pan.1=subset(Pan, Chr==b7[i,"chr"] & Pos>=b7[i,"start5"] & Pos <=b7[i,"stop3"])
    b7[i,"PanHpM"]=sum(Pan.1$TotalHpM)
  }
  
  b7$PanHpM=round(b7$PanHpM,2) # Round PanHpM output to 2 decimals
  
  # Now some extra bits:
  b7$GenomeHpM=round(1/12.01*b7$len_max,2) # What is the expected HpM for an interval this size (assuming uniform Spo11-oligo hits across the 12.01 Mbp genome)?
  b7$Obs_Exp_HpM=round(b7$PanHpM/b7$GenomeHpM,2) # What is the fold difference between observed and expected association with Spo11 hits for each event?
  b7$GenomeHpM=round(1/12.01*b7$len_max,0) #round to nearest whole number before plotting
  #############################################################
  
  ##Auto-classification##
  ##Suggested values only: Adjust manually during visual inspection##
  
  b7$Notes <- "" #blank space to write in any notes during inspection
  
  #min and max number of breaks necessary to explain the event. note: the breaks are not necessarily formed by Spo11.
  b7$CutMin <-ifelse(b7$CO==1 & b7$groupe==3, 2, #COs with incompatible hDNA
           ifelse(b7$groupe==0 & (b7$classe==1 | b7$classe==2 | b7$classe==3), 2, #NCO 2sis, CO 3/4 chromatids, multiple COs
           1))
  
  b7$CutMax <-ifelse(b7$CO ==0 & b7$groupe==3,2, #NCOs with incompatible hDNA
           ifelse(b7$CO==1 & (b7$classe==15 | b7$classe==16 |b7$groupe==3), 2, #COs with gap (could be double cut), and COs with incompatible hDNA
           ifelse(b7$groupe==0 & (b7$chromatids=="2_sis" | b7$chromatids==3) &(b7$classe==1 | b7$classe==2), 2, #NCO 2_sis and CO with 3 chromatids
           ifelse(b7$groupe==0 & b7$chromatids==4 & b7$classe==2, 3, #Co with 4 chromatids
           ifelse(b7$groupe==0 & b7$classe==3, 2, #multiple COs on 2/3/4 chromatids
           1)))))
  
  #Likely Spo11 Breaks: if there are multiple breaks, how many are likely to be caused by Spo11? i.e. how many hit an additional chromo?
  b7$LSB <-ifelse(b7$groupe==0 & (b7$classe==1 | b7$classe==2), 2, #NCO 2_sis, CO with 3 or 4 chromatids, multiple COs
           ifelse(b7$groupe==0 & b7$classe==3 & b7$chromatids=="2_nonsis", 1, 
           ifelse(b7$groupe==0 & b7$classe==3 & (b7$chromatids=="2_sis" | b7$chromatids==3 | b7$chromatids==4), 2, #multiple COs on 3 or 4 chromatids
           1)))
  
  write.table(b7,file=sprintf("%s",filename), col.names =TRUE, row.names = FALSE, quote = FALSE, sep="\t")
  
}

