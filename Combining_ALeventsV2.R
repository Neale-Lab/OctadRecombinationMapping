#Marg, 16-09-16
#Combines annotated events tables, removes false events, adds NCO information
#Adds Spo11 hit information needed for double cut analysis
#V2 04-12-17
#Automatically annotates events with 'B_type' as described by Marsolier-Kergoat et al., 2017

#Read in the Spo11 datatracks from Pan et al and Mohibullah et al
Pan <- read.delim("FullMap.Cer3_Pan_HA_1_4h_c.txt")
NeemanWT <- read.delim("WT_SCFullMap_Neeman.txt")

TotalHits=sum(Pan$Watson+Pan$Crick)
Pan$TotalHpM=(Pan$Watson+Pan$Crick)/TotalHits*1000000 # Convert Pan data to Hits per million reads (HpM)

TotalHits=sum(NeemanWT$Watson+NeemanWT$Crick)
NeemanWT$TotalHpM=(NeemanWT$Watson+NeemanWT$Crick)/TotalHits*1000000 # Convert Pan data to Hits per million reads (HpM)

# Create new empty master dataframe
masterA=NULL

library(stringr)
#library(data.table)

# Sequentially imports each Event Table and combines into a Master Table with additional identified columns for Meiosis, Threshold, Genotype and GenotypeID (e.g. OM2, 1500, OM, 2)
files = list.files(pattern="AEvents_") # import files names with "AEvents_" string into variable "files" 
files1 = length(files) # Count number of files
files2 = read.table(text = files, sep = "_", as.is = TRUE) #Split file names by "_" separator and create table "files2"

for (j in 1:files1) {
 data=NULL 
 data <- read.table(files[j], sep = "\t", header=TRUE)
 if(ncol(data)==35){
 data <- data[,-35]
 data <- data[,-34]}
 data["Meiosis"]=files2[j,2] # Insert Meiosis identifier from files2 into a new column
 data["threshold"]=str_extract(files2[j,3],"[[:digit:]]+") # Extracts digit portion of string
 data["Genotype"]=gsub('([A-z]+).*','\\1',files2[j,2])
 data["GenotypeID"]=gsub('.*([0-9]+)','\\1',files2[j,2])
 data["midpoint"]=(data$start5+data$start3+data$stop5+data$stop3)/4 # Add midpoint column
 data["debut"]=(data$midpoint-(0.5*data$len_mid)) # Add event start column
 data["fin"]=(data$midpoint+(0.5*data$len_mid)) # Add event end column
 masterA <- rbind(masterA, data) # Combine data into master table

}

#Remove mitotic or false events
masterA <- masterA[which(masterA$LSB!=0),]
masterA <- masterA[which(masterA$threshold==1500),]

masterA$CO[masterA$CO == "U"] <- NA
masterA$CO <- as.numeric(as.character(masterA$CO))

masterA$LCO <- masterA$CO
masterA$LNCO <- (masterA$LSB-masterA$CO)
masterA$MaxNCO <- (masterA$CutMax-masterA$CO)
masterA$MinNCO <- (masterA$CutMin-masterA$CO)

masterA$midpoint1 <- 0
masterA$midpoint2 <- 0
masterA$midpoint3 <- 0

masterA$midpoint1 <-ifelse(masterA$LSB ==2, masterA$midpoint-(masterA$len_mid/4),
                           ifelse(masterA$LSB ==3, masterA$midpoint-(masterA$len_mid/3),
                                  masterA$midpoint))

masterA$midpoint2 <-ifelse(masterA$LSB ==2, masterA$midpoint+(masterA$len_mid/4),
                           ifelse(masterA$LSB ==3, masterA$midpoint,
                                  NA))

masterA$midpoint3 <-ifelse(masterA$LSB ==3, masterA$midpoint+(masterA$len_mid/3),
                           NA)


#masterA <- rename(masterA, c("segspaff"="SpAff"))
#masterA <- rename(masterA, c("segnb_snp"="seg_nb_SNP"))

masterA$seg_nb_SNP <- 0
masterA$SpAff <- 0

# Reorder columns
masterA=masterA[c("id","Meiosis","threshold","Genotype","GenotypeID","chr","start5",
                  "start3","stop5","stop3","debut","fin","type",
                  "chromatids","len_min","len_mid","len_max", "nb_seg", "seg_len_mid", "seg_len_min", "seg_len_max",	
                  "seg_start5", "seg_start3",	"seg_stop5", "seg_stop3",  "PanHpM",	"GenomeHpM",	"Obs_Exp_HpM", 
                  "nb_snp", "seg_nb_SNP", "SpAff", "chr_aff", "groupe", "classe", "LSB", "LCO", "LNCO", 
                  "CutMin", "CutMax", "MaxNCO", "MinNCO", "midpoint", "midpoint1", "midpoint2", "midpoint3",  
                  "commande", "Notes")]

masterA$B_type <- 'Other'
masterA$classe <- as.numeric(masterA$classe)
masterA$groupe <- as.numeric(masterA$groupe)
masterB <- NULL

switch <- 'off'
for(x in 1:nrow(masterA)){
  df <- masterA[x,]
  if((df$type) =='05:03'){df$type <- '5:3'}
  if((df$type) =='03:05'){df$type <- '3:5'}
  if((df$type) =='06:02'){df$type <- '6:2'}
  if((df$type) =='02:06'){df$type <- '2:6'}
  if(is.na(df$LCO) ==FALSE){ #if there isn't an NA in LCO (which is the case for U events),
      if(df$LCO == 0 & df$groupe==1 & (df$classe ==10|df$classe ==10.1|df$classe ==7|df$classe ==8|df$classe ==9)){df$B_type <-  'Full Conversion'}     
      if(df$LCO == 0 & df$groupe==1 & (df$classe ==1|df$classe ==2)){df$B_type <- 'One-Sided'} 
      if(df$LCO == 0 & df$groupe==1 & (df$classe ==3|df$classe ==4)){df$B_type <- 'One-Sided + patch'} 
      if(df$LCO == 0 & df$groupe==2 & (df$classe ==11|df$classe ==12|df$classe ==15|df$classe ==16|df$classe ==19|df$classe ==20)){df$B_type <- 'Two-sided'} 
      if(df$LCO == 0 & df$groupe==2 & (df$classe ==13|df$classe ==14|df$classe ==17|df$classe ==18|df$classe ==21|df$classe ==22|df$classe ==23)){df$B_type <- 'Two-sided + patch'} 
      if(df$LCO == 0 & df$groupe==3 & (df$classe ==30|df$classe ==31)){df$B_type <- 'Two Chr'} 
      if(df$LNCO == 0 & df$groupe==1 & df$classe ==100){df$B_type <- 'No Strand Transfer'}     
      if(df$LNCO == 0 & df$groupe==1 & df$classe ==11){df$B_type <- 'Sym 4:4'}      
      if(df$LNCO == 0 & df$groupe==1 & (df$classe ==15|df$classe ==16)){df$B_type <- 'Full Conversion'} 
      if(df$LNCO == 0 & df$groupe==1 & (df$classe ==7|df$classe ==8|df$classe ==9|df$classe ==10|df$classe ==1|df$classe ==2|df$classe ==3|df$classe ==4|df$classe ==5|df$classe ==6|df$classe ==12)){df$B_type <- 'One-Sided Unidirectional'} 
      if(df$LNCO == 0 & df$groupe==2 & (df$classe ==30|df$classe ==31|df$classe ==32)){df$B_type <- 'Two-Sided Unidirectional'}         
      if(df$LNCO == 0 & df$groupe==3 & df$classe ==21){df$B_type <- 'Sym 4:4'} 
      if(df$LNCO == 0 & df$groupe==3 & df$classe ==20){df$B_type <- 'One-Sided Bidirectional'} 
      if(grepl("3:5a|5:3a", df$type)) {switch <- 'on'}
      if(df$LNCO == 0 & df$groupe==3 & df$classe ==20 & switch == 'on'){df$B_type <- 'Two-Sided Bidirectional'} 
      switch<- 'off'
      if(df$LNCO > 0 & df$groupe==0){df$B_type <- '2-sis NCO'}
      if(df$LCO == 1 & df$groupe==0){df$B_type <- '>2 Chr CO'}
      if(df$LCO > 1 & df$groupe==0){df$B_type <- 'dCO'}
      masterB <- rbind(masterB, df) }
  if(is.na(df$LCO) ==TRUE) { #if the event is a U event
      df$B_type <-  'Unknown'
      masterB <- rbind(masterB, df)}
}

#############################################################
# Calculate Spo11 hits local per event here and add column
#Loop through events and calculate Pan HpM for each event, then add to table b6$PanHpM. Looping is slow, but I don't have a better solution.

masterB$extend <- masterB$len_max/2

masterB$LocalEventHpM_PAN <- 0
masterB$LocalEventHpM_NWT <- 0

#################without Flanking####################
for (i in 1:nrow(masterB)){
  Pan.1=subset(Pan, Chr==masterB[i,"chr"] & Pos>=masterB[i,"start5"] & Pos <=masterB[i,"stop3"])
  masterB[i,"LocalEventHpM_PAN"]=sum(Pan.1$TotalHpM)
  NWT.1=subset(NeemanWT, Chr==masterB[i,"chr"] & Pos>=masterB[i,"start5"] & Pos <=masterB[i,"stop3"])
  masterB[i,"LocalEventHpM_NWT"]=sum(NWT.1$TotalHpM)
}
############################################

masterB$LocalEventHpM_PAN=round(masterB$LocalEventHpM_PAN,4) # Round PanHpM output to 2 decimals
masterB$LocalEventHpM_NWT=round(masterB$LocalEventHpM_NWT,4) 

#Write out master files
wd = getwd()
out = paste(wd,"/","masterAEventTable",sep="")
write.table(masterB, out, col.names =TRUE, row.names = FALSE, quote = FALSE, sep="\t", append=TRUE)

