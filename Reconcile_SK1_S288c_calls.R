#Reconcile the variant calls against the S288c and SK1 genomes.
#Change filenames and working directory as appropriate.
##Update 25-5-17 - now reports which spores were affected by each segment change.

library(plyr) 

##Open the files made by the Variant calling scripts.
####IMPORTANT NOTE: Only use the 'Converted HPs' file version when analysing msh2 tetrads. This is for calling mixed reads.#########
#Hash out the appropriate lines to open the correct files.
dfc <- read.table("OM1_combined_S288c.txt", header=T)
#dfc <- read.table("OM1_combined_ConvertedHPs_S288c.txt", header=T)
dfk <- read.table("OM1_combined_SK1.txt", header=T)
#dfk <- read.table("OM1_combined_ConvertedHPs_SK1.txt", header=T)

dfkc <- merge(dfk, dfc, by=c("uID", "sID", "chr", "pos_c", "pos_k", "type_k", "type_c", "Var_len"), sort=F, all=T)
dfkc <- dfkc[order(dfkc$chr, dfkc$pos_c),]
#it's important to sort after rather than during the merge, because otherwise it gets sorted on pos_k and messes up the order.

#there are rarely some rows in one set but not the other, they are kept if all=T but will contain NA
#turn NAs to 0
dfkc[is.na(dfkc)] <- 0
#this will give some warnings, because it doesn't want to convert the b columns to a number when they previously contained a string. 
#It's OK for them to remain NA, as long as the c and k columns are converted to 0.
dfkc[is.na(dfkc)] <- 'U' #Now convert the string columns to U value

dfkc$c1 <- (dfkc$c1.x+dfkc$c1.y)
dfkc$c2 <- (dfkc$c2.x+dfkc$c2.y)
dfkc$c3 <- (dfkc$c3.x+dfkc$c3.y)
dfkc$c4 <- (dfkc$c4.x+dfkc$c4.y)
dfkc$c5 <- (dfkc$c5.x+dfkc$c5.y)
dfkc$c6 <- (dfkc$c6.x+dfkc$c6.y)
dfkc$c7 <- (dfkc$c7.x+dfkc$c7.y)
dfkc$c8 <- (dfkc$c8.x+dfkc$c8.y)

dfkc$k1 <- (dfkc$k1.x+dfkc$k1.y)
dfkc$k2 <- (dfkc$k2.x+dfkc$k2.y)
dfkc$k3 <- (dfkc$k3.x+dfkc$k3.y)
dfkc$k4 <- (dfkc$k4.x+dfkc$k4.y)
dfkc$k5 <- (dfkc$k5.x+dfkc$k5.y)
dfkc$k6 <- (dfkc$k6.x+dfkc$k6.y)
dfkc$k7 <- (dfkc$k7.x+dfkc$k7.y)
dfkc$k8 <- (dfkc$k8.x+dfkc$k8.y)
#adding up k and c reads from the two alignments.
#These were previously divided by 2 to give the average, but that was removed to show the seperation better when plotting on a log scale.

####Special Rules only for Reconstructed Octads (TRM6 and TCMM15)
#Do not run this section for normal samples
#This is for comparing mixed reads with reads from one haplotype

#for(x in 1:nrow(dfkc)){
#if(dfkc[x,'b2'] =='CP' && (dfkc[x,'b1'] =='NP'|dfkc[x,'b1'] =='U')){dfkc[x,'b1'] <- 'CP'}
#if(dfkc[x,'b2'] =='KP' && (dfkc[x,'b1'] =='NP'|dfkc[x,'b1'] =='U')){dfkc[x,'b1'] <- 'KP'}
#if(dfkc[x,'kb2'] =='CP' && (dfkc[x,'kb1'] =='NP'|dfkc[x,'kb1'] =='U')){dfkc[x,'kb1'] <- 'CP'}
#if(dfkc[x,'kb2'] =='KP' && (dfkc[x,'kb1'] =='NP'|dfkc[x,'kb1'] =='U')){dfkc[x,'kb1'] <- 'KP'}
#if(dfkc[x,'b2'] =='CP' && dfkc[x,'b1'] =='HP'){dfkc[x,'b1'] <- 'KP'}
#if(dfkc[x,'b2'] =='KP' && dfkc[x,'b1'] =='HP'){dfkc[x,'b1'] <- 'CP'}
#if(dfkc[x,'kb2'] =='CP' && dfkc[x,'kb1'] =='HP'){dfkc[x,'kb1'] <- 'KP'}
#if(dfkc[x,'kb2'] =='KP' && dfkc[x,'kb1'] =='HP'){dfkc[x,'kb1'] <- 'CP'}

#if(dfkc[x,'b4'] =='CP' && (dfkc[x,'b3'] =='NP'|dfkc[x,'b3'] =='U')){dfkc[x,'b3'] <- 'CP'}
#if(dfkc[x,'b4'] =='KP' && (dfkc[x,'b3'] =='NP'|dfkc[x,'b3'] =='U')){dfkc[x,'b3'] <- 'KP'}
#if(dfkc[x,'kb4'] =='CP' && (dfkc[x,'kb3'] =='NP'|dfkc[x,'kb3'] =='U')){dfkc[x,'kb3'] <- 'CP'}
#if(dfkc[x,'kb4'] =='KP' && (dfkc[x,'kb3'] =='NP'|dfkc[x,'kb3'] =='U')){dfkc[x,'kb3'] <- 'KP'}
#if(dfkc[x,'b4'] =='CP' && dfkc[x,'b3'] =='HP'){dfkc[x,'b3'] <- 'KP'}
#if(dfkc[x,'b4'] =='KP' && dfkc[x,'b3'] =='HP'){dfkc[x,'b3'] <- 'CP'}
#if(dfkc[x,'kb4'] =='CP' && dfkc[x,'kb3'] =='HP'){dfkc[x,'kb3'] <- 'KP'}
#if(dfkc[x,'kb4'] =='KP' && dfkc[x,'kb3'] =='HP'){dfkc[x,'kb3'] <- 'CP'}

#if(dfkc[x,'b6'] =='CP' && (dfkc[x,'b5'] =='NP'|dfkc[x,'b5'] =='U')){dfkc[x,'b5'] <- 'CP'}
#if(dfkc[x,'b6'] =='KP' && (dfkc[x,'b5'] =='NP'|dfkc[x,'b5'] =='U')){dfkc[x,'b5'] <- 'KP'}
#if(dfkc[x,'kb6'] =='CP' && (dfkc[x,'kb5'] =='NP'|dfkc[x,'kb5'] =='U')){dfkc[x,'kb5'] <- 'CP'}
#if(dfkc[x,'kb6'] =='KP' && (dfkc[x,'kb5'] =='NP'|dfkc[x,'kb5'] =='U')){dfkc[x,'kb5'] <- 'KP'}
#if(dfkc[x,'b6'] =='CP' && dfkc[x,'b5'] =='HP'){dfkc[x,'b5'] <- 'KP'}
#if(dfkc[x,'b6'] =='KP' && dfkc[x,'b5'] =='HP'){dfkc[x,'b5'] <- 'CP'}
#if(dfkc[x,'kb6'] =='CP' && dfkc[x,'kb5'] =='HP'){dfkc[x,'kb5'] <- 'KP'}
#if(dfkc[x,'kb6'] =='KP' && dfkc[x,'kb5'] =='HP'){dfkc[x,'kb5'] <- 'CP'}

#if(dfkc[x,'b8'] =='CP' && (dfkc[x,'b7'] =='NP'|dfkc[x,'b7'] =='U')){dfkc[x,'b7'] <- 'CP'}
#if(dfkc[x,'b8'] =='KP' && (dfkc[x,'b7'] =='NP'|dfkc[x,'b7'] =='U')){dfkc[x,'b7'] <- 'KP'}
#if(dfkc[x,'kb8'] =='CP' && (dfkc[x,'kb7'] =='NP'|dfkc[x,'kb7'] =='U')){dfkc[x,'kb7'] <- 'CP'}
#if(dfkc[x,'kb8'] =='KP' && (dfkc[x,'kb7'] =='NP'|dfkc[x,'kb7'] =='U')){dfkc[x,'kb7'] <- 'KP'}
#if(dfkc[x,'b8'] =='CP' && dfkc[x,'b7'] =='HP'){dfkc[x,'b7'] <- 'KP'}
#if(dfkc[x,'b8'] =='KP' && dfkc[x,'b7'] =='HP'){dfkc[x,'b7'] <- 'CP'}
#if(dfkc[x,'kb8'] =='CP' && dfkc[x,'kb7'] =='HP'){dfkc[x,'kb7'] <- 'KP'}
#if(dfkc[x,'kb8'] =='KP' && dfkc[x,'kb7'] =='HP'){dfkc[x,'kb7'] <- 'CP'}
#}

######################################################
##resolve conflicts between the SK1 and S288c calls###
dfkc$b1<- ifelse(dfkc$kb1 == "KP" & dfkc$b1=='CP', NA, ifelse(dfkc$kb1 == "CP" & dfkc$b1=='KP', NA, ifelse(dfkc$kb1 == "KP", 0, ifelse(dfkc$kb1 == "CP", 1, ifelse(dfkc$b1 == "KP", 0, ifelse(dfkc$b1 == "CP", 1, NA))))))
dfkc$b2<- ifelse(dfkc$kb2 == "KP" & dfkc$b2=='CP', NA, ifelse(dfkc$kb2 == "CP" & dfkc$b2=='KP', NA, ifelse(dfkc$kb2 == "KP", 0, ifelse(dfkc$kb2 == "CP", 1, ifelse(dfkc$b2 == "KP", 0, ifelse(dfkc$b2 == "CP", 1, NA))))))
dfkc$b3<- ifelse(dfkc$kb3 == "KP" & dfkc$b3=='CP', NA, ifelse(dfkc$kb3 == "CP" & dfkc$b3=='KP', NA, ifelse(dfkc$kb3 == "KP", 0, ifelse(dfkc$kb3 == "CP", 1, ifelse(dfkc$b3 == "KP", 0, ifelse(dfkc$b3 == "CP", 1, NA))))))
dfkc$b4<- ifelse(dfkc$kb4 == "KP" & dfkc$b4=='CP', NA, ifelse(dfkc$kb4 == "CP" & dfkc$b4=='KP', NA, ifelse(dfkc$kb4 == "KP", 0, ifelse(dfkc$kb4 == "CP", 1, ifelse(dfkc$b4 == "KP", 0, ifelse(dfkc$b4 == "CP", 1, NA))))))
dfkc$b5<- ifelse(dfkc$kb5 == "KP" & dfkc$b5=='CP', NA, ifelse(dfkc$kb5 == "CP" & dfkc$b5=='KP', NA, ifelse(dfkc$kb5 == "KP", 0, ifelse(dfkc$kb5 == "CP", 1, ifelse(dfkc$b5 == "KP", 0, ifelse(dfkc$b5 == "CP", 1, NA))))))
dfkc$b6<- ifelse(dfkc$kb6 == "KP" & dfkc$b6=='CP', NA, ifelse(dfkc$kb6 == "CP" & dfkc$b6=='KP', NA, ifelse(dfkc$kb6 == "KP", 0, ifelse(dfkc$kb6 == "CP", 1, ifelse(dfkc$b6 == "KP", 0, ifelse(dfkc$b6 == "CP", 1, NA))))))
dfkc$b7<- ifelse(dfkc$kb7 == "KP" & dfkc$b7=='CP', NA, ifelse(dfkc$kb7 == "CP" & dfkc$b7=='KP', NA, ifelse(dfkc$kb7 == "KP", 0, ifelse(dfkc$kb7 == "CP", 1, ifelse(dfkc$b7 == "KP", 0, ifelse(dfkc$b7 == "CP", 1, NA))))))
dfkc$b8<- ifelse(dfkc$kb8 == "KP" & dfkc$b8=='CP', NA, ifelse(dfkc$kb8 == "CP" & dfkc$b8=='KP', NA, ifelse(dfkc$kb8 == "KP", 0, ifelse(dfkc$kb8 == "CP", 1, ifelse(dfkc$b8 == "KP", 0, ifelse(dfkc$b8 == "CP", 1, NA))))))

#If the two files conflict, value is NA e.g. if kb1=="KP" & b1=="CP".
#value is also NA if neither call produced a KP or CP
#c=1 k=0

####
#remove some of the unnecessary columns
dfkc<- subset(dfkc, select = c(chr, pos_c, c1,c2,c3,c4,c5,c6,c7,c8,k1,k2,k3,k4,k5,k6,k7,k8,b1,b2,b3,b4,b5,b6,b7,b8, type_k, Var_len))

#Some calculations for quality control##

##work out how many rows contain at least one NA, which will be discarded later##
statistic5 <- c('Total_Rows', 'Removed_Rows', 'Percent_Removed_Rows', 'S_total', 'I_total', 'D_total', 'Removed_S', 'Percent_Removed_S', 'Removed_I', 'Percent_Removed_I', 'Removed_D', 'Percent_Removed_D')
value5 <- c(0,0,0,0,0,0,0,0,0,0,0,0)
statstable5 <- data.frame(statistic5, value5)

#Work out some stats about the number of rows containing at least one NA 
dfkc_noNA <- dfkc[complete.cases(dfkc),]
statstable5[1,2] <- nrow(dfkc)
statstable5[2,2] <- (nrow(dfkc) - nrow(dfkc_noNA))
statstable5[3,2] <- (statstable5[2,2]/statstable5[1,2])*100
s <- dfkc[which(dfkc$type_k=='s'),]
i <- dfkc[which(dfkc$type_k=='i'),]
d <- dfkc[which(dfkc$type_k=='d'),]
s_noNA <- dfkc_noNA[which(dfkc_noNA$type_k=='s'),]
i_noNA <- dfkc_noNA[which(dfkc_noNA$type_k=='i'),]
d_noNA <- dfkc_noNA[which(dfkc_noNA$type_k=='d'),]
statstable5[4,2]  <- nrow(s)
statstable5[5,2]  <- nrow(i)
statstable5[6,2]  <- nrow(d)
statstable5[7,2]  <- nrow(s)-nrow(s_noNA)
statstable5[8,2]  <- ((statstable5[7,2]/statstable5[4,2])*100)
statstable5[9,2]  <- nrow(i)-nrow(i_noNA)
statstable5[10,2]  <- ((statstable5[9,2]/statstable5[5,2])*100)
statstable5[11,2]  <- nrow(d)-nrow(d_noNA)
statstable5[12,2]  <- ((statstable5[11,2]/statstable5[6,2])*100)
statstable5$value5 <- round(statstable5$value5, digits = 2)   

write.table(statstable5,file="OM1_additional_Stats_S288c_SK1.txt", row.names=F, quote=F, sep='\t')

################################################
#makes the gr2 file which can be visualized in R

####Make gr2 with only INDEL 8-0 positions removed (keep SNP 8-0s)####
gr2 <- dfkc
gr2$temp <- (gr2$b1 + gr2$b2 + gr2$b3 + gr2$b4 + gr2$b5 + gr2$b6 + gr2$b7 + gr2$b8)

gr2i <- gr2[which(gr2$type_k=='i' & gr2$temp!=0 & gr2$temp!=8),] #This also removes ones where the total is NA.
gr2d <- gr2[which(gr2$type_k=='d' & gr2$temp!=0 & gr2$temp!=8),]
gr2s <- gr2[which(gr2$type_k=='s'),]
gr2NA <- subset(gr2, is.na(gr2$temp)) #subset the ones where the total is NA
gr2 <- rbind(gr2i, gr2d, gr2s, gr2NA)
gr2 <- gr2[order(gr2$chr, gr2$pos_c),]
#remove the temp column 
gr2<- subset(gr2, select = c(chr, pos_c, c1,c2,c3,c4,c5,c6,c7,c8,k1,k2,k3,k4,k5,k6,k7,k8,b1,b2,b3,b4,b5,b6,b7,b8, type_k, Var_len))
gr2 <- unique(gr2) #remove any duplicates that may have been introduced from SNPs totalling NA

write.table(gr2,file="OM1_gr2.txt", row.names=F, quote=F, sep='\t') #if you want to keep 8-0 SNPs, use this file for visualization

#This section makes the file compatible with the event caller
sub <- subset(dfkc2, select = c(b1,b2,b3,b4,b5,b6,b7,b8))
#isolates the b columns

sub[!!rowSums(is.na(sub)),] <- NA
#if any column contains an NA, then all columns for that row are converted to NA

dfkc2$b1 <- sub$b1
dfkc2$b2 <- sub$b2
dfkc2$b3 <- sub$b3
dfkc2$b4 <- sub$b4
dfkc2$b5 <- sub$b5
dfkc2$b6 <- sub$b6
dfkc2$b7 <- sub$b7
dfkc2$b8 <- sub$b8
#Puts the b columns back into the dataframe

write.table(dfkc2,file="OM1_binarized_CK.txt", row.names=F, quote=F, sep='\t')
df <- read.table("OM1_binarized_CK.txt", header=T)

df$temp <- (df$b1 + df$b2 + df$b3 + df$b4 + df$b5 + df$b6 + df$b7 + df$b8)

#This outputs a table containing 8-0 positions for examination as they may not actually be SNPs
nonsnps8 <-df[which(df$temp==8),]
nonsnps0 <-df[which(df$temp==0),]
total <- rbind(nonsnps0, nonsnps8)

write.table(total,file="EZ-list_OM1_KC.txt", row.names=F, quote=F, sep='\t')

####Make binary with only INDEL 8-0 positions removed (keep SNP 8-0s)####

dfi <- df[which(df$type_k=='i' & df$temp!=0 & df$temp!=8),]
dfd <- df[which(df$type_k=='d' & df$temp!=0 & df$temp!=8),]
dfs <- df[which(df$type_k=='s'),]
df3 <- rbind(dfi, dfd, dfs)
df3 <- df3[order(df3$chr, df3$pos_c),]

statistic6 <- c('Total_Rows', 'Indel 8-0s Removed')
value6 <- c(0,0)
statstable6 <- data.frame(statistic6, value6)
statstable6[1,2] <- nrow(df)
statstable6[2,2] <- (nrow(df) - nrow(df3))
write.table(statstable6,file="OM1_Indel_8-0s_Removed.txt", row.names=F, quote=F, sep='\t')

#remove the temp column and the c and k columns
df3<- subset(df3, select = c(chr, pos_c, b1,b2,b3,b4,b5,b6,b7,b8))
#To remove 'NA's, select only rows containing 0 or 1 
df3 <- df3[df3$b1 %in% c(0, 1), ]

############################################################
##detail which spores were affected by each pattern change##
df3$sporeaffect <- '0'
for(i in 1:nrow(df3)){
  sub1 <- df3[i,]
  sub2 <- df3[(i-1),]
  string <- '0'
  if(nrow(sub2)>0){
    if(sub1$chr==sub2$chr){
      if(sub1$b1 != sub2$b1){ string <- paste(string, '1', sep = "")}
      if(sub1$b2 != sub2$b2){ string <- paste(string, '2', sep = "")}
      if(sub1$b3 != sub2$b3){ string <- paste(string, '3', sep = "")}
      if(sub1$b4 != sub2$b4){ string <- paste(string, '4', sep = "")}
      if(sub1$b5 != sub2$b5){ string <- paste(string, '5', sep = "")}
      if(sub1$b6 != sub2$b6){ string <- paste(string, '6', sep = "")}
      if(sub1$b7 != sub2$b7){ string <- paste(string, '7', sep = "")}
      if(sub1$b8 != sub2$b8){ string <- paste(string, '8', sep = "")}
      df3[(i-1),11] <- string
      }
   }
  
}

write.table(df3,file="OM1_binary.txt", row.names=F, quote=F, sep='\t')
#This file is compatible with the event caller.

##Some more QC calculations###

dfk$b1<- ifelse(dfk$kb1 == "KP", 0, ifelse(dfk$kb1 == "CP", 1, NA))
dfk$b2<- ifelse(dfk$kb2 == "KP", 0, ifelse(dfk$kb2 == "CP", 1, NA))
dfk$b3<- ifelse(dfk$kb3 == "KP", 0, ifelse(dfk$kb3 == "CP", 1, NA))
dfk$b4<- ifelse(dfk$kb4 == "KP", 0, ifelse(dfk$kb4 == "CP", 1, NA))
dfk$b5<- ifelse(dfk$kb5 == "KP", 0, ifelse(dfk$kb5 == "CP", 1, NA))
dfk$b6<- ifelse(dfk$kb6 == "KP", 0, ifelse(dfk$kb6 == "CP", 1, NA))
dfk$b7<- ifelse(dfk$kb7 == "KP", 0, ifelse(dfk$kb7 == "CP", 1, NA))
dfk$b8<- ifelse(dfk$kb8 == "KP", 0, ifelse(dfk$kb8 == "CP", 1, NA))

dfk<- subset(dfk, select = c(chr, pos_c, type_k, c1,c2,c3,c4,c5,c6,c7,c8,k1,k2,k3,k4,k5,k6,k7,k8,b1,b2,b3,b4,b5,b6,b7,b8))

##work out how many rows contain at least one NA, and so will be discarded later##
statistic5 <- c('Total_Rows', 'Removed_Rows', 'Percent_Removed_Rows', 'S_total', 'I_total', 'D_total', 'Removed_S', 'Percent_Removed_S', 'Removed_I', 'Percent_Removed_I', 'Removed_D', 'Percent_Removed_D')
value5 <- c(0,0,0,0,0,0,0,0,0,0,0,0)
statstable5 <- data.frame(statistic5, value5)

#Work out some stats about the number of rows containing at least one NA (which will be discarded)
df_noNA <- dfk[complete.cases(dfk),]
statstable5[1,2] <- nrow(dfk)
statstable5[2,2] <- (nrow(dfk) - nrow(df_noNA))
statstable5[3,2] <- (statstable5[2,2]/statstable5[1,2])*100
s <- dfk[which(dfk$type_k=='s'),]
i <- dfk[which(dfk$type_k=='i'),]
d <- dfk[which(dfk$type_k=='d'),]
s_noNA <- df_noNA[which(df_noNA$type_k=='s'),]
i_noNA <- df_noNA[which(df_noNA$type_k=='i'),]
d_noNA <- df_noNA[which(df_noNA$type_k=='d'),]
statstable5[4,2]  <- nrow(s)
statstable5[5,2]  <- nrow(i)
statstable5[6,2]  <- nrow(d)
statstable5[7,2]  <- nrow(s)-nrow(s_noNA)
statstable5[8,2]  <- ((statstable5[7,2]/statstable5[4,2])*100)
statstable5[9,2]  <- nrow(i)-nrow(i_noNA)
statstable5[10,2]  <- ((statstable5[9,2]/statstable5[5,2])*100)
statstable5[11,2]  <- nrow(d)-nrow(d_noNA)
statstable5[12,2]  <- ((statstable5[11,2]/statstable5[6,2])*100)
statstable5$value5 <- round(statstable5$value5, digits = 2)   

write.table(statstable5,file="OM1_additional_Stats_SK1.txt", row.names=F, quote=F, sep='\t')

dfc$b1<- ifelse(dfc$b1 == "KP", 0, ifelse(dfc$b1 == "CP", 1, NA))
dfc$b2<- ifelse(dfc$b2 == "KP", 0, ifelse(dfc$b2 == "CP", 1, NA))
dfc$b3<- ifelse(dfc$b3 == "KP", 0, ifelse(dfc$b3 == "CP", 1, NA))
dfc$b4<- ifelse(dfc$b4 == "KP", 0, ifelse(dfc$b4 == "CP", 1, NA))
dfc$b5<- ifelse(dfc$b5 == "KP", 0, ifelse(dfc$b5 == "CP", 1, NA))
dfc$b6<- ifelse(dfc$b6 == "KP", 0, ifelse(dfc$b6 == "CP", 1, NA))
dfc$b7<- ifelse(dfc$b7 == "KP", 0, ifelse(dfc$b7 == "CP", 1, NA))
dfc$b8<- ifelse(dfc$b8 == "KP", 0, ifelse(dfc$b8 == "CP", 1, NA))

dfc<- subset(dfc, select = c(chr, pos_c, type_k, c1,c2,c3,c4,c5,c6,c7,c8,k1,k2,k3,k4,k5,k6,k7,k8,b1,b2,b3,b4,b5,b6,b7,b8))

##work out how many rows contain at least one NA, and so will be discarded later##
statistic5 <- c('Total_Rows', 'Removed_Rows', 'Percent_Removed_Rows', 'S_total', 'I_total', 'D_total', 'Removed_S', 'Percent_Removed_S', 'Removed_I', 'Percent_Removed_I', 'Removed_D', 'Percent_Removed_D')
value5 <- c(0,0,0,0,0,0,0,0,0,0,0,0)
statstable5 <- data.frame(statistic5, value5)

#Work out some stats about the number of rows containing at least one NA (which will be discarded)
df_noNA <- dfc[complete.cases(dfc),]
statstable5[1,2] <- nrow(dfc)
statstable5[2,2] <- (nrow(dfc) - nrow(df_noNA))
statstable5[3,2] <- (statstable5[2,2]/statstable5[1,2])*100
s <- dfc[which(dfc$type_k=='s'),]
i <- dfc[which(dfc$type_k=='i'),]
d <- dfc[which(dfc$type_k=='d'),]
s_noNA <- df_noNA[which(df_noNA$type_k=='s'),]
i_noNA <- df_noNA[which(df_noNA$type_k=='i'),]
d_noNA <- df_noNA[which(df_noNA$type_k=='d'),]
statstable5[4,2]  <- nrow(s)
statstable5[5,2]  <- nrow(i)
statstable5[6,2]  <- nrow(d)
statstable5[7,2]  <- nrow(s)-nrow(s_noNA)
statstable5[8,2]  <- ((statstable5[7,2]/statstable5[4,2])*100)
statstable5[9,2]  <- nrow(i)-nrow(i_noNA)
statstable5[10,2]  <- ((statstable5[9,2]/statstable5[5,2])*100)
statstable5[11,2]  <- nrow(d)-nrow(d_noNA)
statstable5[12,2]  <- ((statstable5[11,2]/statstable5[6,2])*100)
statstable5$value5 <- round(statstable5$value5, digits = 2)   

write.table(statstable5,file="OM1_additional_Stats_S288c.txt", row.names=F, quote=F, sep='\t')
