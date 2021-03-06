#Run the variant caller for the SK1 alignment.
#Set working directory and change filenames as required.
#If you are analysing tetrads, simply hash out the file opening and function calls for every other sample.

#load plyr package
library(plyr)

#Source the function script
source("SNP_indel_caller_S288c.R")

#opening the files created by SNPisolator
df1 <- read.table("OM1AOut_SK1.txt",header=T, stringsAsFactors=FALSE)
df2 <- read.table("OM1AlphaOut_SK1.txt",header=T, stringsAsFactors=FALSE)
df3 <- read.table("OM1BOut_SK1.txt",header=T, stringsAsFactors=FALSE)
df4 <- read.table("OM1BetaOut_SK1.txt",header=T, stringsAsFactors=FALSE)
df5 <- read.table("OM1COut_SK1.txt",header=T, stringsAsFactors=FALSE)
df6 <- read.table("OM1GammaOut_SK1.txt",header=T, stringsAsFactors=FALSE)
df7 <- read.table("OM1DOut_SK1.txt",header=T, stringsAsFactors=FALSE)
df8 <- read.table("OM1DeltaOut_SK1.txt",header=T, stringsAsFactors=FALSE)

#Running the function
SNP_indel_caller_SK1(df1, "OM1A", 5, 5, 10)
SNP_indel_caller_SK1(df2, "OM1Alpha", 5,5,10)
SNP_indel_caller_SK1(df3, "OM1B", 5,5,10)
SNP_indel_caller_SK1(df4, "OM1Beta", 5,5,10)
SNP_indel_caller_SK1(df5, "OM1C", 5,5,10)
SNP_indel_caller_SK1(df6, "OM1Gamma", 5,5,10)
SNP_indel_caller_SK1(df7, "OM1D", 5,5,10)
SNP_indel_caller_SK1(df8, "OM1Delta", 5,5,10)

#df - df to use
#filename - name for files produced
#threshold_i -minimum number of reads for a SNP to be called
#threshold_s -minimum number of reads for an indel to be called
#threshold_k - additional threshold for indels - if SK1-type reads are above this limit, they have a more lenient calling threshold

###################################################################################

###################################################################################
#Combine individual stats files for each octad
###SNPS####
stat1 <- read.table("OM1ASNP_Stats_SK1.txt",header=T)
stat2 <- read.table("OM1AlphaSNP_Stats_SK1.txt",header=T)
stat3 <- read.table("OM1BSNP_Stats_SK1.txt",header=T)
stat4 <- read.table("OM1BetaSNP_Stats_SK1.txt",header=T)
stat5 <- read.table("OM1CSNP_Stats_SK1.txt",header=T)
stat6 <- read.table("OM1GammaSNP_Stats_SK1.txt",header=T)
stat7 <- read.table("OM1DSNP_Stats_SK1.txt",header=T)
stat8 <- read.table("OM1DeltaSNP_Stats_SK1.txt",header=T)

statmerged <- merge(stat1, stat2, by=c("statistic"), sort=FALSE)
statmerged <- merge(statmerged, stat3, by=c("statistic"), sort=FALSE)
statmerged <- merge(statmerged, stat4, by=c("statistic"), sort=FALSE)
statmerged <- merge(statmerged, stat5, by=c("statistic"), sort=FALSE)
statmerged <- merge(statmerged, stat6, by=c("statistic"), sort=FALSE)
statmerged <- merge(statmerged, stat7, by=c("statistic"), sort=FALSE)
statmerged <- merge(statmerged, stat8, by=c("statistic"), sort=FALSE)

statmerged$Octad_averages <- rowMeans(statmerged[,2:9])
statmerged$Octad_averages <- round(statmerged$Octad_averages, digits = 3)

write.table(statmerged,file="OM1SNP_Stats_SK1.txt",row.names=F, quote=F, sep='\t')

file.remove("OM1ASNP_Stats_SK1.txt")
file.remove("OM1AlphaSNP_Stats_SK1.txt")
file.remove("OM1BSNP_Stats_SK1.txt")
file.remove("OM1BetaSNP_Stats_SK1.txt")
file.remove("OM1CSNP_Stats_SK1.txt")
file.remove("OM1GammaSNP_Stats_SK1.txt")
file.remove("OM1DSNP_Stats_SK1.txt")
file.remove("OM1DeltaSNP_Stats_SK1.txt")

###insertions###
stat1 <- read.table("OM1Ainsertion_Stats_SK1.txt",header=T)
stat2 <- read.table("OM1Alphainsertion_Stats_SK1.txt",header=T)
stat3 <- read.table("OM1Binsertion_Stats_SK1.txt",header=T)
stat4 <- read.table("OM1Betainsertion_Stats_SK1.txt",header=T)
stat5 <- read.table("OM1Cinsertion_Stats_SK1.txt",header=T)
stat6 <- read.table("OM1Gammainsertion_Stats_SK1.txt",header=T)
stat7 <- read.table("OM1Dinsertion_Stats_SK1.txt",header=T)
stat8 <- read.table("OM1Deltainsertion_Stats_SK1.txt",header=T)

statmerged <- merge(stat1, stat2, by=c("statistic2"), sort=FALSE)
statmerged <- merge(statmerged, stat3, by=c("statistic2"), sort=FALSE)
statmerged <- merge(statmerged, stat4, by=c("statistic2"), sort=FALSE)
statmerged <- merge(statmerged, stat5, by=c("statistic2"), sort=FALSE)
statmerged <- merge(statmerged, stat6, by=c("statistic2"), sort=FALSE)
statmerged <- merge(statmerged, stat7, by=c("statistic2"), sort=FALSE)
statmerged <- merge(statmerged, stat8, by=c("statistic2"), sort=FALSE)

statmerged$Octad_averages <- rowMeans(statmerged[,2:9])
statmerged$Octad_averages <- round(statmerged$Octad_averages, digits = 3)                             

write.table(statmerged,file="OM1insertion_Stats_SK1.txt",row.names=F, quote=F, sep='\t')

file.remove("OM1Ainsertion_Stats_SK1.txt")
file.remove("OM1Alphainsertion_Stats_SK1.txt")
file.remove("OM1Binsertion_Stats_SK1.txt")
file.remove("OM1Betainsertion_Stats_SK1.txt")
file.remove("OM1Cinsertion_Stats_SK1.txt")
file.remove("OM1Gammainsertion_Stats_SK1.txt")
file.remove("OM1Dinsertion_Stats_SK1.txt")
file.remove("OM1Deltainsertion_Stats_SK1.txt")

###deletions###

stat1 <- read.table("OM1Adeletion_Stats_after_SK1.txt",header=T)
stat2 <- read.table("OM1Alphadeletion_Stats_after_SK1.txt",header=T)
stat3 <- read.table("OM1Bdeletion_Stats_after_SK1.txt",header=T)
stat4 <- read.table("OM1Betadeletion_Stats_after_SK1.txt",header=T)
stat5 <- read.table("OM1Cdeletion_Stats_after_SK1.txt",header=T)
stat6 <- read.table("OM1Gammadeletion_Stats_after_SK1.txt",header=T)
stat7 <- read.table("OM1Ddeletion_Stats_after_SK1.txt",header=T)
stat8 <- read.table("OM1Deltadeletion_Stats_after_SK1.txt",header=T)

statmerged <- merge(stat1, stat2, by=c("statistic3"), sort=FALSE)
statmerged <- merge(statmerged, stat3, by=c("statistic3"), sort=FALSE)
statmerged <- merge(statmerged, stat4, by=c("statistic3"), sort=FALSE)
statmerged <- merge(statmerged, stat5, by=c("statistic3"), sort=FALSE)
statmerged <- merge(statmerged, stat6, by=c("statistic3"), sort=FALSE)
statmerged <- merge(statmerged, stat7, by=c("statistic3"), sort=FALSE)
statmerged <- merge(statmerged, stat8, by=c("statistic3"), sort=FALSE)

statmerged$Octad_averages <- rowMeans(statmerged[,2:9])
statmerged$Octad_averages <- round(statmerged$Octad_averages, digits = 3)                             

write.table(statmerged,file="OM1deletion_Stats_after_SK1.txt",row.names=F, quote=F, sep='\t')

file.remove("OM1Adeletion_Stats_after_SK1.txt")
file.remove("OM1Alphadeletion_Stats_after_SK1.txt")
file.remove("OM1Bdeletion_Stats_after_SK1.txt")
file.remove("OM1Betadeletion_Stats_after_SK1.txt")
file.remove("OM1Cdeletion_Stats_after_SK1.txt")
file.remove("OM1Gammadeletion_Stats_after_SK1.txt")
file.remove("OM1Ddeletion_Stats_after_SK1.txt")
file.remove("OM1Deltadeletion_Stats_after_SK1.txt")

stat1 <- read.table("OM1Adeletion_Stats_before_SK1.txt",header=T)
stat2 <- read.table("OM1Alphadeletion_Stats_before_SK1.txt",header=T)
stat3 <- read.table("OM1Bdeletion_Stats_before_SK1.txt",header=T)
stat4 <- read.table("OM1Betadeletion_Stats_before_SK1.txt",header=T)
stat5 <- read.table("OM1Cdeletion_Stats_before_SK1.txt",header=T)
stat6 <- read.table("OM1Gammadeletion_Stats_before_SK1.txt",header=T)
stat7 <- read.table("OM1Ddeletion_Stats_before_SK1.txt",header=T)
stat8 <- read.table("OM1Deltadeletion_Stats_before_SK1.txt",header=T)

statmerged <- merge(stat1, stat2, by=c("statistic4"), sort=FALSE)
statmerged <- merge(statmerged, stat3, by=c("statistic4"), sort=FALSE)
statmerged <- merge(statmerged, stat4, by=c("statistic4"), sort=FALSE)
statmerged <- merge(statmerged, stat5, by=c("statistic4"), sort=FALSE)
statmerged <- merge(statmerged, stat6, by=c("statistic4"), sort=FALSE)
statmerged <- merge(statmerged, stat7, by=c("statistic4"), sort=FALSE)
statmerged <- merge(statmerged, stat8, by=c("statistic4"), sort=FALSE)

statmerged$Octad_averages <- rowMeans(statmerged[,2:9])
statmerged$Octad_averages <- round(statmerged$Octad_averages, digits = 3)                             

write.table(statmerged,file="OM1deletion_Stats_before_SK1.txt",row.names=F, quote=F, sep='\t')

#delete the individual files
file.remove("OM1Adeletion_Stats_before_SK1.txt")
file.remove("OM1Alphadeletion_Stats_before_SK1.txt")
file.remove("OM1Bdeletion_Stats_before_SK1.txt")
file.remove("OM1Betadeletion_Stats_before_SK1.txt")
file.remove("OM1Cdeletion_Stats_before_SK1.txt")
file.remove("OM1Gammadeletion_Stats_before_SK1.txt")
file.remove("OM1Ddeletion_Stats_before_SK1.txt")
file.remove("OM1Deltadeletion_Stats_before_SK1.txt")
###################################################################################

###################################################################################
#Called SNPs Combiner
#Combines the 8 members of an octad into 1 file
table1 <- read.table("OM1ACalled_SK1.txt",header=T)
table2 <- read.table("OM1AlphaCalled_SK1.txt",header=T)
table3 <- read.table("OM1BCalled_SK1.txt",header=T)
table4 <- read.table("OM1BetaCalled_SK1.txt",header=T)
table5 <- read.table("OM1CCalled_SK1.txt",header=T)
table6 <- read.table("OM1GammaCalled_SK1.txt",header=T)
table7 <- read.table("OM1DCalled_SK1.txt",header=T)
table8 <- read.table("OM1DeltaCalled_SK1.txt",header=T)

library(plyr)
#plyr is needed for the renaming function so all the c, k and b columns don't have the same name

table1 <- rename(table1, c("c"="c1", "k"="k1", "kb"="kb1"))
table2 <- rename(table2, c("c"="c2", "k"="k2", "kb"="kb2"))
table3 <- rename(table3, c("c"="c3", "k"="k3", "kb"="kb3"))
table4 <- rename(table4, c("c"="c4", "k"="k4", "kb"="kb4"))
table5 <- rename(table5, c("c"="c5", "k"="k5", "kb"="kb5"))
table6 <- rename(table6, c("c"="c6", "k"="k6", "kb"="kb6"))
table7 <- rename(table7, c("c"="c7", "k"="k7", "kb"="kb7"))
table8 <- rename(table8, c("c"="c8", "k"="k8", "kb"="kb8"))

#merges the 8 dataframes into 1
merged <- merge(table1, table2, by=c("uID", "sID",  "chrom", "pos_c", "pos_k", "type_k", "type_c", "Var_len"), all=FALSE, sort=FALSE)
merged <- merge(merged, table3, by=c("uID", "sID",  "chrom", "pos_c", "pos_k", "type_k", "type_c", "Var_len"), all=FALSE, sort=FALSE)
merged <- merge(merged, table4, by=c("uID", "sID",  "chrom", "pos_c", "pos_k", "type_k", "type_c", "Var_len"), all=FALSE, sort=FALSE)
merged <- merge(merged, table5, by=c("uID", "sID",  "chrom", "pos_c", "pos_k", "type_k", "type_c", "Var_len"), all=FALSE, sort=FALSE)
merged <- merge(merged, table6, by=c("uID", "sID",  "chrom", "pos_c", "pos_k", "type_k", "type_c", "Var_len"), all=FALSE, sort=FALSE)
merged <- merge(merged, table7, by=c("uID", "sID",  "chrom", "pos_c", "pos_k", "type_k", "type_c", "Var_len"), all=FALSE, sort=FALSE)
merged <- merge(merged, table8, by=c("uID", "sID",  "chrom", "pos_c", "pos_k", "type_k", "type_c", "Var_len"), all=FALSE, sort=FALSE)

#Puts the columns in the right order
merged <- merged[c("uID", "sID", "chrom", "pos_c", "pos_k", "type_c", "type_k", "Var_len", "c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8", 
                   "k1", "k2", "k3", "k4", "k5", "k6", "k7", "k8", "kb1", "kb2", "kb3", "kb4", "kb5", "kb6", "kb7", "kb8")]

#renames these two columns to match MC's specification
merged <- rename(merged, c("chrom" = "chr"))

write.table(merged,file="OM1_combined_SK1.txt", row.names=F, quote=F, sep='\t')

file.remove("OM1ACalled_SK1.txt")
file.remove("OM1AlphaCalled_SK1.txt")
file.remove("OM1BCalled_SK1.txt")
file.remove("OM1BetaCalled_SK1.txt")
file.remove("OM1CCalled_SK1.txt")
file.remove("OM1GammaCalled_SK1.txt")
file.remove("OM1DCalled_SK1.txt")
file.remove("OM1DeltaCalled_SK1.txt")

###################################################################################
#call heteroduplex in tetrads

#Only for msh2 tetrads and tetradified octads
#tetrad data is adapted to the octad pipeline by duplicating each spore.
#If a msh2 tetrad has a position called 'HP', one copy of the spore should have it converted to KP and the other to CP. 
#Thus 'octadifying' the tetrad and making it a 7:1 event.
#This would keep the output as if it were an artificial (and slightly incorrect) octad.
#Issue: Commonly, heteroduplex occurs near the end of chromosomes which is probably actually just from repetitive sequences. 
#This creates small 'U' events entirely composed of hDNA. These can be removed from the event table.

df <- read.table("OM1_combined_SK1.txt", header=T)
#only if the type==S (indels seem to never be called HP anyway, and dels get converted to NP earlier - but may as well be on the safe side)
df$kb1 <- as.character(df$kb1)
df$kb2 <- as.character(df$kb2)
df$kb3 <- as.character(df$kb3)
df$kb4 <- as.character(df$kb4)
df$kb5 <- as.character(df$kb5)
df$kb6 <- as.character(df$kb6)
df$kb7 <- as.character(df$kb7)
df$kb8 <- as.character(df$kb8)
#ifelse statements don't work properly on factors - they return the factor level - so make sure these columns contain characters
df$kb1 <- ifelse(df$type_k=='s' & df$kb1 == 'HP', 'KP', df$kb1)
df$kb2 <- ifelse(df$type_k=='s' & df$kb2 == 'HP', 'CP', df$kb2)
df$kb3 <- ifelse(df$type_k=='s' & df$kb3 == 'HP', 'KP', df$kb3)
df$kb4 <- ifelse(df$type_k=='s' & df$kb4 == 'HP', 'CP', df$kb4)
df$kb5 <- ifelse(df$type_k=='s' & df$kb5 == 'HP', 'KP', df$kb5)
df$kb6 <- ifelse(df$type_k=='s' & df$kb6 == 'HP', 'CP', df$kb6)
df$kb7 <- ifelse(df$type_k=='s' & df$kb7 == 'HP', 'KP', df$kb7)
df$kb8 <- ifelse(df$type_k=='s' & df$kb8 == 'HP', 'CP', df$kb8)

write.table(df,file="OM1_combined_ConvertedHPs_SK1.txt", row.names=F, quote=F, sep='\t')

#############################################################################################
###Make HPs File - table containing all positions called as HP in at least 1 spore, for QC###
df <- read.table("OM1_combined_SK1.txt", header=T)
HPs1 <-df[which(df$kb1=='HP'),]
HPs2 <-df[which(df$kb2=='HP'),]
HPs3 <-df[which(df$kb3=='HP'),]
HPs4 <-df[which(df$kb4=='HP'),]
HPs5 <-df[which(df$kb5=='HP'),]
HPs6 <-df[which(df$kb6=='HP'),]
HPs7 <-df[which(df$kb7=='HP'),]
HPs8 <-df[which(df$kb8=='HP'),]

total <- rbind(HPs1, HPs2, HPs3,HPs4, HPs5, HPs6, HPs7, HPs8)
total <- total[!duplicated(total), ]
write.table(total,file="HP-list_OM1_SK1.txt", row.names=F, quote=F, sep='\t')
####End###






