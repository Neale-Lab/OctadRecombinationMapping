#Use:
#Isolate variant reads (listed in a variant table) from an octad.
#Change filenames with find and replace - replace the meiosis identifier.
#Set working directories as appropriate
#Requirements: Variant table, SNP_indel_isolator function scripts, PySAMStat output for all 8 members of Octad
#Install the plyr package

library(plyr)

#opening the polymorphisms table
#Change the file as appropriate for your samples
polys <- read.table("Variant_Table_6.txt", header=T, stringsAsFactors=FALSE)

#Source the isolator functions
source("SNP_indel_isolator_S288c.R")
source("SNP_indel_isolator_SK1.R")

# opening the S288c files, you will need to change the filenames.
dfc <- read.table("OM1A_S288c.txt",header=T, stringsAsFactors=FALSE)
SNP_indel_isolator_S288c(dfc, "OM1A", polys)
tm1 <-dfc[which(dfc$chrom=='2u'),] #recording the number of reads mapping to 2-micron plasmid for QC
mt1 <-dfc[which(dfc$chrom=='mt'),] #recording the number of reads mapping to mitochondrial DNA for QC
totaltm1 <- sum(tm1$reads_all)
totalmt1 <- sum(mt1$reads_all)
reads1c <- sum(dfc$reads_all)-sum(mt1$reads_all)-sum(tm1$reads_all)

dfc <- read.table("OM1Alpha_S288c.txt",header=T, stringsAsFactors=FALSE)
SNP_indel_isolator_S288c(dfc, "OM1Alpha", polys)
tm2 <-dfc[which(dfc$chrom=='2u'),]
mt2 <-dfc[which(dfc$chrom=='mt'),]
totaltm2 <- sum(tm2$reads_all)
totalmt2 <- sum(mt2$reads_all)
reads2c <- sum(dfc$reads_all)-sum(mt2$reads_all)-sum(tm2$reads_all)

dfc <- read.table("OM1B_S288c.txt",header=T, stringsAsFactors=FALSE)
SNP_indel_isolator_S288c(dfc, "OM1B", polys)
tm3 <-dfc[which(dfc$chrom=='2u'),]
mt3 <-dfc[which(dfc$chrom=='mt'),]
totaltm3 <- sum(tm3$reads_all)
totalmt3 <- sum(mt3$reads_all)
reads3c <- sum(dfc$reads_all)-sum(mt3$reads_all)-sum(tm3$reads_all)

dfc <- read.table("OM1Beta_S288c.txt",header=T, stringsAsFactors=FALSE)
SNP_indel_isolator_S288c(dfc, "OM1Beta", polys)
tm4 <-dfc[which(dfc$chrom=='2u'),]
mt4 <-dfc[which(dfc$chrom=='mt'),]
totaltm4 <- sum(tm4$reads_all)
totalmt4 <- sum(mt4$reads_all)
reads4c <- sum(dfc$reads_all)-sum(mt4$reads_all)-sum(tm4$reads_all)

dfc <- read.table("OM1C_S288c.txt",header=T, stringsAsFactors=FALSE)
SNP_indel_isolator_S288c(dfc, "OM1C", polys)
tm5 <-dfc[which(dfc$chrom=='2u'),]
mt5 <-dfc[which(dfc$chrom=='mt'),]
totaltm5 <- sum(tm5$reads_all)
totalmt5 <- sum(mt5$reads_all)
reads5c <- sum(dfc$reads_all)-sum(mt5$reads_all)-sum(tm5$reads_all)

dfc <- read.table("OM1Gamma_S288c.txt",header=T, stringsAsFactors=FALSE)
SNP_indel_isolator_S288c(dfc, "OM1Gamma", polys)
tm6 <-dfc[which(dfc$chrom=='2u'),]
mt6 <-dfc[which(dfc$chrom=='mt'),]
totaltm6 <- sum(tm6$reads_all)
totalmt6 <- sum(mt6$reads_all)
reads6c <- sum(dfc$reads_all)-sum(mt6$reads_all)-sum(tm6$reads_all)

dfc <- read.table("OM1D_S288c.txt",header=T, stringsAsFactors=FALSE)
SNP_indel_isolator_S288c(dfc, "OM1D", polys)
tm7 <-dfc[which(dfc$chrom=='2u'),]
mt7 <-dfc[which(dfc$chrom=='mt'),]
totaltm7 <- sum(tm7$reads_all)
totalmt7 <- sum(mt7$reads_all)
reads7c <- sum(dfc$reads_all)-sum(mt7$reads_all)-sum(tm7$reads_all)

dfc <- read.table("OM1Delta_S288c.txt",header=T, stringsAsFactors=FALSE)
SNP_indel_isolator_S288c(dfc, "OM1Delta", polys)
tm8 <-dfc[which(dfc$chrom=='2u'),]
mt8 <-dfc[which(dfc$chrom=='mt'),]
totaltm8 <- sum(tm8$reads_all)
totalmt8 <- sum(mt8$reads_all)
reads8c <- sum(dfc$reads_all)-sum(mt8$reads_all)-sum(tm8$reads_all)

# Opening SK1 files 

dfk <- read.table("OM1A_SK1.txt",header=T, stringsAsFactors=FALSE)
SNP_indel_isolator_SK1(dfk, "OM1A", polys)
reads1k <- sum(dfk$reads_all)

dfk <- read.table("OM1Alpha_SK1.txt",header=T, stringsAsFactors=FALSE)
SNP_indel_isolator_SK1(dfk, "OM1Alpha", polys)
reads2k <- sum(dfk$reads_all)

dfk <- read.table("OM1B_SK1.txt",header=T, stringsAsFactors=FALSE)
SNP_indel_isolator_SK1(dfk, "OM1B", polys)
reads3k <- sum(dfk$reads_all)

dfk <- read.table("OM1Beta_SK1.txt",header=T, stringsAsFactors=FALSE)
SNP_indel_isolator_SK1(dfk, "OM1Beta", polys)
reads4k <- sum(dfk$reads_all)

dfk <- read.table("OM1C_SK1.txt",header=T, stringsAsFactors=FALSE)
SNP_indel_isolator_SK1(dfk, "OM1C", polys)
reads5k <- sum(dfk$reads_all)

dfk <- read.table("OM1Gamma_SK1.txt",header=T, stringsAsFactors=FALSE)
SNP_indel_isolator_SK1(dfk, "OM1Gamma", polys)
reads6k <- sum(dfk$reads_all)

dfk <- read.table("OM1D_SK1.txt",header=T, stringsAsFactors=FALSE)
SNP_indel_isolator_SK1(dfk, "OM1D", polys)
reads7k <- sum(dfk$reads_all)

dfk <- read.table("OM1Delta_SK1.txt",header=T, stringsAsFactors=FALSE)
SNP_indel_isolator_SK1(dfk, "OM1Delta", polys)
reads8k <- sum(dfk$reads_all)

##Compile a table of total reads for quality control##

Spore <- c("A", "alpha", "B", "beta", "C", "gamma", "D", "delta")
reads_c <- c(reads1c, reads2c, reads3c, reads4c, reads5c, reads6c, reads7c, reads8c)
reads_k <- c(reads1k, reads2k, reads3k, reads4k, reads5k, reads6k, reads7k, reads8k)
reads_mito <- c(totalmt1, totalmt2, totalmt3, totalmt4, totalmt5, totalmt6, totalmt7, totalmt8)
reads_plasmid <- c(totaltm1, totaltm2, totaltm3, totaltm4, totaltm5, totaltm6, totaltm7, totaltm8)

table <- data.frame(Spore, reads_c, reads_k, reads_mito, reads_plasmid)

write.table(table,file="OM1TotalreadsKC.txt",row.names=F, quote=F, sep='\t')

