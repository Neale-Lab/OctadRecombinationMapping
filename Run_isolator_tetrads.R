#Use:
#Isolate variant reads (listed in a variant table) from a tetrad.
#Change filenames with find and replace - replace the meiosis identifier.
#Set working directories as appropriate
#Requirements: Variant table, SNP_indel_isolator function scripts, PySAMStat output for all 4 members of tetrad
#Install the plyr package

library(plyr)

#opening the polymorphisms table
#Change the file as appropriate for your samples
polys <- read.table("Variant_Table_6.txt", header=T, stringsAsFactors=FALSE)

#Source the isolator functions
source("SNP_indel_isolator_S288c.R")
source("SNP_indel_isolator_SK1.R")

# opening the S288c files, you will need to change the filenames.
dfc <- read.table("WT1A_S288c.txt",header=T, stringsAsFactors=FALSE)
SNP_indel_isolator_S288c(dfc, "WT1A", polys)
tm1 <-dfc[which(dfc$chrom=='2u'),] #recording the number of reads mapping to 2-micron plasmid for QC
mt1 <-dfc[which(dfc$chrom=='mt'),] #recording the number of reads mapping to mitochondrial DNA for QC
totaltm1 <- sum(tm1$reads_all)
totalmt1 <- sum(mt1$reads_all)
reads1c <- sum(dfc$reads_all)-sum(mt1$reads_all)-sum(tm1$reads_all)

dfc <- read.table("WT1B_S288c.txt",header=T, stringsAsFactors=FALSE)
SNP_indel_isolator_S288c(dfc, "WT1B", polys)
tm2 <-dfc[which(dfc$chrom=='2u'),]
mt2 <-dfc[which(dfc$chrom=='mt'),]
totaltm2 <- sum(tm2$reads_all)
totalmt2 <- sum(mt2$reads_all)
reads2c <- sum(dfc$reads_all)-sum(mt2$reads_all)-sum(tm2$reads_all)

dfc <- read.table("WT1C_S288c.txt",header=T, stringsAsFactors=FALSE)
SNP_indel_isolator_S288c(dfc, "WT1C", polys)
tm3 <-dfc[which(dfc$chrom=='2u'),]
mt3 <-dfc[which(dfc$chrom=='mt'),]
totaltm3 <- sum(tm3$reads_all)
totalmt3 <- sum(mt3$reads_all)
reads3c <- sum(dfc$reads_all)-sum(mt3$reads_all)-sum(tm3$reads_all)

dfc <- read.table("WT1D_S288c.txt",header=T, stringsAsFactors=FALSE)
SNP_indel_isolator_S288c(dfc, "WT1D", polys)
tm4 <-dfc[which(dfc$chrom=='2u'),]
mt4 <-dfc[which(dfc$chrom=='mt'),]
totaltm4 <- sum(tm4$reads_all)
totalmt4 <- sum(mt4$reads_all)
reads4c <- sum(dfc$reads_all)-sum(mt4$reads_all)-sum(tm4$reads_all)

# Opening SK1 files 
dfk <- read.table("WT1A_SK1.txt",header=T, stringsAsFactors=FALSE)
SNP_indel_isolator_SK1(dfk, "WT1A", polys)
reads1k <- sum(dfk$reads_all)

dfk <- read.table("WT1B_SK1.txt",header=T, stringsAsFactors=FALSE)
SNP_indel_isolator_SK1(dfk, "WT1B", polys)
reads2k <- sum(dfk$reads_all)

dfk <- read.table("WT1C_SK1.txt",header=T, stringsAsFactors=FALSE)
SNP_indel_isolator_SK1(dfk, "WT1C", polys)
reads3k <- sum(dfk$reads_all)

dfk <- read.table("WT1D_SK1.txt",header=T, stringsAsFactors=FALSE)
SNP_indel_isolator_SK1(dfk, "WT1D", polys)
reads4k <- sum(dfk$reads_all)

######################################################
##Compile a table of total reads for quality control##

Spore <- c("A", "B", "C", "D")
reads_mito <- c(totalmt1, totalmt3, totalmt5, totalmt7)
reads_plasmid <- c(totalWT1, totaltm3, totaltm5, totaltm7)
reads_c <- c(reads1c, reads2c, reads3c, reads4c)
reads_k <- c(reads1k, reads2k, reads3k, reads4k)

table <- data.frame(Spore, reads_c, reads_k, reads_mito, reads_plasmid)
write.table(table,file="WT1totalreadsKC.txt",row.names=F, quote=F, sep='\t')


