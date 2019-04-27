SNP_indel_isolator_SK1 <- function(df, filename, polys) {
  
  #Chromosome names must be converted to numbers in order to match those in the variant table.
  #(Chromosome names need to be roman numerals for PySamStat)
  df$chrom <- as.character(df$chrom)
  df$chrom[df$chrom == "XVI"] <- "16"
  df$chrom[df$chrom == "XV"] <- "15"
  df$chrom[df$chrom == "XIV"] <- "14"
  df$chrom[df$chrom == "XIII"] <- "13"
  df$chrom[df$chrom == "XII"] <- "12"
  df$chrom[df$chrom == "XI"] <- "11"
  df$chrom[df$chrom == "IX"] <- "9"
  df$chrom[df$chrom == "X"] <- "10"
  df$chrom[df$chrom == "VIII"] <- "8"
  df$chrom[df$chrom == "VII"] <- "7"
  df$chrom[df$chrom == "VI"] <- "6"
  df$chrom[df$chrom == "IV"] <- "4"
  df$chrom[df$chrom == "V"] <- "5"
  df$chrom[df$chrom == "III"] <- "3"
  df$chrom[df$chrom == "II"] <- "2"
  df$chrom[df$chrom == "I"] <- "1"
  
  #making a subtable from the coverage file containing only the important columns
  df<- subset(df, select = c(chrom, pos, reads_all, deletions, insertions, A, C, T, G))
  df <-rename(df, c("pos"="pos_k"))
 
  merged_data <- merge(polys, df, by=c("chrom", "pos_k"), all=FALSE, sort=FALSE)
  #This merges the coverage dataset with the SNP table according to chromosome and position (this is why it's important to have the same chromosome names). 
  #all=FALSE means that any positions in the coverage file that are not in the SNP table are discarded.
  
  merged_data<- subset(merged_data, select = c(uID, sID, chrom, pos_c, pos_k, seq_c, seq_k, type_c,	type_k, Var_len, reads_all, deletions, insertions, A, C, T, G))
  
  write.table(merged_data,file=sprintf("%sOut_SK1.txt",filename), row.names=F, quote=F, sep='\t')
  
  ##Subtables for genotype testing##
  ##note, it is important to do this on the unmerged table, 
  ##as if there are no SNPs in a gene it will give 0 reads
  
  dftel11 <- df[which(df$chrom==2),]
  dftel12 <- df[which(df$pos_k<59382),]
  dftel13 <- df[which(df$pos_k>51019),]
  dftel1 <- merge(dftel11, dftel12, by=c("chrom", "pos_k", "reads_all"), all=FALSE, sort=FALSE)
  dftel1 <- merge(dftel1, dftel13, by=c("chrom", "pos_k", "reads_all"), all=FALSE, sort=FALSE)
  
  dfmsh21 <- df[which(df$pos_k>147382),] 
  dfmsh22 <- df[which(df$pos_k<150276),] 
  dfmsh23 <- df[which(df$chrom==15),] 
  dfmsh2 <- merge(dfmsh21, dfmsh22, by=c("chrom", "pos_k", "reads_all"), all=FALSE, sort=FALSE)
  dfmsh2 <- merge(dfmsh2, dfmsh23, by=c("chrom", "pos_k", "reads_all"), all=FALSE, sort=FALSE)
  
  dfndt801 <- df[which(df$pos_k>356561),] 
  dfndt802 <- df[which(df$pos_k<358444),] 
  dfndt803 <- df[which(df$chrom==8),] 
  dfndt80 <- merge(dfndt801, dfndt802, by=c("chrom", "pos_k", "reads_all"), all=FALSE, sort=FALSE)
  dfndt80 <- merge(dfndt80, dfndt803, by=c("chrom", "pos_k", "reads_all"), all=FALSE, sort=FALSE)
  
  dfrad241 <- df[which(df$pos_k<538279),]
  dfrad242 <- df[which(df$pos_k>536300),]
  dfrad243 <- df[which(df$chrom==5),]
  dfrad24 <- merge(dfrad241, dfrad242, by=c("chrom", "pos_k", "reads_all"), all=FALSE, sort=FALSE)
  dfrad24 <- merge(dfrad24, dfrad243, by=c("chrom", "pos_k", "reads_all"), all=FALSE, sort=FALSE)
  
  dfsml11 <- df[which(df$pos_k>159383),] 
  dfsml12 <- df[which(df$pos_k<159697),] 
  dfsml13 <- df[which(df$chrom==13),] 
  dfsml1 <- merge(dfsml11, dfsml12, by=c("chrom", "pos_k", "reads_all"), all=FALSE, sort=FALSE)
  dfsml1 <- merge(dfsml1, dfsml13, by=c("chrom", "pos_k", "reads_all"), all=FALSE, sort=FALSE)
  
  gene <- c('tel1', 'msh2', 'ndt80', 'rad24', 'sml1')
  reads <- c(0,0,0,0,0)
  gentable <- data.frame(gene, reads)
  gentable[1,2] <- sum(dftel1$reads_all)
  gentable[2,2] <- sum(dfmsh2$reads_all)
  gentable[3,2] <- sum(dfndt80$reads_all)
  gentable[4,2] <- sum(dfrad24$reads_all)
  gentable[5,2] <- sum(dfsml1$reads_all)
  
  write.table(gentable,file=sprintf("%sGenotype_SK1.txt",filename), row.names=F, quote=F, sep='\t')
  
}
