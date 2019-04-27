#This combines all sorted Event tables and gr2 files into large master tables.#

#Version 2: 18 July 2016
#Accepts new style of event tables (events pre-sorted by complexity)
#Version 3: 13 June 2016
#Accepts new event tables which contain more information about individual segments (length, position, SNPs, spores affected)
#"spore affected" indicates which spores were affected by each pattern change
#04-12-17: changed 'debut' from 'start5' to 'start3' for proper imaging of events with few SNPs

#Set wd as required
#the program will make files in a folder named "Output_Files" in your current working directory. Create this folder first if needed.
#Note: if you run the program multiple times, it will write into a pre-existing master table, making duplicated entries. 
#Avoid this by changing the filename of the master table.

# Create new empty master dataframe
master=NULL

library(stringr)

# Sequentially imports each Event Table and combines into a Master Table with additional identified columns for Meiosis, Threshold, Genotype and GenotypeID (e.g. OM2, 1500, OM, 2)
files = list.files(pattern="NEvents_") # import files names with "Event_" string into variable "files" 
files1 = length(files) # Count number of files
files2 = read.table(text = files, sep = "_", as.is = TRUE) #Split file names by "_" separator and create table "files2"

for (j in 1:files1) {
 data=NULL 
 data <- read.table(files[j], sep = "\t", header=TRUE) #Import datatable from files number 1 to j
 data["Meiosis"]=files2[j,2] # Insert Meiosis identifier from files2 into a new column
 data["threshold"]=files2[j,3] # Insert threshold value from files2 into a new column
 data["Genotype"]=gsub('([A-z]+).*','\\1',files2[j,2]) #Extracts character portion of string
 data["GenotypeID"]=gsub('.*([0-9]+)','\\1',files2[j,2])  # Extracts digit portion of string
 data["midpoint"]=(data$start5+data$start3+data$stop5+data$stop3)/4 # Add midpoint column
 data["debut"]=(data$start3) # Add event start column
 data["fin"]=(data$stop3) # Add event end column
 master <- rbind(master, data) # Combine data into master table

}

master2 <- master

# Reorder columns
master2=master2[c("id","Meiosis","threshold","Genotype","GenotypeID","chr","start5",
                  "start3","stop5","stop3","midpoint","debut","fin","type","CO",
                  "chromatids","len_min","len_mid","len_max", "nb_seg", "seg_len_mid", "seg_len_min", "seg_len_max",	
                  "seg_start5", "seg_start3",	"seg_stop5", "seg_stop3", "PanHpM",	"GenomeHpM",	"Obs_Exp_HpM", 
                  "nb_snp", "segnb_snp",	"segspaff", "chr_aff", "groupe", "classe", "Notes","CutMin", "CutMax", "LSB")]

#Write out master files
wd = getwd()
out = paste(wd,"/","/","MasterEventTable",sep="")
write.table(master2, out, col.names =TRUE, row.names = FALSE, quote = FALSE, sep="\t", append=TRUE)


######################

# Create new empty master dataframe
masterGR=NULL
require(stringr)

# Sequentially imports each gr2 Table and combines into a Master Table with additional identified columns for Meiosis, Threshold, Genotype and GenotypeID (e.g. OM2, 1500, OM, 2)
# This is rather slow due to size of gr2 files
filesGR = list.files(pattern="_gr2") # import files names with "_gr2" string into variable "files" 
filesGR1 = length(filesGR) # Count number of files
filesGR2 = read.table(text = filesGR, sep = "_", as.is = TRUE) #Split file names by "_" separator and create table "files2"

for (j in 1:filesGR1)
  
{data <- read.table(filesGR[j], sep = "\t", header=TRUE, stringsAsFactors=FALSE) #Import datatable from files number 1 to j
 data["Meiosis"]=filesGR2[j,1] # Insert Meiosis identifier from filesGR2 into a new column
 data["Genotype"]=gsub('([A-Z]+).*','\\1',filesGR2[j,1]) # Extracts character portion of string
 data["GenotypeID"]=gsub('.*([0-9]+)','\\1',filesGR2[j,1])  # Extracts digit portion of string
 masterGR <- rbind(masterGR, data) # Combine data into masterGR table
}

##############################
# Create deriviative masterGR_noNA table that removes any line containing an NA call in the binary columns.
# WHY? Removes positions with read depth that was too low to make a call.
# This allows haplotype calls to be contiguous across such positions (rather than having gaps when drawn).

masterGR_noNA=na.omit(masterGR)
rownames(masterGR_noNA)=NULL

##############################

# Add inter-SNP interval start and stop locations
masterGR_noNA[,"offset+1"]=masterGR_noNA[c(2:nrow(masterGR_noNA),NA),"pos_c"]
masterGR_noNA[,"offset-1"]=masterGR_noNA[c(NA,1:nrow(masterGR_noNA)-1),"pos_c"]

masterGR_noNA[,"startdif"]=(((masterGR_noNA[,"pos_c"]-masterGR_noNA[,"offset-1"])/2)+abs((masterGR_noNA[,"pos_c"]-masterGR_noNA[,"offset-1"])/2))/2 # Adding together real and absolute values and dividing by 2 creates zeros for negative values
masterGR_noNA[,"stopdif"]=(((masterGR_noNA[,"offset+1"]-masterGR_noNA[,"pos_c"])/2)+abs((masterGR_noNA[,"offset+1"]-masterGR_noNA[,"pos_c"])/2))/2 # Adding together real and absolute values and dividing by 2 creates zeros for negative values
masterGR_noNA[,"startdif"][masterGR_noNA[,"startdif"] == 0] <- NA # Syntax to convert all zeros to NA
masterGR_noNA[,"stopdif"][masterGR_noNA[,"stopdif"] == 0] <- NA # Syntax to convert all zeros to NA

masterGR_noNA[,"startSNP"]=masterGR_noNA[,"startdif"]+masterGR_noNA[,"offset-1"] # Calculate startSNP 
masterGR_noNA[,"stopSNP"]=masterGR_noNA[,"stopdif"]+masterGR_noNA[,"pos_c"] # Calculate stopSNP

##############################

# Add inter-SNP interval start and stop locations to masterGR table
masterGR[,"offset+1"]=masterGR[c(2:nrow(masterGR),NA),"pos_c"]
masterGR[,"offset-1"]=masterGR[c(NA,1:nrow(masterGR)-1),"pos_c"]

masterGR[,"startdif"]=(((masterGR[,"pos_c"]-masterGR[,"offset-1"])/2)+abs((masterGR[,"pos_c"]-masterGR[,"offset-1"])/2))/2 # Adding together real and absolute values and dividing by 2 creates zeros for negative values
masterGR[,"stopdif"]=(((masterGR[,"offset+1"]-masterGR[,"pos_c"])/2)+abs((masterGR[,"offset+1"]-masterGR[,"pos_c"])/2))/2 # Adding together real and absolute values and dividing by 2 creates zeros for negative values
masterGR[,"startdif"][masterGR[,"startdif"] == 0] <- NA # Syntax to convert all zeros to NA
masterGR[,"stopdif"][masterGR[,"stopdif"] == 0] <- NA # Syntax to convert all zeros to NA

masterGR[,"startSNP"]=masterGR[,"startdif"]+masterGR[,"offset-1"] # Calculate startSNP 
masterGR[,"stopSNP"]=masterGR[,"stopdif"]+masterGR[,"pos_c"] # Calculate stopSNP

##############################


#Write out masterGR and masterGR_noNA files

wd = getwd()
out = paste(wd,"/","Output_Files","/","MasterGRTable",sep="")
write.table(masterGR, out, col.names =TRUE, row.names = FALSE, quote = FALSE, sep="\t", append=TRUE)

wd = getwd()
out = paste(wd,"/","Output_Files","/","MasterGRnoNATable",sep="")
write.table(masterGR_noNA, out, col.names =TRUE, row.names = FALSE, quote = FALSE, sep="\t", append=TRUE)

######################

