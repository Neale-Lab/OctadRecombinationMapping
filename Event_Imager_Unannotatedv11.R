#Program to make images of each recombination event

#Change notes:
# Builds on FonctionsAnalyseSpo11RMM written by Marie-Claude
# Version 2: 14th August 2015: Requires updated EventsM1p2MATTv1 for Event Table generation.
  # Specifically, the values sent to graph_bin_annot3 command are now the "mid" coordinates of where the events start and finish.
  # FonctionsAnalyse v2 now extends by an arbitrary amount both sides of each event (e.g. 5kb) and highlights the event being queried in the centre.
# Version 3: 15th August
  #Now polls the linked Event Table in order to draw all events that may be visible in the field of view.
  # Note: I currently can't handle chromosomes with NA values in the CO-NCO column. These should be altered to a different value: "U"? for unclassified
  # Update: Now expects and accepts unclassifiable events types "U" in the CO-NCO column
# Version 4: 18th August
  # Now expects and requires a combined Events table (master) and combined gr2 table (masterGR). Use Combine Events Tables.R script for this job FIRST
  # Function call:graph_bin_annot4<-function(meiosis,chromo,debut,fin,ordre=1:8) requires the "Meiosis" identifier to be included
# Version 5: 18th August
  # Rather than accepting a function call, it now steps through the specified draw table, sequentially drawing each event
# Version 6: 19th August
  # Now includes plots for the called event lengths at three specified merge thresholds (500bp, 1500bp, 5000bp)
  # Note: These merge threshold event files must be generated and combined into the master event table
  # Now includes event midpoints in bars below the binary calls
# Version 7: 26th August
  # Now  allows the main overlaying of the chosen event threshold (prior versions alwyas overlaid 1500 threshold) independently of event threshold specification.
  # i.e. You can now choose to look at all events called with 5000 threshold, but specify the main overlay as 1500 or 500 for example.
# Version 8: 10th September
  # Now draws inter-SNP bars, creating joined up plots. Requires updated masterGR and masterGR_noNA tables that contain calculations for inter-SNP midpoints with and without removing NA rows respectively.
  # Depending on which line is hashed out on lines 137/138 (i.e. which GR table is referenced) NA-containing SNPs are either not drawn, are are "joined across".
  # Script now also avoids drawing errors when rarely a region of interest contains no SNPs. Solved by while loop on line 88
  # Now also plots Red1 profile
#Version 9: Marg, 27 May 2016 
  #Plots variant positions above main graph
  #plots indel locations and lengths. Triangle point down - deletion, triangle point up - insertion (relative to S288c)
  #Requires new-style gr2 files which contain variant length and event type columns.
#Version 10: Marg, 18 July 2016
  #Accepts new style of event tables (events pre-sorted by complexity)
  #Has more information in the figure legends
  #Plots thresholds of 1kb, 1.5kb and 2kb.
  #Plots genes, hotspot strength values and Rec8 peaks 
#Version 11: Marg, 12 July 2017
  #has smaller bars representing chromatids, which makes the events easier to read
  #uses Neeman's Spo11-Oligo data (smoothed on the fly)

######################
#Note: If you have just run the Combining tables script, master2, masterGR and masterGR_noNA will be in the environment already 
#so there is no need to reopen the files.
master2<-read.delim("MasterEventTable")
master2$Meiosis <- as.character(master2$Meiosis)
master2$debut <- as.numeric(master2$debut)
master2$fin <- as.numeric(master2$fin)
masterGR_noNA <- read.delim("MasterGRnoNATable")
masterGR <- read.delim("MasterGRTable")
masterGR$Var_len2 <- (masterGR$Var_len-1)
masterGR_noNA$Var_len2 <- (masterGR_noNA$Var_len-1)
#the variant lengths stored in this column are actually 1bp too long because they include the first base upstream. 
#Subtract 1 to make it more accurate.

master2 <- master2[complete.cases(master2$CO),]

#Open the datatracks
rmm_file<-read.delim("RMMSubSamp_simple")
AllElementsDUB = read.table("AllElementsDUB.txt", sep = "\t", header=TRUE) #Import datatable
rec_file<-read.delim("Rec8Marg.txt")
spo_file <- read.delim("WT_SCFullMap_Neeman.txt")
TotalHits=sum(spo_file$Watson+spo_file$Crick)
spo_file$TotalHits=(spo_file$Watson+spo_file$Crick)
spo_file$TotalHpM=(spo_file$Watson+spo_file$Crick)/TotalHits*1000000 # Convert Pan data to Hits per million reads (HpM)
hotspots<- read.delim("WT_Hotspot.Table_Neeman.txt")

# Subset the datatable
draw<-master2 # Temp holding table for events that will be drawn. 

# amend the following lines to subset the datatable as required
draw <- draw[which(draw$threshold==1500),] # Calling threshold of events to draw.
mainoverlaythreshold=1500 # The variable that specifies the threshold to use for main event painting/overlay
draw=subset(draw,Meiosis =="OM1" ) # (refers to a specific meiosis)
#you can also subset this on other variables e.g. if you want to draw all crossovers from one genotype, etc.

#######################

ordre=1:8 # Order of plotting spores
extend=5000 # How many bp to extend the plotting either side of the event

# Automated PDF generation
wd = getwd(); out = paste(wd,"/","Output_Files","/","OM1_Images_±", extend,"bp_",".pdf",sep=""); pdf(file=out, width=21,height=9);

# Replace funtion call with a loop that steps through the draw subtable

if(nrow(draw)==0){print("no lines in draw")}

if(nrow(draw)>0){
  for (j in 1: nrow(draw) ) {
    meiosis=draw[j,"Meiosis"]
    chromo=draw[j,"chr"]
    debut=draw[j,"debut"]
    fin=draw[j,"fin"]
    
    # How much to extend the plotting either side of the event
    debut=debut-extend
    fin=fin+extend
  	
  # This sets up the graphing parameters: mfrow calculates number of rows of grpahics to be drawn. mar sets the margins of the graph
  layout(matrix(c(0,1,1,1,2:17,18,18,19,19,20,20,21,22,22,23,23,24,24,24,0),35, 1, byrow = T))
  par(mar=c(0,6,0,2),oma = c(4,0,0,0),las=1) # Sets margins per graph and outside margins per grouped set
  
  # Subset the masterGR tables for the region of interest
  a<-subset(masterGR,Meiosis==meiosis & chr == chromo & (debut) <= pos_c & pos_c <= (fin)) 
  a1<-subset(masterGR_noNA,Meiosis==meiosis & chr == chromo & (debut) <= pos_c & pos_c <= (fin)) 
  
  while (nrow(a1)<20){ # Check whether the region of interest contains sufficient SNPs to draw - this is important to prevent an error when a/a1 contain no rows due to rare zero SNP events (i.e. COs) that also lack SNPs in flanking thresholded region.
    debut=debut-1000
    fin=fin+1000
    a<-subset(masterGR,Meiosis==meiosis & chr == chromo & (debut) <= pos_c & pos_c <= (fin)) 
    a1<-subset(masterGR_noNA,Meiosis==meiosis & chr == chromo & (debut) <= pos_c & pos_c <= (fin)) 
  }
  
if(nrow(a)>0){

  # Subset master based on Meiosis and Chromosome being drawn
  e2=subset(master2, Meiosis==meiosis & chr == chromo)
  
  # Specify event colours based on type (CO or NCO) - fill in master table with this information
  for (i in 1:nrow(e2)){
    if(e2[i,"CO"]==1){e2[i,"eventcolour"]="wheat"}
    if(e2[i,"CO"]==0){e2[i,"eventcolour"]="thistle2"}
    if(e2[i,"CO"]==2){e2[i,"eventcolour"]="wheat"} 
    if(e2[i,"CO"]==3){e2[i,"eventcolour"]="wheat"}
    if(e2[i,"CO"]=="U"){e2[i,"eventcolour"]="powderblue"}
  }
  
  #Subset e2 based on threshold values so that each can be plotted
  e1000=subset(e2, Meiosis==meiosis & chr == chromo & threshold==1000)
  e1500=subset(e2, Meiosis==meiosis & chr == chromo & threshold==1500)
  e2000=subset(e2, Meiosis==meiosis & chr == chromo & threshold==2000)
  mainoverlay=subset(e2,Meiosis==meiosis & chr == chromo & threshold==mainoverlaythreshold)
  
  ######add SNP coordinates as text######
  plot(a$pos_c,a[,1],type="n",ylim=c(0,300),ylab=paste("SNPs"),cex.lab=1,font=1,xlab="",xlim=c(debut,fin), xaxt="n",yaxt="n", axes=FALSE, xaxs="i",yaxs="i")
  #type=n for no plotting. axes=FALSE gets rid of border. 
  text(a$pos_c, y=100, labels = a$pos_c, srt=90, pos=3, cex=1) 
  #cex=1 controls the size of the text; however, making it smaller doesn't really help with overlap
  arrows(a$pos_c,0,a$pos_c,50,col="black",length=0,lwd=1, code=3, lend=1)
  
  #######################################
  
  # Graphs of raw data
  for (i in (ordre+2)){  
  plot(a$pos_c,a[,i],type="l",lwd=2,col=2,cex.axis=0.8,ylim=c(0.6,300),yaxp=c(1,100,1),ylab=paste("S",i-2,sep=""),cex.lab=1,font=1,xaxs="i",yaxs="i",xlab="",xlim=c(debut,fin), xaxt="n", log="y")
  
  # We impose xlim=c(start,end) to wedge with spo data and rmm
  
  points(a$pos_c,a[,i],pch=4,col=2,cex=2)
  points(a$pos_c,a[,(i+8)],pch=4,col=4,cex=2)
  lines(a$pos_c,a[,(i+8)],col=4,lwd=1)
  points(a$pos_c,rep(100,length(a$pos_c)),pch=-1*a[,(i+16)]+17,cex=1.5,col=-2*a[,(i+16)]+4)
  } #End of i first loop
  
  # Graphs of binary data
  position=c(85,95,85,95,85,95,85,95)
  haut<-0
  for (i in (ordre+2)){
  plot(a$pos_c,rep(100,length(a$pos_c)),xaxt="n",yaxt="n",type="n",ylab=paste("S",i-2,sep=""),cex.lab=1,font=2,xaxs="i",yaxs="i",ylim=c(50,130),axes=FALSE,xlim=c(debut,fin))
  
  # On impose xlim=c(debut,fin) pour caler avec les données de spo et rmm
  alt<-ifelse(haut==0,120,120)
  bas<-ifelse(haut==0,60,60)
  arrows(mainoverlay$debut,90,mainoverlay$fin,90,col=mainoverlay$eventcolour,length=0,lwd=30, lend=1) # Underlay a highlighting box across all event regions
  arrows(a$pos_c,position[i-2]-30,a$pos_c,position[i-2]+30,col=-2*a[,(i+16)]+4,length=0,lwd=1)
  
  #To alter the thickness of the chromatid bars, change lwd on this line
  #arrows(a$startSNP,position[i-2],a$stopSNP,position[i-2],col=-2*a[,(i+16)]+4,length=0,lwd=15, lend=1) # Add SNP start/stop boundary overlays
  arrows(a1$startSNP,position[i-2],a1$stopSNP,position[i-2],col=-2*a1[,(i+16)]+4,length=0,lwd=10, lend=1) # Add SNP start/stop boundary overlays
  
  haut<-1-haut
  
  } #End of i second loop
  
  # Top <-1 Top toggles the lengths of the lines between top = 0 and top = 1
  
  # Addition of plots that indicate the events specified by 1000 and 2000bp event-calling thresholds
  plot(a$pos_c,rep(100,length(a$pos_c)),xaxt="n",yaxt="n",type="n",ylab=paste("Events"),cex.lab=1.25,font=2,xaxs="i",yaxs="i",ylim=c(0,165),axes=FALSE,xlim=c(debut,fin))
  arrows(e1000$debut,140,e1000$fin,140,col=e1000$eventcolour,length=0,lwd=12, code=3, lend=1); text(debut,140, "1.0 kb merge",pos=4) # Underlay a highlighting box across all event regions
  arrows(e1500$debut,100,e1500$fin,100,col=e1500$eventcolour,length=0,lwd=12, code=3, lend=1); text(debut,100, "1.5 kb merge",pos=4) # Underlay a highlighting box across all event regions
  arrows(e2000$debut,60,e2000$fin,60,col=e2000$eventcolour,length=0,lwd=12, code=3, lend=1); text(debut,60, "2.0 kb merge",pos=4) # Underlay a highlighting box across all event regions
  
  # Plot event midpoints
  arrows(e1000$midpoint,130,e1000$midpoint,150,col="black",length=0,lwd=1, code=3, lend=1)
  arrows(e1500$midpoint,90,e1500$midpoint,110,col="black",length=0,lwd=1, code=3, lend=1)
  arrows(e2000$midpoint,50,e2000$midpoint,70,col="black",length=0,lwd=1, code=3, lend=1)
  
  ####################################################################
  ###plot indels as triangles with length underneath###
  plot(a$pos_c,a[,1],type="n",ylim=c(0,300),ylab=paste("Indels"),cex.lab=1.25,font=2,xlab="",xlim=c(debut,fin), xaxt="n",yaxt="n", axes=FALSE, xaxs="i",yaxs="i")
  #type=n for no plotting. axes=FALSE gets rid of border. 
  aDEL=subset(a, type_k=="d")
  if(nrow(aDEL)>0){
  yDEL =rep(200, nrow(aDEL))
  points(aDEL$pos_c, yDEL, type = "p", pch=25, col="dimgrey", cex=2) 
  text(aDEL$pos_c, (yDEL-145), labels = aDEL$Var_len2, pos=3, cex=1) 
  }
  
  aINS=subset(a, type_k=="i")
  if(nrow(aINS)>0){
  yINS =rep(200, nrow(aINS))
  points(aINS$pos_c, yINS, type = "p", pch=24, col="dimgrey", cex=2) 
  text(aINS$pos_c, (yINS-145), labels = aINS$Var_len2, pos=3, cex=1) 
  }
  
  ####################################################################
  
  # Addition des données de Spo11 and RMM
  rbPal <- colorRampPalette(c('black','blue','red','pink')) # Creates a colour palette if chosing to plot Spo11 as a coloured histogram "h"
  spo_file$Col <- rbPal(10)[as.numeric(cut(spo_file$TotalHits,breaks = 10))]
  
  map<-subset(spo_file,Chr == chromo & debut <= Pos & Pos <= fin)
  map2<-subset(rmm_file,chr == chromo & debut <= pos & pos <= fin)
  map5<-subset(hotspots,Chr == chromo & debut <= Start & End <= fin) #Note: Pan file variables are called CHROM, HS_START, HS_END
  map6<-subset(rec_file,Chr == chromo & debut <= Position & Position <= fin)
  
  # Hanninng window code:
  require("e1071") # This package permits calculation of hanning function
  
  #Decompression code here. Because data is sparse, this is needed in order to populate all the missing rows prior to smoothing
  spo11.0=subset(spo_file,Chr == chromo & debut <= Pos & Pos <= fin) #Make a sub-table of the DSB data that only contains those rows where chr = 1 in range of interest
  spo11.1 <- data.frame(Chr=chromo, Pos=(debut:fin)) # Creates expanded dataframe with ONLY Chr and Pos locations
  spo11.1 <- merge(spo11.1,spo11.0, all=TRUE) # Merge expanded empty dataframe with compressed spo11.0 dataframe
  spo11.1[is.na(spo11.1)] <- 0 # Convert all NA values to zero
  
  win=101 #Smoothing Window length
  hw=hanning.window(win) #create hanning window (require package e1071 to be loaded)
  temp=c(rep(0,win),spo11.1$TotalHpM, rep(0,win)) # Extend by the length of the sliding window with zeros at both ends
  spo11.s=filter(temp,hw) # smooth the temp vector using the hann window and the filter function
  spo11.s=spo11.s[(win+1):(length(spo11.s)-win)] # trim smooth back to correct length
  
  spo11.1$Spo11smoothed=spo11.s #Write into starting spo11.0 tabel as new column
  
  plot(spo11.1$Pos,spo11.1$Spo11smoothed,type="l",xlim=c(debut,fin),ylim=c(10,max(spo11.1$Spo11smoothed)),yaxp=c(100,round(max(spo11.1$Spo11smoothed),-2),2),ylab="Spo11",bty="u",xaxs="i",lwd=2,xaxt="n",cex.axis=0.9,cex.lab=1.25)
  
  ###plot hotspots###
  plot(map5$MedianPoint,map5$Total_HpM,type="n",ylim=c(0,300),ylab="",cex.lab=1.25,font=2,xlab="",xlim=c(debut,fin), xaxt="n",yaxt="n", axes=FALSE, xaxs="i",yaxs="i")
  #type=n for no plotting. axes=FALSE gets rid of border.
  if(nrow(map5)>0){
    arrows(map5$Start,140,map5$End,140,col='lightsteelblue',length=0,lwd=12, code=3, lend=1); text="" # Underlay a highlighting box across all event regions
    text(map5$MedianPoint, y=0, labels = map5$TotalHpM, pos=3, cex=1) 
  }
  
  ##########################################################################################################################################################
  #Now plot the gene datatrack
  #First subset the relevant data
  genes=AllElementsDUB #First make a copy of the ALLElements table
  genes=subset(genes,chr==chromo & start>(debut-10000) & stop<(fin+10000)) #Make a sub-table of ALLElements where chr = 1 and has limits just beyond plot range
  genes=subset(genes,type=="gene") #Make a sub-table of ALLElements
  #Now perform the plot
  plot(genes$start,genes$start, xaxt="n",yaxt="n",type="n", ylab=paste("Genes"),cex.lab=1.5,font=2, xlim=c(debut,fin), ylim=c(0,120),axes=FALSE) #set up empty plot
  # Following module draws arrows for each element
  xrange=fin-debut
  ahead=xrange/25 #make arrowhead length proportional to plot range
  ahead[(ahead>500)]=500 #limit max length to 500
  av=75 #arrow vertical location relative to plot dimensions
  ahw=15 #arrow/head width
  
  genesW=subset(genes,genename !="Dubious_ORF" & orientation =="+") #Make a sub-table of ALLElements
  if(nrow(genesW)>0){
  for (i in 1:nrow(genesW)){
    polygon(c(genesW[i,"start"], genesW[i,"stop"]-ahead, genesW[i,"stop"]-ahead,genesW[i,"stop"], genesW[i,"stop"]-ahead, genesW[i,"stop"]-ahead, genesW[i,"start"]),c(av+ahw,av+ahw,av+ahw+ahw,av,av-ahw-ahw,av-ahw, av-ahw), col="palegreen", border="palegreen4")
    text((genesW[i,"start"]+genesW[i,"stop"])/2,av, font=3, genesW[i,"genename"], cex=0.9) }
  }
    
  genesW=subset(genes,genename=="Dubious_ORF" & orientation =="+") #Make a sub-table of ALLElements
  if(nrow(genesW)>0){
  for (i in 1:nrow(genesW)){
    polygon(c(genesW[i,"start"], genesW[i,"stop"]-ahead, genesW[i,"stop"]-ahead,genesW[i,"stop"], genesW[i,"stop"]-ahead, genesW[i,"stop"]-ahead, genesW[i,"start"]),c(av+ahw,av+ahw,av+ahw+ahw,av,av-ahw-ahw,av-ahw, av-ahw), col="palegreen", border="palegreen4", lty=2)
    text((genesW[i,"start"]+genesW[i,"stop"])/2,av, font=3, genesW[i,"sysname"], cex=0.9) }
  }
  
  av=25 #arrow vertical location for Crick genes relative to plot dimensions
  genesC=subset(genes,genename !="Dubious_ORF" & orientation =="-") #Make a sub-table of ALLElements
  if(nrow(genesC)>0){
  for (i in 1:nrow(genesC)){
    polygon(c(genesC[i,"stop"], genesC[i,"start"]+ahead, genesC[i,"start"]+ahead,genesC[i,"start"], genesC[i,"start"]+ahead, genesC[i,"start"]+ahead, genesC[i,"stop"]),c(av+ahw,av+ahw,av+ahw+ahw,av,av-ahw-ahw,av-ahw, av-ahw), col="lightpink", border ="lightpink4")
    text((genesC[i,"start"]+genesC[i,"stop"])/2,av, font=3, genesC[i,"genename"], cex=0.9) }
  }
  
  genesC=subset(genes,genename=="Dubious_ORF" & orientation =="-") #Make a sub-table of ALLElements
  if(nrow(genesC)>0){
  for (i in 1:nrow(genesC)){
    polygon(c(genesC[i,"stop"], genesC[i,"start"]+ahead, genesC[i,"start"]+ahead,genesC[i,"start"], genesC[i,"start"]+ahead, genesC[i,"start"]+ahead, genesC[i,"stop"]),c(av+ahw,av+ahw,av+ahw+ahw,av,av-ahw-ahw,av-ahw, av-ahw), col="lightpink", border ="lightpink4", lty=2)
    text((genesC[i,"start"]+genesC[i,"stop"])/2,av, font=3, genesC[i,"sysname"], cex=0.9) }
  }
  ####Plot RMM track###
  plot(map2$pos,map2$rmm,type="l",xlim=c(debut,fin),ylim=c(1,max(map2$rmm)),yaxp=c(1,round(max(map2$rmm),0),2),ylab="RMM",bty="u",xaxs="i",lwd=4,col='darkgoldenrod',
       xaxp=c((round(debut,-3)),(round(fin,-3)),(round(round(fin,-3)-round(debut,-3)) / 500)),cex.axis=0.9, cex.lab=1.25)
  if(nrow(map6)>0){
  arrows(map6$Position,0,map6$Position,50,col="black",length=0,lwd=1, code=3, lend=1)
  }
  ###Make figure legend###
  
  title(xlab = paste(sep="","Red=S288c, Blue=SK1",
                     "\n",
                     meiosis, "  Id=", draw[j,"id"],
                     "  Chromosome=",chromo,
                     "  Threshold=", draw[j,"threshold"], " bp",
                     "  Mid Length=", draw[j,"len_mid"]," bp",
                     "\n",
                     "  Group=", draw[j,"groupe"],
                     "  Class=", draw[j,"classe"],
                     "  No. of COs=", draw[j,"CO"],
                     "  Variants=", draw[j,"nb_snp"],
                     "  Spo11 HpM=", draw[j,"PanHpM"],
                     "  Type=", draw[j,"type"],
                     "  Sp.Aff=", draw[j,"segspaff"]), 
        outer = T, line = 2, cex.lab=1.4) # Cex.lab controls the labelling size
}
  } #End of j loop
}
dev.off()

