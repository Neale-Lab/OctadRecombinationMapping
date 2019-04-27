#This makes a small modification to the Annot file, moving the values in the sporeaffect column down

#the function
annotmodifier <- function(df, filename){
  
dfN<- NULL
for(i in 1:16){
  subset1 <- df[which(df$chr==i),]
  subset2 <- df[which(df$chr==i),]
  for(j in 1:nrow(subset1)){
    if(j==1){subset1[j,9] <- 0}
    if(j>1){subset1[j,9] <- subset2[(j-1),9]}
  }
  dfN <- rbind(dfN, subset1)
  
}
write.table(dfN,file=sprintf("NAnnot_%s",filename),row.names=F, quote=F, sep='\t')
}

#running the function (replace filename for each sample you have)
df1<-read.delim("Annot_OM1")
annotmodifier(df1, 'OM1')

df2<-read.delim("Annot_OM2")
annotmodifier(df2, 'OM2')

df1<-read.delim("Annot_OM3")
annotmodifier(df1, 'OM3')

df1<-read.delim("Annot_OM4")
annotmodifier(df1, 'OM4')

df1<-read.delim("Annot_OM5")
annotmodifier(df1, 'OM5')

df1<-read.delim("Annot_OM6")
annotmodifier(df1, 'OM6')

df1<-read.delim("Annot_OM7")
annotmodifier(df1, 'OM7')

df1<-read.delim("Annot_OM8")
annotmodifier(df1, 'OM8')

df1<-read.delim("Annot_OM9")
annotmodifier(df1, 'OM9')

df1<-read.delim("Annot_OM10")
annotmodifier(df1, 'OM10')

df1<-read.delim("Annot_OM11")
annotmodifier(df1, 'OM11')

df1<-read.delim("Annot_OM12")
annotmodifier(df1, 'OM12')

df1<-read.delim("Annot_OM13")
annotmodifier(df1, 'OM13')



df1<-read.delim("Annot_OMN4")
annotmodifier(df1, 'OMN4')

df1<-read.delim("Annot_OMN5")
annotmodifier(df1, 'OMN5')

df1<-read.delim("Annot_OMN6")
annotmodifier(df1, 'OMN6')

df1<-read.delim("Annot_OM1")
annotmodifier(df1, 'OM1')
