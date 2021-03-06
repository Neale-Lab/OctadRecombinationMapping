SNP_indel_caller_SK1 <- function(df, filename, threshold_s, threshold_i, threshold_k){
  
  #remove mitochondrial reads
  df <- df[which(df$chrom!='mt'),]
  
  #creates two new columns c and k containing the number of SK1 and SK1 reads respectively for each position
  #Because the reads were aligned to SK1, any indels will be relative to SK1. 
  #Thus for Sk1 reads I take insertions or deletions, but do not need to do the same for SK1.
  
  df$c<-ifelse(df$type_c=='d', df$deletions, 
        ifelse(df$type_c=='i', df$insertions,
        ifelse(df$type_c=='s' & df$seq_c == "T", df$T, 
        ifelse(df$type_c=='s' & df$seq_c == "A", df$A, 
        ifelse(df$type_c=='s' & df$seq_c == "C", df$C, 
        ifelse(df$type_c=='s' & df$seq_c == "G", df$G, 
        0))))))
  
  df$ktemp<-ifelse(df$seq_k == "T", df$T, 
        ifelse(df$seq_k == "A", df$A, 
        ifelse(df$seq_k == "C", df$C, 
        ifelse(df$seq_k == "G", df$G, 
        0))))

  
  df$k<-ifelse(df$type_c=='i', df$ktemp-df$c, 
        df$ktemp)
  df$k <- ifelse(df$k < 0, 0, df$k)
  
  #Note: the first base of an in/del will always be called as reference now since it is the same in both backgrounds
  #For insertions relative to the reference, insertion hits are subtracted from reference hits 
  #because the first base of an insertion is common to both parents, meaning it will always get 'reference' hits.
  
  df$SNPreads <- (df$c+df$k) #total number of SNP (C and K) reads
  
  df$reads_all <-ifelse(df$type_c=='s', (df$reads_all+df$insertions), df$reads_all)
  
  #insertions don't count seperately to 'reads_all'.
  #When the type is S, we are expecting a SNP.
  #if there are insertion or deletion reads, they should count against it being called as a SNP.
  #deletion reads already count towards reads_all.
  
  ##Calling the genotypes##
  
  #CP for SK1 (higher confidence = > 5 reads)
  #KP for SK1 (higher confidence = > 5 reads)
  #CF for SK1 with fewer than 5 reads
  #KF for SK1 with fewer than 5 reads
  #HP for het with = > 5 reads
  #HF for het with < 5 reads
  #NP for ambiguous with = >5 reads
  #NF for ambiguous with <5 reads

df$kb <-  ifelse(df$type_c=='i' & df$reads_all<threshold_i & (df$c/(df$reads_all))>=0.3, "CF",
          ifelse(df$type_c=='i' & df$reads_all>=threshold_i & (df$c/(df$reads_all))>=0.3, "CP",
          ifelse(df$type_c=='i' & df$reads_all>=threshold_i & df$c>=threshold_k & (df$c/(df$reads_all))>=0.2, "CP",
          ifelse(df$type_c=='i' & df$reads_all>=threshold_i & df$Var_len>=20 & (df$c/(df$reads_all))>=0.2, "CP",
#lower pass threshold for long insertions and positions with at least 10 insertion hits
          ifelse(df$type_c=='i' & df$reads_all>=threshold_i & df$c<2 & (df$k/(df$reads_all))>=0.95, "KP",
          ifelse(df$type_c=='i' & df$reads_all<threshold_i & df$c<2 & (df$k/(df$reads_all))>=0.95, "KF",
#if there aren't any insertion reads, it's probably SK1 (unless read depth is low).
#But need to check that the majority of reads match the reference as well, else it might be a variant table error
          ifelse(df$type_c=='d' & df$reads_all<threshold_i & (df$c/(df$reads_all))>=0.3, "CF",
          ifelse(df$type_c=='d' & df$reads_all>=threshold_i & (df$c/(df$reads_all))>=0.3, "CP",
          ifelse(df$type_c=='d' & df$reads_all>=threshold_i & df$c>=threshold_k & (df$c/(df$reads_all))>=0.2, "CP",
          ifelse(df$type_c=='d' & df$reads_all>=threshold_i & df$Var_len>=20 & (df$c/(df$reads_all))>=0.2, "CP",
#as with insertions there is a reduced threshold for positions with at least 10 deletion hits, and for longer deletions.
#deletion thresholds are higher than insertion thresholds because they aren't affected by the first base being the same,
#but they still have a lower theshold than SNPs because they don't align as well   
          ifelse(df$type_c=='d' & df$reads_all<threshold_i & df$c<2 &  (df$k/(df$reads_all))>=0.95, "KF",
          ifelse(df$type_c=='d' & df$reads_all>=threshold_i & df$c<2 & (df$k/(df$reads_all))>=0.95, "KP",
#if there aren't any deletion reads, it's probably SK1 (unless read depth is low).
#But need to check that the majority of reads match the reference as well, else it might be a variant table error
#SK1 hits are only evaluated if insertions =0 or 1
          ifelse(df$type_c=='s' & df$reads_all<threshold_s & (df$c/(df$reads_all))>=0.75, "CF",
          ifelse(df$type_c=='s' & df$reads_all<threshold_s & (df$k/(df$reads_all))>=0.9, "KF",
          ifelse(df$type_c=='s' & df$reads_all>=threshold_s & (df$c/(df$reads_all))>=0.75, "CP",
          ifelse(df$type_c=='s' & df$reads_all>=threshold_s & (df$k/(df$reads_all))>=0.9, "KP",
          ifelse((df$reads_all)<threshold_s & ((df$SNPreads)/df$reads_all)>0.9 & (df$c/(df$SNPreads))<0.7 & (df$k/(df$SNPreads))<0.7, "HF",
          ifelse((df$reads_all)>=threshold_s & ((df$SNPreads)/df$reads_all)>0.9 & (df$c/(df$SNPreads))<0.7 & (df$k/(df$SNPreads))<0.7, "HP",
#heteroduplex calls for all variant types. Note, for indels, these should NOT be later converted - only used for statistical purposes.
#note, it is important to do the H calls before the N calls because the H calls also fit the criteria for N calls
#And also to do the H calls after the deletion calls, because some deletion bases get called HP which results in an overall U call, allowing the insertion to dominate
#insertions never get called H with these criteria, (because insertions don't count seperately to 'reads_all'?)
          ifelse(df$type_c=='i' & df$reads_all<threshold_i & (df$k/(df$reads_all))<0.95 & (df$c/(df$reads_all))<0.3, "NF",
          ifelse(df$type_c=='i' & df$reads_all>=threshold_i & (df$k/(df$reads_all))<0.95 & (df$c/(df$reads_all))<0.3, "NP",  
          ifelse(df$type_c=='d' & df$reads_all<threshold_i & (df$k/(df$reads_all))<0.95 & (df$c/(df$reads_all))<0.3, "NF",
          ifelse(df$type_c=='d' & df$reads_all>=threshold_i & (df$k/(df$reads_all))<0.95 & (df$c/(df$reads_all))<0.3, "NP", 
          ifelse(df$type_c=='s' & df$reads_all<threshold_s & (df$k/(df$reads_all))<0.9 & (df$c/(df$reads_all))<0.75, "NF",
          ifelse(df$type_c=='s' & df$reads_all>=threshold_s & (df$k/(df$reads_all))<0.9 & (df$c/(df$reads_all))<0.75, "NP",  
          "U"))))))))))))))))))))))))
  
  ###Collapsing deletions###
  df_before_collapse <- df
  dels <- df[which(df$type_c=='d'),]
  delsFB <-dels[which(dels$seq_c!="-"),] #keep only the first base of each deletion
  delsOB <-dels[which(dels$seq_c=="-"),] #the other bases of the deletions
  
  #1bp = 100% (1 call)
  #5bp = 40% (2 calls)
  #10bp= 33% (3 calls)
  #20bp= 25% (5 calls)
  #50bp= 16% (8 calls)
  #100bp=10% (10 calls)
  #100+ =5%
  
  # Save delsFB to a new data frame so results can be compared afterward. 
  delsFB_new <- delsFB
  
  # Loop through the unique uID values in delsOB.
  for (unique_id in unique(delsFB$uID)) {   #rarely, there may be occurences of uIDs in delsOB that are not in delsFB, which causes an error when combining the tables. So only run the loop for uIDs in delsFB.
    # Subset delsOB for the current unique uID.
    current_delsOB_subset <- delsOB[delsOB$uID == unique_id, ]
    current_length <- current_delsOB_subset[1,10]-1  #store length of current indel for working out the average reads
    # Assign desired value to result based on the conditions.
    totalc = sum(current_delsOB_subset$c) #total c reads
    totalk = sum(current_delsOB_subset$k) #total k reads
    avgc <- totalc/current_length #average c reads
    avgk <- totalk/current_length #average k reads
    # Assign desired value to result based on the conditions.
    if ("KP" %in% current_delsOB_subset$kb & "CP" %in% current_delsOB_subset$kb) {
      result <- "HP"
    }
    else if ("KP" %in% current_delsOB_subset$kb) {
      # Nested conditional logic to ensure a minimum proportion of "KP" values, according to the number of instances.
      #Only if the number of KPs or CPs >0, or it will give /0=infinite
      if (length(current_delsOB_subset$kb) == 1 & length(current_delsOB_subset$kb[current_delsOB_subset$kb == "KP"]) == 1) {
        result <- "KP"
      }
      else if (length(current_delsOB_subset$kb) > 1 & length(current_delsOB_subset$kb) <= 5 & length(current_delsOB_subset$kb[current_delsOB_subset$kb == "KP"]) >= 1 & length(current_delsOB_subset$kb[current_delsOB_subset$kb == "KP"])/length(current_delsOB_subset$kb) >= 0.4) {
        result <- "KP"
      }
      else if (length(current_delsOB_subset$kb) > 5 & length(current_delsOB_subset$kb) <= 10 & length(current_delsOB_subset$kb[current_delsOB_subset$kb == "KP"]) >= 1 & length(current_delsOB_subset$kb[current_delsOB_subset$kb == "KP"])/length(current_delsOB_subset$kb) >= 0.33) {
        result <- "KP"
      }
      else if (length(current_delsOB_subset$kb) > 10 & length(current_delsOB_subset$kb) <= 20 & length(current_delsOB_subset$kb[current_delsOB_subset$kb == "KP"]) >= 1 & length(current_delsOB_subset$kb[current_delsOB_subset$kb == "KP"])/length(current_delsOB_subset$kb) >= 0.25) {
        result <- "KP"
      }
      else if (length(current_delsOB_subset$kb) > 20 & length(current_delsOB_subset$kb) <= 50 & length(current_delsOB_subset$kb[current_delsOB_subset$kb == "KP"]) >= 1 & length(current_delsOB_subset$kb[current_delsOB_subset$kb == "KP"])/length(current_delsOB_subset$kb) >= 0.16) {
        result <- "KP"
      }
      else if (length(current_delsOB_subset$kb) > 50 & length(current_delsOB_subset$kb) <= 100 & length(current_delsOB_subset$kb[current_delsOB_subset$kb == "KP"]) >= 1 & length(current_delsOB_subset$kb[current_delsOB_subset$kb == "KP"])/length(current_delsOB_subset$kb) >= 0.1) {
        result <- "KP"
      }
      else if (length(current_delsOB_subset$kb) > 100 & length(current_delsOB_subset$kb[current_delsOB_subset$kb == "KP"]) >= 1 & length(current_delsOB_subset$kb[current_delsOB_subset$kb == "KP"])/length(current_delsOB_subset$kb) >= 0.05) {
        result <- "KP"
      }
    }
    else if ("CP" %in% current_delsOB_subset$kb) {
      # Nested conditional logic to ensure a minimum proportion of "CP" values, according to the number of instances.
      #Only if the number of KPs or CPs >0, or it will give /0=infinite
      if (length(current_delsOB_subset$kb) == 1 & length(current_delsOB_subset$kb[current_delsOB_subset$kb == "CP"]) == 1) {
        result <- "CP"
      }
      else if (length(current_delsOB_subset$kb) > 1 & length(current_delsOB_subset$kb) <= 5 & length(current_delsOB_subset$kb[current_delsOB_subset$kb == "CP"]) >= 1 & length(current_delsOB_subset$kb[current_delsOB_subset$kb == "CP"])/length(current_delsOB_subset$kb) >= 0.4) {
        result <- "CP"
      }
      else if (length(current_delsOB_subset$kb) > 5 & length(current_delsOB_subset$kb) <= 10 & length(current_delsOB_subset$kb[current_delsOB_subset$kb == "CP"]) >= 1 & length(current_delsOB_subset$kb[current_delsOB_subset$kb == "CP"])/length(current_delsOB_subset$kb) >= 0.33) {
        result <- "CP"
      }
      else if (length(current_delsOB_subset$kb) > 10 & length(current_delsOB_subset$kb) <= 20 & length(current_delsOB_subset$kb[current_delsOB_subset$kb == "CP"]) >= 1 & length(current_delsOB_subset$kb[current_delsOB_subset$kb == "CP"])/length(current_delsOB_subset$kb) >= 0.25) {
        result <- "CP"
      }
      else if (length(current_delsOB_subset$kb) > 20 & length(current_delsOB_subset$kb) <= 50 & length(current_delsOB_subset$kb[current_delsOB_subset$kb == "CP"]) >= 1 & length(current_delsOB_subset$kb[current_delsOB_subset$kb == "CP"])/length(current_delsOB_subset$kb) >= 0.16) {
        result <- "CP"
      }
      else if (length(current_delsOB_subset$kb) > 50 & length(current_delsOB_subset$kb) <= 100 & length(current_delsOB_subset$kb[current_delsOB_subset$kb == "CP"]) >= 1 & length(current_delsOB_subset$kb[current_delsOB_subset$kb == "CP"])/length(current_delsOB_subset$kb) >= 0.1) {
        result <- "CP"
      }
      else if (length(current_delsOB_subset$kb) > 100 & length(current_delsOB_subset$kb[current_delsOB_subset$kb == "CP"]) >= 1 & length(current_delsOB_subset$kb[current_delsOB_subset$kb == "CP"])/length(current_delsOB_subset$kb) >= 0.05) {
        result <- "CP"
      }
    }
    else {
      result <- "U"
    }
    # Insert the result into delsFB$kb.
    delsFB_new[delsFB_new$uID == unique_id, ]$kb <- result
    delsFB_new[delsFB_new$uID == unique_id, ]$c <- avgc #replace the c reads for the first base with the average c reads from the other bases
    delsFB_new[delsFB_new$uID == unique_id, ]$k <- avgk
  }
  

  #put delsFB rows back into dataframe#
  df <- df[which(df$type_c!='d'),] #remove old del entries
  df <- rbind(df, delsFB_new) #adding the converted dels back in
  df <- df[order(df$chrom, df$pos_k),] #sort back into correct order
  
  sub<- subset(df, select = c(uID, sID, chrom, pos_k, pos_c, type_c, type_k, Var_len, c, k, kb))
  write.table(sub,file=sprintf("%sCalled_SK1.txt",filename),row.names=F, quote=F, sep='\t')
  
  ##Make a dataframe to hold statistics##
  S_reads<- df[which(df$type_c=='s'),]
  I_reads<- df[which(df$type_c=='i'),]
  D_reads_after <- df[which(df$type_c=='d'),]
  D_reads_before <-df_before_collapse[which(df_before_collapse$type_c=='d'),]
  
  statistic <- c('Median_SNP_reads','Mean_SNP_reads','Total_SNP_reads','Min_SNP_Reads','Max_SNP_reads', 'read_depth_threshold', 'Percent_failing_read_depth_threshold', 'Number_failing_read_depth_threshold', 'Number_CF', 'Number_CP', 'Number_KF', 'Number_KP', 'Number_NF', 'Number_NP', 'Number_HF', 'Number_HP', 'Percent_CF', 'Percent_CP', 'Percent_KF', 'Percent_KP', 'Percent_NF', 'Percent_NP', 'Percent_HF', 'Percent_HP', 'Total_C', 'Total_K', 'Total_N', 'Total_H', 'Percent_C', 'Percent_K', 'Percent_N', 'Percent_H', 'Number_U', 'Percent_U')
  value <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  
  statstable <- data.frame(statistic, value)
  #Stats about the number of reads - Median_SNP_reads,Mean_SNP_reads, Total_SNP_reads, Min_SNP_Reads, 
  #Max_SNP_reads, read_depth_threshold, Percent_failing_read_depth_threshold, Number_failing_read_depth_threshold
  statstable[1,2] <- median(S_reads$reads_all)
  statstable[2,2] <- mean(S_reads$reads_all)
  statstable[3,2] <- sum(S_reads$reads_all)
  statstable[4,2] <- min(S_reads$reads_all)
  statstable[5,2] <- max(S_reads$reads_all)
  statstable[6,2] <- threshold_s
  statstable[7,2] <- 100*(sum(S_reads$reads_all < threshold_s)/sum(S_reads$reads_all))
  statstable[8,2] <- sum(S_reads$reads_all<threshold_s)
  
  #More stats - 'Number_CF', 'Number_CP', 'Number_KF', 'Number_KP', 'Number_NF', 'Number_NP', 'Number_HF', 'Number_HP'
  statstable[9,2] <- sum(S_reads$kb=='CF')
  statstable[10,2] <- sum(S_reads$kb=='CP')
  statstable[11,2] <- sum(S_reads$kb=='KF')
  statstable[12,2] <- sum(S_reads$kb=='KP')
  statstable[13,2] <- sum(S_reads$kb=='NF')
  statstable[14,2] <- sum(S_reads$kb=='NP')
  statstable[15,2] <- sum(S_reads$kb=='HF')
  statstable[16,2] <- sum(S_reads$kb=='HP')
  
  #'Percent_CF', 'Percent_CP', 'Percent_KF', 'Percent_KP', 'Percent_NF', 'Percent_NP', 'Percent_HF', 'Percent_HP'
  statstable[17,2] <-100*(statstable[9,2]/(statstable[9,2]+statstable[10,2]+statstable[11,2]+statstable[12,2]+statstable[13,2]+statstable[14,2]+statstable[15,2]+statstable[16,2])) 
  statstable[18,2] <-100*(statstable[10,2]/(statstable[9,2]+statstable[10,2]+statstable[11,2]+statstable[12,2]+statstable[13,2]+statstable[14,2]+statstable[15,2]+statstable[16,2])) 
  statstable[19,2] <-100*(statstable[11,2]/(statstable[9,2]+statstable[10,2]+statstable[11,2]+statstable[12,2]+statstable[13,2]+statstable[14,2]+statstable[15,2]+statstable[16,2])) 
  statstable[20,2] <-100*(statstable[12,2]/(statstable[9,2]+statstable[10,2]+statstable[11,2]+statstable[12,2]+statstable[13,2]+statstable[14,2]+statstable[15,2]+statstable[16,2])) 
  statstable[21,2] <-100*(statstable[13,2]/(statstable[9,2]+statstable[10,2]+statstable[11,2]+statstable[12,2]+statstable[13,2]+statstable[14,2]+statstable[15,2]+statstable[16,2])) 
  statstable[22,2] <-100*(statstable[14,2]/(statstable[9,2]+statstable[10,2]+statstable[11,2]+statstable[12,2]+statstable[13,2]+statstable[14,2]+statstable[15,2]+statstable[16,2])) 
  statstable[23,2] <-100*(statstable[15,2]/(statstable[9,2]+statstable[10,2]+statstable[11,2]+statstable[12,2]+statstable[13,2]+statstable[14,2]+statstable[15,2]+statstable[16,2])) 
  statstable[24,2] <-100*(statstable[16,2]/(statstable[9,2]+statstable[10,2]+statstable[11,2]+statstable[12,2]+statstable[13,2]+statstable[14,2]+statstable[15,2]+statstable[16,2])) 
  
  #'Total_C', 'Total_K', 'Total_N', 'Total_H', 'Percent_C', 'Percent_K', 'Percent_N', 'Percent_H'
  statstable[25,2] <- (statstable[9,2] + statstable[10,2])
  statstable[26,2] <- (statstable[11,2] + statstable[12,2])
  statstable[27,2] <- (statstable[13,2] + statstable[14,2])
  statstable[28,2] <- (statstable[15,2] + statstable[16,2])
  statstable[33,2] <- sum(S_reads$kb=='U')
  statstable[29,2] <- 100*(statstable[25,2]/(statstable[25,2] + statstable[26,2]+ statstable[27,2]+ statstable[28,2]+statstable[33,2]))
  statstable[30,2] <- 100*(statstable[26,2]/(statstable[25,2] + statstable[26,2]+ statstable[27,2]+ statstable[28,2]+statstable[33,2]))
  statstable[31,2] <- 100*(statstable[27,2]/(statstable[25,2] + statstable[26,2]+ statstable[27,2]+ statstable[28,2]+statstable[33,2]))
  statstable[32,2] <- 100*(statstable[28,2]/(statstable[25,2] + statstable[26,2]+ statstable[27,2]+ statstable[28,2]+statstable[33,2]))
  statstable[34,2] <- 100*(statstable[33,2]/(statstable[25,2] + statstable[26,2]+ statstable[27,2]+ statstable[28,2]+statstable[33,2]))
  
  statstable$value <- round(statstable$value, digits = 3)
  
  statstable <- rename(statstable, c("value"=paste(match.call()[3])))
  write.table(statstable,file=sprintf("%sSNP_Stats_SK1.txt",filename),row.names=F, quote=F, sep='\t')
  
  ###############################################################################################
  #Insertion stats#
  
  statistic2 <- c('Median_insertion_reads','Mean_insertion_reads','Total_insertion_reads','Min_insertion_Reads','Max_insertion_reads', 'read_depth_threshold', 'Percent_failing_read_depth_threshold', 'Number_failing_read_depth_threshold', 'Number_CF', 'Number_CP', 'Number_KF', 'Number_KP', 'Number_NF', 'Number_NP', 'Number_HF', 'Number_HP', 'Percent_CF', 'Percent_CP', 'Percent_KF', 'Percent_KP', 'Percent_NF', 'Percent_NP', 'Percent_HF', 'Percent_HP', 'Total_C', 'Total_K', 'Total_N', 'Total_H', 'Percent_C', 'Percent_K', 'Percent_N', 'Percent_H', 'Number_U', 'Percent_U')
  value2 <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  
  statstable2 <- data.frame(statistic2, value2)
  #Stats about the number of reads - Median_SNP_reads,Mean_SNP_reads, Total_SNP_reads, Min_SNP_Reads, 
  #Max_SNP_reads, read_depth_threshold, Percent_failing_read_depth_threshold, Number_failing_read_depth_threshold
  statstable2[1,2] <- median(I_reads$reads_all)
  statstable2[2,2] <- mean(I_reads$reads_all)
  statstable2[3,2] <- sum(I_reads$reads_all)
  statstable2[4,2] <- min(I_reads$reads_all)
  statstable2[5,2] <- max(I_reads$reads_all)
  statstable2[6,2] <- threshold_i
  statstable2[7,2] <- 100*(sum(I_reads$reads_all < threshold_i)/sum(I_reads$reads_all))
  statstable2[8,2] <- sum(I_reads$reads_all<threshold_i)
  
  #More stats - 'Number_CF', 'Number_CP', 'Number_KF', 'Number_KP', 'Number_NF', 'Number_NP', 'Number_HF', 'Number_HP'
  statstable2[9,2] <- sum(I_reads$kb=='CF')
  statstable2[10,2] <- sum(I_reads$kb=='CP')
  statstable2[11,2] <- sum(I_reads$kb=='KF')
  statstable2[12,2] <- sum(I_reads$kb=='KP')
  statstable2[13,2] <- sum(I_reads$kb=='NF')
  statstable2[14,2] <- sum(I_reads$kb=='NP')
  statstable2[15,2] <- sum(I_reads$kb=='HF')
  statstable2[16,2] <- sum(I_reads$kb=='HP')
  
  #'Percent_CF', 'Percent_CP', 'Percent_KF', 'Percent_KP', 'Percent_NF', 'Percent_NP', 'Percent_HF', 'Percent_HP'
  statstable2[17,2] <-100*(statstable2[9,2]/(statstable2[9,2]+statstable2[10,2]+statstable2[11,2]+statstable2[12,2]+statstable2[13,2]+statstable2[14,2]+statstable2[15,2]+statstable2[16,2])) 
  statstable2[18,2] <-100*(statstable2[10,2]/(statstable2[9,2]+statstable2[10,2]+statstable2[11,2]+statstable2[12,2]+statstable2[13,2]+statstable2[14,2]+statstable2[15,2]+statstable2[16,2])) 
  statstable2[19,2] <-100*(statstable2[11,2]/(statstable2[9,2]+statstable2[10,2]+statstable2[11,2]+statstable2[12,2]+statstable2[13,2]+statstable2[14,2]+statstable2[15,2]+statstable2[16,2])) 
  statstable2[20,2] <-100*(statstable2[12,2]/(statstable2[9,2]+statstable2[10,2]+statstable2[11,2]+statstable2[12,2]+statstable2[13,2]+statstable2[14,2]+statstable2[15,2]+statstable2[16,2])) 
  statstable2[21,2] <-100*(statstable2[13,2]/(statstable2[9,2]+statstable2[10,2]+statstable2[11,2]+statstable2[12,2]+statstable2[13,2]+statstable2[14,2]+statstable2[15,2]+statstable2[16,2])) 
  statstable2[22,2] <-100*(statstable2[14,2]/(statstable2[9,2]+statstable2[10,2]+statstable2[11,2]+statstable2[12,2]+statstable2[13,2]+statstable2[14,2]+statstable2[15,2]+statstable2[16,2])) 
  statstable2[23,2] <-100*(statstable2[15,2]/(statstable2[9,2]+statstable2[10,2]+statstable2[11,2]+statstable2[12,2]+statstable2[13,2]+statstable2[14,2]+statstable2[15,2]+statstable2[16,2])) 
  statstable2[24,2] <-100*(statstable2[16,2]/(statstable2[9,2]+statstable2[10,2]+statstable2[11,2]+statstable2[12,2]+statstable2[13,2]+statstable2[14,2]+statstable2[15,2]+statstable2[16,2])) 
  
  #'Total_C', 'Total_K', 'Total_N', 'Total_H', 'Percent_C', 'Percent_K', 'Percent_N', 'Percent_H'
  statstable2[25,2] <- (statstable2[9,2] + statstable2[10,2])
  statstable2[26,2] <- (statstable2[11,2] + statstable2[12,2])
  statstable2[27,2] <- (statstable2[13,2] + statstable2[14,2])
  statstable2[28,2] <- (statstable2[15,2] + statstable2[16,2])
  statstable2[33,2] <- sum(I_reads$kb=='U')
  statstable2[29,2] <- 100*(statstable2[25,2]/(statstable2[25,2] + statstable2[26,2]+ statstable2[27,2]+ statstable2[28,2]+statstable2[33,2]))
  statstable2[30,2] <- 100*(statstable2[26,2]/(statstable2[25,2] + statstable2[26,2]+ statstable2[27,2]+ statstable2[28,2]+statstable2[33,2]))
  statstable2[31,2] <- 100*(statstable2[27,2]/(statstable2[25,2] + statstable2[26,2]+ statstable2[27,2]+ statstable2[28,2]+statstable2[33,2]))
  statstable2[32,2] <- 100*(statstable2[28,2]/(statstable2[25,2] + statstable2[26,2]+ statstable2[27,2]+ statstable2[28,2]+statstable2[33,2]))
  statstable2[34,2] <- 100*(statstable2[33,2]/(statstable2[25,2] + statstable2[26,2]+ statstable2[27,2]+ statstable2[28,2]+statstable2[33,2]))
  
  statstable2$value2 <- round(statstable2$value2, digits = 3)
  
  statstable2 <- rename(statstable2, c("value2"=paste(match.call()[3])))
  #rename this column so it will have the sample name for merging
  write.table(statstable2,file=sprintf("%sinsertion_Stats_SK1.txt",filename),row.names=F, quote=F, sep='\t')
  
  ###############################################################################################
  #Deletion stats - After Collapsing#
  
  statistic3 <- c('Median_deletion_reads','Mean_deletion_reads','Total_deletion_reads','Min_deletion_Reads','Max_deletion_reads', 'read_depth_threshold', 'Percent_failing_read_depth_threshold', 'Number_failing_read_depth_threshold', 'Number_CF', 'Number_CP', 'Number_KF', 'Number_KP', 'Number_NF', 'Number_NP', 'Number_HF', 'Number_HP', 'Percent_CF', 'Percent_CP', 'Percent_KF', 'Percent_KP', 'Percent_NF', 'Percent_NP', 'Percent_HF', 'Percent_HP', 'Total_C', 'Total_K', 'Total_N', 'Total_H', 'Percent_C', 'Percent_K', 'Percent_N', 'Percent_H', 'Number_U', 'Percent_U')
  value3 <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  
  statstable3 <- data.frame(statistic3, value3)
  #Stats about the number of reads - Median_SNP_reads,Mean_SNP_reads, Total_SNP_reads, Min_SNP_Reads, 
  #Max_SNP_reads, read_depth_threshold, Percent_failing_read_depth_threshold, Number_failing_read_depth_threshold
  statstable3[1,2] <- median(D_reads_after$reads_all)
  statstable3[2,2] <- mean(D_reads_after$reads_all)
  statstable3[3,2] <- sum(D_reads_after$reads_all)
  statstable3[4,2] <- min(D_reads_after$reads_all)
  statstable3[5,2] <- max(D_reads_after$reads_all)
  statstable3[6,2] <- threshold_i
  statstable3[7,2] <- 100*(sum(D_reads_after$reads_all < threshold_i)/sum(D_reads_after$reads_all))
  statstable3[8,2] <- sum(D_reads_after$reads_all<threshold_i)
  
  #More stats - 'Number_CF', 'Number_CP', 'Number_KF', 'Number_KP', 'Number_NF', 'Number_NP', 'Number_HF', 'Number_HP'
  statstable3[9,2] <- sum(D_reads_after$kb=='CF')
  statstable3[10,2] <- sum(D_reads_after$kb=='CP')
  statstable3[11,2] <- sum(D_reads_after$kb=='KF')
  statstable3[12,2] <- sum(D_reads_after$kb=='KP')
  statstable3[13,2] <- sum(D_reads_after$kb=='NF')
  statstable3[14,2] <- sum(D_reads_after$kb=='NP')
  statstable3[15,2] <- sum(D_reads_after$kb=='HF')
  statstable3[16,2] <- sum(D_reads_after$kb=='HP')
  
  #'Percent_CF', 'Percent_CP', 'Percent_KF', 'Percent_KP', 'Percent_NF', 'Percent_NP', 'Percent_HF', 'Percent_HP'
  statstable3[17,2] <-100*(statstable3[9,2]/(statstable3[9,2]+statstable3[10,2]+statstable3[11,2]+statstable3[12,2]+statstable3[13,2]+statstable3[14,2]+statstable3[15,2]+statstable3[16,2])) 
  statstable3[18,2] <-100*(statstable3[10,2]/(statstable3[9,2]+statstable3[10,2]+statstable3[11,2]+statstable3[12,2]+statstable3[13,2]+statstable3[14,2]+statstable3[15,2]+statstable3[16,2])) 
  statstable3[19,2] <-100*(statstable3[11,2]/(statstable3[9,2]+statstable3[10,2]+statstable3[11,2]+statstable3[12,2]+statstable3[13,2]+statstable3[14,2]+statstable3[15,2]+statstable3[16,2])) 
  statstable3[20,2] <-100*(statstable3[12,2]/(statstable3[9,2]+statstable3[10,2]+statstable3[11,2]+statstable3[12,2]+statstable3[13,2]+statstable3[14,2]+statstable3[15,2]+statstable3[16,2])) 
  statstable3[21,2] <-100*(statstable3[13,2]/(statstable3[9,2]+statstable3[10,2]+statstable3[11,2]+statstable3[12,2]+statstable3[13,2]+statstable3[14,2]+statstable3[15,2]+statstable3[16,2])) 
  statstable3[22,2] <-100*(statstable3[14,2]/(statstable3[9,2]+statstable3[10,2]+statstable3[11,2]+statstable3[12,2]+statstable3[13,2]+statstable3[14,2]+statstable3[15,2]+statstable3[16,2])) 
  statstable3[23,2] <-100*(statstable3[15,2]/(statstable3[9,2]+statstable3[10,2]+statstable3[11,2]+statstable3[12,2]+statstable3[13,2]+statstable3[14,2]+statstable3[15,2]+statstable3[16,2])) 
  statstable3[24,2] <-100*(statstable3[16,2]/(statstable3[9,2]+statstable3[10,2]+statstable3[11,2]+statstable3[12,2]+statstable3[13,2]+statstable3[14,2]+statstable3[15,2]+statstable3[16,2])) 
  
  #'Total_C', 'Total_K', 'Total_N', 'Total_H', 'Percent_C', 'Percent_K', 'Percent_N', 'Percent_H'
  statstable3[25,2] <- (statstable3[9,2] + statstable3[10,2])
  statstable3[26,2] <- (statstable3[11,2] + statstable3[12,2])
  statstable3[27,2] <- (statstable3[13,2] + statstable3[14,2])
  statstable3[28,2] <- (statstable3[15,2] + statstable3[16,2])
  statstable[33,2] <- sum(D_reads_after$kb=='U')
  statstable3[29,2] <- 100*(statstable3[25,2]/(statstable3[25,2] + statstable3[26,2]+ statstable3[27,2]+ statstable3[28,2]+statstable3[33,2]))
  statstable3[30,2] <- 100*(statstable3[26,2]/(statstable3[25,2] + statstable3[26,2]+ statstable3[27,2]+ statstable3[28,2]+statstable3[33,2]))
  statstable3[31,2] <- 100*(statstable3[27,2]/(statstable3[25,2] + statstable3[26,2]+ statstable3[27,2]+ statstable3[28,2]+statstable3[33,2]))
  statstable3[32,2] <- 100*(statstable3[28,2]/(statstable3[25,2] + statstable3[26,2]+ statstable3[27,2]+ statstable3[28,2]+statstable3[33,2]))
  statstable3[34,2] <- 100*(statstable3[33,2]/(statstable3[25,2] + statstable3[26,2]+ statstable3[27,2]+ statstable3[28,2]+statstable3[33,2]))
  
  statstable3$value3 <- round(statstable3$value3, digits = 3)
  
  statstable3 <- rename(statstable3, c("value3"=paste(match.call()[3])))
  #rename this column so it will have the sample name for merging
  write.table(statstable3,file=sprintf("%sdeletion_Stats_after_collapsing_SK1.txt",filename),row.names=F, quote=F, sep='\t')
  
  ###############################################################################################
  #Deletion stats - Before Collapsing#
  
  statistic4 <- c('Median_deletion_reads','Mean_deletion_reads','Total_deletion_reads','Min_deletion_Reads','Max_deletion_reads', 'read_depth_threshold', 'Percent_failing_read_depth_threshold', 'Number_failing_read_depth_threshold', 'Number_CF', 'Number_CP', 'Number_KF', 'Number_KP', 'Number_NF', 'Number_NP', 'Number_HF', 'Number_HP', 'Percent_CF', 'Percent_CP', 'Percent_KF', 'Percent_KP', 'Percent_NF', 'Percent_NP', 'Percent_HF', 'Percent_HP', 'Total_C', 'Total_K', 'Total_N', 'Total_H', 'Percent_C', 'Percent_K', 'Percent_N', 'Percent_H', 'Number_U', 'Percent_U')
  value4 <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  
  statstable4 <- data.frame(statistic4, value4)
  #Stats about the number of reads - Median_SNP_reads,Mean_SNP_reads, Total_SNP_reads, Min_SNP_Reads, 
  #Max_SNP_reads, read_depth_threshold, Percent_failing_read_depth_threshold, Number_failing_read_depth_threshold
  statstable4[1,2] <- median(D_reads_before$reads_all)
  statstable4[2,2] <- mean(D_reads_before$reads_all)
  statstable4[3,2] <- sum(D_reads_before$reads_all)
  statstable4[4,2] <- min(D_reads_before$reads_all)
  statstable4[5,2] <- max(D_reads_before$reads_all)
  statstable4[6,2] <- threshold_i
  statstable4[7,2] <- 100*(sum(D_reads_before$reads_all < threshold_i)/sum(D_reads_before$reads_all))
  statstable4[8,2] <- sum(D_reads_before$reads_all<threshold_i)
  
  #More stats - 'Number_CF', 'Number_CP', 'Number_KF', 'Number_KP', 'Number_NF', 'Number_NP', 'Number_HF', 'Number_HP'
  statstable4[9,2] <- sum(D_reads_before$kb=='CF')
  statstable4[10,2] <- sum(D_reads_before$kb=='CP')
  statstable4[11,2] <- sum(D_reads_before$kb=='KF')
  statstable4[12,2] <- sum(D_reads_before$kb=='KP')
  statstable4[13,2] <- sum(D_reads_before$kb=='NF')
  statstable4[14,2] <- sum(D_reads_before$kb=='NP')
  statstable4[15,2] <- sum(D_reads_before$kb=='HF')
  statstable4[16,2] <- sum(D_reads_before$kb=='HP')
  
  #'Percent_CF', 'Percent_CP', 'Percent_KF', 'Percent_KP', 'Percent_NF', 'Percent_NP', 'Percent_HF', 'Percent_HP'
  statstable4[17,2] <-100*(statstable4[9,2]/(statstable4[9,2]+statstable4[10,2]+statstable4[11,2]+statstable4[12,2]+statstable4[13,2]+statstable4[14,2]+statstable4[15,2]+statstable4[16,2])) 
  statstable4[18,2] <-100*(statstable4[10,2]/(statstable4[9,2]+statstable4[10,2]+statstable4[11,2]+statstable4[12,2]+statstable4[13,2]+statstable4[14,2]+statstable4[15,2]+statstable4[16,2])) 
  statstable4[19,2] <-100*(statstable4[11,2]/(statstable4[9,2]+statstable4[10,2]+statstable4[11,2]+statstable4[12,2]+statstable4[13,2]+statstable4[14,2]+statstable4[15,2]+statstable4[16,2])) 
  statstable4[20,2] <-100*(statstable4[12,2]/(statstable4[9,2]+statstable4[10,2]+statstable4[11,2]+statstable4[12,2]+statstable4[13,2]+statstable4[14,2]+statstable4[15,2]+statstable4[16,2])) 
  statstable4[21,2] <-100*(statstable4[13,2]/(statstable4[9,2]+statstable4[10,2]+statstable4[11,2]+statstable4[12,2]+statstable4[13,2]+statstable4[14,2]+statstable4[15,2]+statstable4[16,2])) 
  statstable4[22,2] <-100*(statstable4[14,2]/(statstable4[9,2]+statstable4[10,2]+statstable4[11,2]+statstable4[12,2]+statstable4[13,2]+statstable4[14,2]+statstable4[15,2]+statstable4[16,2])) 
  statstable4[23,2] <-100*(statstable4[15,2]/(statstable4[9,2]+statstable4[10,2]+statstable4[11,2]+statstable4[12,2]+statstable4[13,2]+statstable4[14,2]+statstable4[15,2]+statstable4[16,2])) 
  statstable4[24,2] <-100*(statstable4[16,2]/(statstable4[9,2]+statstable4[10,2]+statstable4[11,2]+statstable4[12,2]+statstable4[13,2]+statstable4[14,2]+statstable4[15,2]+statstable4[16,2])) 
  
  #'Total_C', 'Total_K', 'Total_N', 'Total_H', 'Percent_C', 'Percent_K', 'Percent_N', 'Percent_H'
  statstable4[25,2] <- (statstable4[9,2] + statstable4[10,2])
  statstable4[26,2] <- (statstable4[11,2] + statstable4[12,2])
  statstable4[27,2] <- (statstable4[13,2] + statstable4[14,2])
  statstable4[28,2] <- (statstable4[15,2] + statstable4[16,2])
  statstable4[33,2] <- sum(D_reads_before$kb=='U')
  statstable4[29,2] <- 100*(statstable4[25,2]/(statstable4[25,2] + statstable4[26,2]+ statstable4[27,2]+ statstable4[28,2]+statstable4[33,2]))
  statstable4[30,2] <- 100*(statstable4[26,2]/(statstable4[25,2] + statstable4[26,2]+ statstable4[27,2]+ statstable4[28,2]+statstable4[33,2]))
  statstable4[31,2] <- 100*(statstable4[27,2]/(statstable4[25,2] + statstable4[26,2]+ statstable4[27,2]+ statstable4[28,2]+statstable4[33,2]))
  statstable4[32,2] <- 100*(statstable4[28,2]/(statstable4[25,2] + statstable4[26,2]+ statstable4[27,2]+ statstable4[28,2]+statstable4[33,2]))
  statstable4[34,2] <- 100*(statstable4[33,2]/(statstable4[25,2] + statstable4[26,2]+ statstable4[27,2]+ statstable4[28,2]+statstable4[33,2]))
  
  statstable4$value4 <- round(statstable4$value4, digits = 3)
  
  statstable4 <- rename(statstable4, c("value4"=paste(match.call()[3])))
  #rename this column so it will have the sample name for merging
  write.table(statstable4,file=sprintf("%sdeletion_Stats_before_collapsing_SK1.txt",filename),row.names=F, quote=F, sep='\t')
  
}

