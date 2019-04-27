
a = open("OMN6_binary.txt").read().split("\n")

out = open("Annot_OMN6","w")
out.write("chr\tstart5\tstart3\tstop5\tstop3\tpattern\ttype\tnb_snp\tsporeaffect\n")

c = open("ChrLen2011_Rep3").read().split("\n") #this file just contains the length for each chromosome
chrlen = [0]
for i in range(1,len(c)-1):
        chrlen.append(c[i].split("\t")[1])

def classif(octuplet):  ##define the function 'classif' which acts on an octuplet
 

    pattern = str(sum(octuplet))+":"+str(8-sum(octuplet)) ##to work out the pattern, sum up the total of the 1 and 0 binary calls

    dif = 0
    for m in range(4):
        if octuplet[2*m] != octuplet[2*m+1]:  #when overall pattern is the same, but chromos involved are different (trans patterns)
            dif = 1
            break
            
    if pattern in ["6:2","2:6","4:4"] and dif == 1:
        pattern+="*"

    return pattern
        


for chromo in range(1,17): #start of big loop: for each chromosome,
    cpos = []
    cgen = []
    for i in range(1,len(a)-1): #subloop 1: for each line of the binary file,
        line = a[i].split("\t") #split the line into components
      # sporeaffect =line[10]
        if int(line[0]) == chromo : #if the chromosome is the one currently being worked on
            cpos.append(int(line[1])) #update cpos with pos_c from the binary file  (int: convert a string to an integer)
            gen = []
            for j in range(8):
                gen.append(int(line[2+j])) #combine the binary calls for each spore into an 8-digit number (i think)
            cgen.append(gen)


    start5 = 1  #start5 position is 1 for the first segment on the chromosome
    start3 = cpos[0] #start3 position is location of next variant
    ref = cgen[0] #reference: the pattern of the current line
    nb_snp = 1 #one SNP per line
    #line2 = a[i].split("\t")
    #sporeline =line2[10]

    ######################################################
    #########work out count here##########################
    ################################
    count = 0
    for j in range(1,len(a)-1): #subloop 2: for each line of the binary file,
        line = a[j].split("\t") #split the line into components
        if int(line[0]) < chromo : #count lines from previous chromosomes
             count += 1
   ######################################################
            
    #    fin = 1
    for m in range(1,len(cpos)-1): #subloop 2: for each base of the chromosome,
        if cgen[m] != ref : #if the pattern does not match the previous,
            stop5 = cpos[m-1] #update coordinates for end of segment
            stop3 = cpos[m]
  
            ######################################################
            q = m+count
            line = a[q].split("\t") #split the line into components
            sporeaffect =int(line[10])
            ######################################################
            out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(chromo,start5,start3,
                            stop5,stop3,ref,classif(ref),nb_snp,sporeaffect))
            debut = m
    #            fin = 0
            break
        else :
            nb_snp += 1 #if the pattern does match the previous, add another SNP to the total
           
    debut = m
    start5 = cpos[debut-1]
    start3 = cpos[debut]
    ref = cgen[debut]
    nb_snp = 1
    for x in range(debut+1,len(cpos)-1):
        if cgen[x] != ref:
            stop5 = cpos[x-1]
            stop3 = cpos[x]
            ######################################################
            r = x+count
            line = a[r].split("\t") #split the line into components
            sporeaffect =int(line[10])
            ######################################################
            out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(chromo,start5,start3,
                            stop5,stop3,ref,classif(ref),nb_snp,sporeaffect))
            start5 = stop5
            start3 = stop3
            ref = cgen[x]
            nb_snp = 1
        else :
            nb_snp += 1

            
    if cgen[len(cpos)-1] == ref :
        nb_snp += 1
        out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(chromo,start5,start3,
                    cpos[len(cpos)-1],chrlen[chromo],ref,classif(ref),nb_snp))
        ###issue: if there are no events on chromosome 1, sporeaffect is never defined for the chr1 segment?
            
    elif cgen[len(cpos)-1] != ref :
        stop5 = cpos[len(cpos)-2]
        stop3 = cpos[len(cpos)-1]
        #line = a[i].split("\t") #split the line into components
        #sporeaffect =line[10]
        out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(chromo,start5,start3,
                    stop5,stop3,ref,classif(ref),nb_snp,sporeaffect))
        start5 = cpos[len(cpos)-2]
        start3 = cpos[len(cpos)-1]
        ref = cgen[len(cpos)-1]
            
        out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(chromo,start5,start3,
                    cpos[len(cpos)-1],chrlen[chromo],ref,classif(ref),1)) #, sporeaffect

out.close()
       
