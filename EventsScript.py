#!/usr/bin/python3
# -*- coding: utf-8 -*-

a = open("NAnnot_OMN6").read().split("\n")

out = open("Events_OMN6_1500","w")
out.write("id\tchr\tstart5\tstart3\tstop5\tstop3\ttype\tCO\t\
chromatids\tlen_min\tlen_mid\tlen_max\tcommande\tgroupe\tseg_len_mid\t\
seg_len_min\tseg_len_max\tseg_start5\tseg_start3\tseg_stop5\tseg_stop3\tnb_snp\tsegnb_snp\tsegspaff\tchr_aff\n")
##add new columns for snps per segment and spores affected per segment

out3 = open("UEvents_OMN6_1500","w")
out3.write("id\tchr\tstart5\tstart3\tstop5\tstop3\ttype\tCO\t\
chromatids\tclasse\tlen_min\tlen_mid\tlen_max\tnb_snp\tcommande\n")

limite = 1500  # distance limite entre 2 events
nb_id = 1
arret = 0

for chromo in range(1,17):

    if arret : break

# Chargement des données
# Loading data

    s_start5 = [] ; s_start3 = [] ; s_stop5 = [] ; s_stop3 = []
    s_pat = [] ; s_type = [] ; s_nb_snp = [] ; s_spaff = [] ; s_chr_aff = []   # s = segment

#takes data from the Annot file
    for i in range(1,len(a)-1):
        if int(a[i].split("\t")[0]) != chromo :
            continue
        line = a[i].split("\t")
        s_start5.append(int(line[1]))
        s_start3.append(int(line[2]))
        s_stop5.append(int(line[3]))
        s_stop3.append(int(line[4]))
        pattern = line[5].replace("[","")
        pattern = pattern.replace("]","")
        pattern = pattern.split(", ")   
        seq = ""
        for c in range(8):
            seq += pattern[c]
        s_pat.append(seq)
        s_type.append(line[6])
        #Fill in the number of variants from each line
        s_nb_snp.append(int(line[7]))
        s_spaff.append(int(line[8]))
        #if there is a previous line, take from that line
        #else = 0

# Search for the first segment 4:4 chr outside the set limit:

    dernier44 = 0

    for k in range(len(s_pat)):
        if s_type[k] == "4:4" and \
           0.5*(s_stop3[k]+s_stop5[k] - (s_start3[k]+s_start5[k])) > limite :
            dernier44 = k
            break
 #dernier44 =latest 4:4

#+= Add AND : adds right operand to the left operand and assign the result to left operand

#-----------------------------------------------------------------
# Boucle d'analyse
# Analysis loop

    for k in range(dernier44+1,len(s_type)):
                
        if s_type[k] == "4:4" and \
           0.5*(s_stop3[k]+s_stop5[k] - (s_start3[k]+s_start5[k])) > limite :
            
            nb_snp = 0
            for j in range(dernier44+1,k):
                nb_snp += s_nb_snp[j]
                
            le44suivant = k
            
            #suivant: the following 4:4


# Determination of chromatid number modified

# String_pat matches the pattern of the 2 strands of one chromatid , e.g. 00
# All string_pat for all segments of the event is chromatid_pattern list .
# If there is > 1 pattern in chromatid_pattern , then the chromatid was affected.

            chr_involved = []
            chromatid_identity = []
            for chromatid in range(4):
                chromatid_pattern = []
                for segm in range(dernier44,le44suivant+1):
                    chromatid_pattern.append(s_pat[segm][2*chromatid:2*chromatid+2])
                        
                if len(set(chromatid_pattern)) > 1: #recording how many chromosomes are involved, but not their identity.
                    chr_involved.append(chromatid) 
                    
                #print('CHR:', chromatid+1)
                #print('CP:',chromatid_pattern)
                #if chromatid_pattern[0] == chromatid_pattern[-1]:
                #    print('DE')
                if chromatid_pattern[0] != chromatid_pattern[-1]:
                #    print('DNE')
                    chromatid_identity.append(chromatid+1) #record which chromatids changed in a crossover (nothing will be recorded for NCOs)
                #    print('CI:',chromatid_identity)

            nb_chromatid = len(chr_involved)

# If nb_chromatid = 2, determining whether the chromatids are involved sisters =
# Identical pattern for dernier44

            if nb_chromatid == 2 :

                if s_pat[dernier44][2*chr_involved[0]:2*chr_involved[0]+2] == \
                   s_pat[dernier44][2*chr_involved[1]:2*chr_involved[1]+2] :
                    nb_chromatid = "2_sis"
                else :
                    nb_chromatid = "2_nonsis"
                    
            #########################################       
            #s_chr_aff.append(chr_involved)
            #print(chr_involved)
            ######################################### 


# Recodage systématique de g_type pour les 5:3a, 6:2*, 4:4* et 4:4co
# Systematic Recoding g_type for 5: 3a, 6: 2 * , 4 * 4 and 4: 4co

            g_type = "" ; g_start5 = "" ; g_start3 = "" ; g_stop5 = "" ; g_stop3 = ""
            list44 = [s_pat[dernier44]]
            list44sh = [s_pat[dernier44]]  # 4:4 sans hétéroduplex
            list53 = []
            list35 = []
            list62 = []
            list26 = []
            list71 = []
            list17 = []
            list80 = []
            list08 = []
            dic_list = {"5:3":list53,"3:5":list35,"6:2":list62,"2:6":list26,
                        "7:1":list71,"1:7":list17,"4:4":list44,"8:0":list80,"0:8":list08}
            letter = ["","a","b","c","d","e","f","g","h","i","j"]
            
            for segm in range(dernier44+1,le44suivant+1):
            # Inclusion du dernier 4:4 = le44suivant mais pas du premier
            # Inclusion of the last 4: 4 = le44suivant but not the first
            
                segm_count = str(s_pat[segm].count("1")) + ":" + str(s_pat[segm].count("0"))

            # Indication des formes alternatives dans tous les cas, même 4:4
            # Indication alternative forms in all cases even 4: 4

                found = 0
                for model in range(len(dic_list[segm_count])):
                    if s_pat[segm] == dic_list[segm_count][model]:
                        segm_type = segm_count + letter[model]
                        found = 1
                        break
                if not found :
                    segm_type = segm_count + letter[len(dic_list[segm_count])]
                    dic_list[segm_count].append(s_pat[segm])

            # Indication des hDNA dans les formes paires
            # Indication of the forms in pairs HDNA
            
                if segm_count in ["4:4","6:2","2:6"]:
                    hDNA = 0
                    for chromatide in range(4):
                        if s_pat[segm][2*chromatide] != s_pat[segm][2*chromatide + 1]:
                            hDNA = 1
                            break
                    if hDNA ==1:
                        segm_type += "i"

            # Indication des 4:4co par rapport au 4:4 précédent
            # Indication of 4: 4co compared to the 4: 4 previous
            
                if segm_count == "4:4" and "i" not in segm_type \
                   and s_pat[segm] != list44sh[-1] :  # dernier élément de list44sh - # Last element of list44sh
                    segm_type += "CO"
                    list44sh.append(s_pat[segm]) # addition uniquement des 4:4 sans hDNA - addition of only 4: 4 HDNA

            # Additions à g_type, g_start, g_stop
            
                if segm == dernier44+1 : 
                    g_type += "(" + segm_type
                    g_start5 += str(s_start5[segm])
                    g_start3 += str(s_start3[segm])
                    g_stop5 += str(s_stop5[segm])
                    g_stop3 += str(s_stop3[segm])
                else :
                    g_type += ")_(" + segm_type
                    g_start5 += "_" + str(s_start5[segm])
                    g_start3 += "_" + str(s_start3[segm])
                    g_stop5 += "_" + str(s_stop5[segm])
                    g_stop3 += "_" + str(s_stop3[segm])
                

# Construction de g_len (longueur de tous les segments avec hDNA)
# Construction g_len (length of all segments with HDNA )
 
### as well as seg lens, record the no of SNPs per seg and the spores affects
            g_type += ")"
            if g_type != "(4:4aCO)" :
                g_len_mid = ""
                g_len_min = ""
                g_len_max = ""
                g_SNPs = ""
                g_spaff = ""
                seg_len_mid = 0.5*(s_stop3[dernier44+1]+s_stop5[dernier44+1]-\
                           (s_start3[dernier44+1]+s_start5[dernier44+1]))
                g_len_mid += str(round(seg_len_mid))
                
                seg_len_min = s_stop5[dernier44+1]-s_start3[dernier44+1]
                g_len_min += str(round(seg_len_min))
                
                seg_len_max = s_stop3[dernier44+1]-s_start5[dernier44+1]
                g_len_max += str(round(seg_len_max))

                seg_SNP = s_nb_snp[dernier44+1]
                g_SNPs += str(seg_SNP)

                seg_spaff = s_spaff[dernier44+1]     #[take from prev row; first occurence change to 0?]
                g_spaff += str(seg_spaff)
                
                for j in range(dernier44+2,k):
                    seg_len_mid = 0.5*(s_stop3[j]+s_stop5[j]-\
                           (s_start3[j]+s_start5[j]))
                    g_len_mid += "_" + str(round(seg_len_mid))
                    seg_len_min = s_stop5[j]-s_start3[j]
                    g_len_min += "_" + str(round(seg_len_min))
                    seg_len_max = s_stop3[j]-s_start5[j]
                    g_len_max += "_" + str(round(seg_len_max))
                    seg_SNP = s_nb_snp[j]
                    g_SNPs += "_" + str(seg_SNP)
                    seg_spaff = s_spaff[j] 
                    g_spaff += "_" + str(seg_spaff)
            else:
                g_len_mid = 0
                g_len_min = 0
                g_len_max = 0
                g_SNPs =0
                g_spaff =0
  
                            
# Determination number of CO , len_mid , len_min , len_max and order

            nb_co = g_type.count("CO")

            len_mid = round(0.5*(s_start3[le44suivant]+s_start5[le44suivant]-\
                                 (s_stop3[dernier44]+s_stop5[dernier44])))
            len_min = s_start5[le44suivant] - s_stop3[dernier44] + 1
            len_max = s_start3[le44suivant] - s_stop5[dernier44] - 1
            if g_type == "(4:4aCO)" :
                len_min = 0
                len_mid = round(0.5*len_max)
                
            ordre = ""
            if nb_chromatid != 4:
                ordre = ",c(" + str(2*chr_involved[0]+1) + "," + str(2*chr_involved[0]+2)
                if nb_chromatid != 1:
                    for inv in range(1,len(chr_involved)):
                        ordre += "," + str(2*chr_involved[inv]+1) + "," + str(2*chr_involved[inv]+2)
                ordre += ")"
                    
                    
# Determining the group if nb_chromatid == 1
            groupe = "NA"
            if nb_chromatid == 1:
                groupe = 1
                if "3:5a" in g_type or "5:3a" in g_type:
                    groupe = 2
                    

# Ecriture
            
            out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\
                      graph_bin_annot3(%s,%s,%s%s)\t%s\t%s\t\
                      %s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(
                nb_id,chromo,
                s_stop5[dernier44],s_stop3[dernier44],
                s_start5[le44suivant],s_start3[le44suivant],g_type,nb_co,
                nb_chromatid,len_min,len_mid,len_max,
                chromo,int(((s_stop5[dernier44]-1000)//1000)*1000),
                int(((s_start3[le44suivant]+1000)//1000)*1000),ordre,groupe,
                g_len_mid,g_len_min,g_len_max,g_start5,g_start3,g_stop5,g_stop3,nb_snp,g_SNPs,g_spaff,chromatid_identity))

#out.write("id\tchr\tstart5\tstart3\tstop5\tstop3\ttype\tCO\t\
#chromatids\tlen_min\tlen_mid\tlen_max\tcommande\tgroupe\tseg_len_mid\t\
#seg_len_min\tseg_len_max\tseg_start5\tseg_start3\tseg_stop5\tseg_stop3\tnb_snp\n")

# Réinitialisation des variables si nb_chromatid != 3 or nb_co != 1
# Resetting variables if nb_chromatid ! = 3 or nb_co ! = 1

            if nb_chromatid != 3 or nb_co != 1:

                nb_id += 1
                dernier44 = le44suivant
                continue
            
    
#-------------------------------------------------------------------
# Analyse des cas 3chr = 3 chromatides involved
 
            chr_in_CO = []
            chr_in_NCO = "v"  # vide
            for chromatid in chr_involved :   # soit 3
                debut = s_pat[dernier44][2*chromatid:2*chromatid+2]
                fin = s_pat[le44suivant][2*chromatid:2*chromatid+2]
                if debut != fin :
                    chr_in_CO.append(chromatid)
                else :
                    chr_in_NCO = chromatid

            if len(chr_in_CO) != 2 :    # signal d'alerte
                print("Pb pour ",nb_id)
                
# Il faut redéfinir les segments sans tenir compte des autres chromatides.

#----------------------
# Chromatide du NCO

# On commence par la chromatide avec le NCO.
# On recode d'abord le pattern de la chr in NCO.

            if s_pat[dernier44][2*chr_in_NCO:2*chr_in_NCO+2] == "00":
                autres_chr = "001111"
            else :
                autres_chr = "110000"

            nco_s_pat = 1000*["a"] # pseudo-liste
            for segm in range(dernier44,le44suivant+1):
                doublet = s_pat[segm][2*chr_in_NCO:2*chr_in_NCO+2]
                nco_s_pat[segm] = doublet+autres_chr

# On regroupe les segments identiques avec la nouvelle définition de pattern
# en incluant dernier44 = élément 0 de new_start, etc.
# et k = dernier 4:4
            
            plage = dernier44
            new_start5 = [] ; new_start3 = []
            new_stop5 = [] ; new_stop3 = []
            new_pat = []
            for segm in range(plage+1,le44suivant+1):
                ref_start5 = s_start5[plage]
                ref_start3 = s_start3[plage]
                ref_pattern = nco_s_pat[plage]

                if nco_s_pat[segm] == ref_pattern :
                    if segm == le44suivant :
                        new_start5.append(ref_start5)
                        new_start3.append(ref_start3)
                        new_stop5.append(s_stop5[segm])
                        new_stop3.append(s_stop3[segm])
                        new_pat.append(ref_pattern)

                elif nco_s_pat[segm] != ref_pattern :
                    new_start5.append(ref_start5)
                    new_start3.append(ref_start3)
                    new_stop5.append(s_stop5[segm-1])
                    new_stop3.append(s_stop3[segm-1])
                    new_pat.append(ref_pattern)
                    plage = segm
                    if segm == le44suivant :
                        new_start5.append(s_start5[segm])
                        new_start3.append(s_start3[segm])
                        new_stop5.append(s_stop5[segm])
                        new_stop3.append(s_stop3[segm])
                        new_pat.append(nco_s_pat[segm])

# Recodage systématique de g_type pour les 5:3a, 6:2*, 4:4* et 4:4co

            g_type = ""
            list44 = [new_pat[0]]
            list44sh = [new_pat[0]]  # 4:4 sans hétéroduplex
            list53 = []
            list35 = []
            list62 = []
            list26 = []
            list71 = []
            list17 = []
            list80 = []
            list08 = []
            dic_list = {"5:3":list53,"3:5":list35,"6:2":list62,"2:6":list26,
                        "7:1":list71,"1:7":list17,"4:4":list44,"8:0":list80,"0:8":list08}
            letter = ["","a","b","c","d","e","f","g","h","i","j"]
            
            for segm in range(1,len(new_pat)):
            # Inclusion du dernier 4:4 = k = le44suivant, mais pas du premier = dernier44
            
                segm_count = str(new_pat[segm].count("1")) + ":" + \
                             str(new_pat[segm].count("0"))

            # Indication des formes alternatives dans tous les cas, même 4:4                   

                found = 0
                for model in range(len(dic_list[segm_count])):
                    if new_pat[segm] == dic_list[segm_count][model]:
                        segm_type = segm_count + letter[model]
                        found = 1
                        break
                if not found :
                    segm_type = segm_count + letter[len(dic_list[segm_count])]
                    dic_list[segm_count].append(new_pat[segm])

            # Indication des hDNA dans les formes paires
            
                if segm_count in ["4:4","6:2","2:6"]:
                    hDNA = 0
                    for chromatide in range(4):
                        if new_pat[segm][2*chromatide] != new_pat[segm][2*chromatide + 1]:
                            hDNA = 1
                            break
                    if hDNA :
                        segm_type += "i"

            # Indication des 4:4co par rapport au 4:4 précédent : non pertinent dans ce cas
            
                if segm_count == "4:4" and "i" not in segm_type \
                   and new_pat[segm] != list44sh[-1] :  # dernier élément de list44sh
                    segm_type += "CO"
                    list44sh.append(new_pat[segm]) # addition uniquement des 4:4 sans hDNA

            # Addition à g_type
            
                if segm == 1 : 
                    g_type += "(" + segm_type
                else :
                    g_type += ")_(" + segm_type
                

# Construction de g_len (longueur de tous les segments en excluant les 4:4 aux extrémités)

            g_type += ")"
            if g_type != "(4:4aCO)" :
                g_len = ""
                seg_len = 0.5*(new_stop3[1]+new_stop5[1]-\
                              (new_start3[1]+new_start5[1]))
                g_len += str(round(seg_len))
                
                for j in range(2,len(new_pat)-1):
                    seg_len = 0.5*(new_stop3[j]+new_stop5[j]-\
                           (new_start3[j]+new_start5[j]))
                    g_len += "_" + str(round(seg_len))
            else:
                g_len = 0

                            
# Détermination CO/NCO et longueur totale

            nb_co = g_type.count("CO")

            longueur = round(0.5*(new_start3[-1]+new_start5[-1]-\
                                 (new_stop3[0]+new_stop5[0])))  # -1 = dernier 4:4 de new_*
            

# Détermination du groupe puisque nb_chromatid == 1

            groupe = 1
            if "3:5a" in g_type or "5:3a" in g_type:
                groupe = 2

# Ecriture de tous les paramètres sauf de overlap déterminé plus tard
# On garde debut_nco et fin_nco pour référence

            debut_nco = 0.5*(new_stop5[0]+new_stop3[0])
            fin_nco = 0.5*(new_start5[-1]+new_start3[-1])

#----------------------
# Chromatides du CO

# On finit par les chromatides du CO.
# On recode d'abord le pattern des chr in CO, soit chr_in_CO.

# Comme il s'agit d'un CO, on a forcément 1 chr rouge et 1 chr bleue
# (sinon le CO n'aurait pu être détecté).

            co_s_pat = 1000*["a"] # pseudo-liste
            for segm in range(dernier44,le44suivant+1):
                quadruplet = str(s_pat[segm][2*chr_in_CO[0]:2*chr_in_CO[0]+2]) + \
                             str(s_pat[segm][2*chr_in_CO[1]:2*chr_in_CO[1]+2])
                
                co_s_pat[segm] = quadruplet + "0011"

# On regroupe les segments identiques avec la nouvelle définition de pattern
# en incluant dernier44 = élément 0 de new_start, etc.
# et k = dernier 4:4 = élément -1 de new_start, etc.
            
            plage = dernier44
            new_start5 = [] ; new_start3 = []
            new_stop5 = [] ; new_stop3 = []
            new_pat = []
            for segm in range(plage+1,le44suivant+1):
                ref_start5 = s_start5[plage]
                ref_start3 = s_start3[plage]
                ref_pattern = co_s_pat[plage]

                if co_s_pat[segm] == ref_pattern and segm == le44suivant :
                    new_start5.append(ref_start5)
                    new_start3.append(ref_start3)
                    new_stop5.append(s_stop5[segm])
                    new_stop3.append(s_stop3[segm])
                    new_pat.append(ref_pattern)

                elif co_s_pat[segm] != ref_pattern :
                    new_start5.append(ref_start5)
                    new_start3.append(ref_start3)
                    new_stop5.append(s_stop5[segm-1])
                    new_stop3.append(s_stop3[segm-1])
                    new_pat.append(ref_pattern)
                    plage = segm
                    if segm == le44suivant :
                        new_start5.append(s_start5[segm])
                        new_start3.append(s_start3[segm])
                        new_stop5.append(s_stop5[segm])
                        new_stop3.append(s_stop3[segm])
                        new_pat.append(co_s_pat[segm])

# Recodage systématique de g_type pour les 5:3a, 6:2*, 4:4* et 4:4co

            g_type = ""
            list44 = [new_pat[0]]
            list44sh = [new_pat[0]]  # 4:4 sans hétéroduplex
            list53 = []
            list35 = []
            list62 = []
            list26 = []
            list71 = []
            list17 = []
            dic_list = {"5:3":list53,"3:5":list35,"6:2":list62,"2:6":list26,
                        "7:1":list71,"1:7":list17,"4:4":list44}
            letter = ["","a","b","c","d","e","f","g","h"]
            
            for segm in range(1,len(new_pat)):
            # Inclusion du dernier 4:4 = k, mais pas du premier = dernier44
            
                segm_count = str(new_pat[segm].count("1")) + ":" + \
                             str(new_pat[segm].count("0"))

            # Indication des formes alternatives dans tous les cas, même 4:4                   

                found = 0
                for model in range(len(dic_list[segm_count])):
                    if new_pat[segm] == dic_list[segm_count][model]:
                        segm_type = segm_count + letter[model]
                        found = 1
                        break
                if not found :
                    segm_type = segm_count + letter[len(dic_list[segm_count])]
                    dic_list[segm_count].append(new_pat[segm])

            # Indication des hDNA dans les formes paires
            
                if segm_count in ["4:4","6:2","2:6"]:
                    hDNA = 0
                    for chromatide in range(4):
                        if new_pat[segm][2*chromatide] != new_pat[segm][2*chromatide + 1]:
                            hDNA = 1
                            break
                    if hDNA :
                        segm_type += "i"

            # Indication des 4:4co par rapport au 4:4 précédent
            
                if segm_count == "4:4" and "i" not in segm_type \
                   and new_pat[segm] != list44sh[-1] :  # dernier élément de list44sh
                    segm_type += "CO"
                    list44sh.append(new_pat[segm]) # addition uniquement des 4:4 sans hDNA

            # Addition à g_type
            
                if segm == 1 : 
                    g_type += "(" + segm_type
                else :
                    g_type += ")_(" + segm_type
                

# Construction de g_len (longueur de tous les segments en excluant les 4:4 aux extrémités)

            g_type += ")"
            if g_type != "(4:4aCO)" :
                g_len = ""
                seg_len = 0.5*(new_stop3[1]+new_stop5[1]-\
                              (new_start3[1]+new_start5[1]))
                g_len += str(round(seg_len))
                
                for j in range(2,len(new_pat)-1):
                    seg_len = 0.5*(new_stop3[j]+new_stop5[j]-\
                           (new_start3[j]+new_start5[j]))
                    g_len += "_" + str(round(seg_len))
            else:
                g_len = 0

                            
# Détermination nb_co et longueur totale

            nb_co = g_type.count("CO")

            longueur = round(0.5*(new_start3[-1]+new_start5[-1]-\
                                 (new_stop3[0]+new_stop5[0]))) # -1 = dernier 4:4 de new_*

# Détermination de overlap

#            debut_nco = 0.5*(new_stop5[0]+new_stop3[0])
#            fin_nco = 0.5*(new_start5[-1]+new_start3[-1])

            debut_co = 0.5*(new_stop5[0]+new_stop3[0])
            fin_co = 0.5*(new_start5[-1]+new_start3[-1])   

            overlap = 0
            if debut_co <= fin_nco <= fin_co or \
               debut_co <= debut_nco <= fin_co or \
               (debut_nco < debut_co and fin_co < fin_nco):
                overlap = 1

            if overlap == 1 :
                distance = 0
            else :
                if fin_nco < debut_co :
                    distance = round(debut_co - fin_nco)
                else :
                    distance = round(debut_nco - fin_co)
            
# Réinitialisation des variables avec système pour arrêter le programme
# après un id particulier pour vérifications

            if nb_id != 1000 :
                nb_id += 1
                dernier44 = le44suivant
            else :
                arret = 1
                break

###Run old script for the U events###
for chromo in range(1,17):

    s_start5 = [] ; s_start3 = [] ; s_stop5 = [] ; s_stop3 = []
    s_pat = [] ; s_type = [] ; s_nb_snp = []  # s = segment

    for i in range(1,len(a)-1):
        if int(a[i].split("\t")[0]) != chromo :
            continue
        line = a[i].split("\t")
        s_start5.append(int(line[1]))
        s_start3.append(int(line[2]))
        s_stop5.append(int(line[3]))
        s_stop3.append(int(line[4]))
        pattern = line[5].replace("[","")
        pattern = pattern.replace("]","")
        pattern = pattern.split(", ")
        s_pat.append(pattern)
        s_type.append(line[6])
        s_nb_snp.append(int(line[7]))

    #----------------------------------------------------------------------

    # Gestion du premier event

    dernier44 = 0
    nb_snp = 0
    g_type = ""
        
    for k in range(1,len(s_type)):
            
        if s_type[k] == "4:4" and \
            0.5*(s_stop3[k]+s_stop5[k] - (s_start3[k]+s_start5[k])) > limite :

            if k > dernier44 + 1:
                g_type = s_type[dernier44+1]
                nb_snp += s_nb_snp[dernier44+1]
                
            for j in range(2,k):
                g_type = g_type + "_" + s_type[j]
                nb_snp += s_nb_snp[j]

                    
            chr_involved = []
            chr_involved_pat = []
            for chromatid in range(4):
                chromatid_pattern = []
                for segment in range(k+1):
                    string_pat = str(s_pat[segment][2*chromatid]) + \
                                 str(s_pat[segment][2*chromatid + 1])
                    chromatid_pattern.append(string_pat)
                if len(set(chromatid_pattern)) > 1:     # indent
                    chr_involved.append(chromatid)
                    chr_involved_pat.append(chromatid_pattern)

            if len(chr_involved) < 3 :
                    
                strand_involved = []
                for strand in range(8):
                    strand_pattern = []
                    for segment in range(k+1):
                        strand_pattern.append(s_pat[segment][strand])
                    if len(set(strand_pattern)) > 1:
                        strand_involved.append(strand)

            if len(chr_involved) == 2 :

                pattern_sister = []
                for chromatid in range(4):
                    pattern_sister.append(str(s_pat[dernier44][2*chromatid]) + \
                                    str(s_pat[dernier44][2*chromatid + 1]))
                if pattern_sister[chr_involved[0]] == \
                    pattern_sister[chr_involved[1]] :
                    sister = 1
                else :
                    sister = 0
                        
    #------------------------------------------------------
                        
            if s_type[dernier44] == "4:4":
                pat_before = s_pat[dernier44]
                pat_after = s_pat[k]
                if pat_after == pat_before:
                    co = "NCO"
                else :
                    co = "CO"
            else :
                co = "U"
                        
    #------------------------------------------------------
                        
            if len(chr_involved) in [3,4]:
                classe = len(chr_involved)
                    
            elif len(chr_involved) == 1:
                if len(strand_involved) == 1:
                    classe = "1_1brin"
                else:
                    classe = "1_2brins"
                        
            elif len(chr_involved) == 2:
                if sister == 1:
                    classe = "2_sis"
                else :
                    classe = "2_nonsis"
                        
    #------------------------------------------------------
            len_mid = round(0.5*(s_stop3[k-1]+s_stop5[k-1]-(s_start3[1]+s_start5[1])))
            len_min = s_stop5[k-1]-s_start3[1] + 1
            len_max = s_stop3[k-1]-s_start5[1] - 1
            if nb_snp == 0:
                len_min = 0
                len_mid = round(0.5*len_max)
                    
                
            out3.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\
                        graph_bin_annot3(%s,%s,%s)\n" %(
                nb_id,chromo,
                s_start5[1],s_start3[1],
                s_stop5[k-1],s_stop3[k-1],g_type,co,
                len(chr_involved),classe,len_min,len_mid,len_max,nb_snp,
                chromo,((s_start5[1]-1000)//1000)*1000,
                ((s_stop3[k-1]+1000)//1000)*1000))

            nb_id += 1
            dernier44 = k
            nb_snp = 0
            g_type = ""
            break

    #----------------------------------------------------------------
            
    # Gestion des events suivants

    for k in range(dernier44+1,len(s_pat)):            
                    
        if (s_type[k] == "4:4" and 0.5*(s_stop3[k]+s_stop5[k] - (s_start3[k]+s_start5[k])) > limite) or \
            k == len(s_type) - 1:

            if k == len(s_type) - 1 and s_type[k] != "4:4" :
                fin = k
            else :
                fin = k-1

            if k > dernier44 + 1:
                g_type = s_type[dernier44+1]
                nb_snp += s_nb_snp[dernier44+1]
                
            for j in range(dernier44+2,fin+1):
                g_type = g_type + "_" + s_type[j]

            for j in range(dernier44+2,k):
                nb_snp += s_nb_snp[j]
                    
            chr_involved = []
            for chromatid in range(4):
                chromatid_pattern = []
                for segment in range(dernier44,k+1):
                    string_pat = str(s_pat[segment][2*chromatid]) + \
                                    str(s_pat[segment][2*chromatid + 1])
                    chromatid_pattern.append(string_pat)
                if len(set(chromatid_pattern)) > 1:
                    chr_involved.append(chromatid)

            if len(chr_involved) < 3 :
                    
                strand_involved = []
                for strand in range(8):
                    strand_pattern = []
                    for segment in range(dernier44,k+1):
                        strand_pattern.append(s_pat[segment][strand])
                    if len(set(strand_pattern)) > 1:
                        strand_involved.append(strand)

            if len(chr_involved) == 2 :

                pattern_sister = []
                for chromatid in range(4):
                    pattern_sister.append(str(s_pat[dernier44][2*chromatid]) + \
                                    str(s_pat[dernier44][2*chromatid + 1]))
                if pattern_sister[chr_involved[0]] == \
                    pattern_sister[chr_involved[1]] :
                    sister = 1
                else :
                    sister = 0

    #-------------------------------------------------

            if s_type[k] == "4:4" :
                pat_before = s_pat[dernier44]
                pat_after = s_pat[k]
                if pat_after == pat_before:
                    co = "NCO"
                else :
                    co = "CO"
            else :
                co = "U"

    #------------------------------------------------------
                        
            if len(chr_involved) in [3,4]:
                classe = len(chr_involved)
                    
            elif len(chr_involved) == 1:
                if len(strand_involved) == 1:
                    classe = "1_1brin"
                else:
                    classe = "1_2brins"
                        
            elif len(chr_involved) == 2:
                if sister == 1:
                    classe = "2_sis"
                else :
                    classe = "2_nonsis"
                        
    #------------------------------------------------------

            len_mid = round(0.5*(s_stop3[fin]+s_stop5[fin]-\
                            (s_start3[dernier44+1]+s_start5[dernier44+1])))
            len_min = s_stop5[fin]-s_start3[dernier44+1] + 1
            len_max = s_stop3[fin]-s_start5[dernier44+1] - 1
            if nb_snp == 0:
                len_min = 0
                len_mid = round(0.5*len_max)
                
            out3.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\
                        graph_bin_annot3(%s,%s,%s)\n" %(
                nb_id,chromo,
                s_start5[dernier44+1],s_start3[dernier44+1],
                s_stop5[fin],s_stop3[fin],g_type,co,
                len(chr_involved),classe,len_min,len_mid,len_max,nb_snp,
                chromo,((s_start5[dernier44+1]-1000)//1000)*1000,
                ((s_stop3[fin]+1000)//1000)*1000))

            nb_id += 1
            dernier44 = k
            nb_snp = 0
            g_type = ""

       
out.close()
out3.close()


