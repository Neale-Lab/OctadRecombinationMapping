# OctadRecombinationMapping
A collection of R and Python scripts used to map recombination events in yeast octads and tetrads.

This repository accompanies the publication "Separable roles of the DNA damage response kinase Mec1(ATR) and its activator Rad24(RAD17) within the regulation of meiotic recombination":
https://www.researchgate.net/publication/329741247_Separable_roles_of_the_DNA_damage_response_kinase_Mec1ATR_and_its_activator_Rad24RAD17_within_the_regulation_of_meiotic_recombination

Program Requirements
•	R for Windows 3.5.1: https://www.r-project.org/
•	R Packages: plyr, e1071
•	Python 3.7.0: https://www.python.org/
•	Bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
•	SAMtools: http://samtools.sourceforge.net/
•	PySAMStats: https://github.com/alimanfoo/pysamstats

For all R scripts, there are additional comments written inside the script which explain how to use them.

Datatracks from external publications are used to plot alongside recombination events. The datatracks are included in this repository because they have been modified in order to work with the scripts presented here. The tracks are the Spo11 profile and hotspot strength values (Mohibullah and Keeney, 2017; Pan et al., 2011), Rec114/Mer2/Mei1 (RMM) profile (Panizza et al., 2011), Rec8 peaks (Ito et al., 2014), and transcribed regions (Nagalakshmi et al., 2008).

Also included in this repository are some example files for key stages in the pipeline:
AEvents_OM1_1500.txt
Annot_OM1
Events_OM1_1500
NAnnot_OM1
NEvents_OM1_1500
OM1_Annotated_Images.pdf
OM1_binary.txt
OM1_gr2.txt
OM1_Images.pdf
Orientation_table.txt
SEvents_OM1_1500
UEvents_OM1_1500

Instructions

Step 1. Alignment.
Input: Paired FASTQ files, appropriate reference genome files. 
Example FastQ files are publicly available at the NCBI Sequence Read Archive (accession numbers SRP151982, SRP152540, SRP152953).
Example reference files are provided. There are two versions of the SK1 genome, one for NDT80-AR samples and one for any other genotypes.
Each sample is aligned to an SK1 and an S288c genome. Files are aligned using bowtie2.
Example command:
bowtie2 -p 8 -X 1000 --local --dovetail -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 --mp 5,1 --rg-id 1 --rg PU:1 --rg LB:1  --rg SM:1 --rg PL:1 -x S288c_2u_mt --trim3 0 --trim5 0 -1 OM1A_R1.fastq -2 OM1A_R2.fastq -S OM1A_S288c.SAM
Output: Aligned SAM files for S288c and SK1 genomes.

Step 2. Sort SAM file.
Input: Aligned SAM files for S288c and SK1 genomes.
To sort the SAM file, use SAMtools “view” and “index” commands. 
Example commands:
samtools sort OM1A_S288c.SAM -o OM1A_S288c_sorted.BAM
samtools index OM1A_S288c_sorted.BAM
Output: Sorted BAM file and BAM index file for S288c and SK1 genomes.

Step 3. Find number of reads per base.
Input: Sorted BAM files, BAM index files, appropriate reference genome files.
To find the number of reads that matched each base or indel at each position, run the PySAMstats “Variation” module, against the SK1 and S288c genomes.
Example command:
pysamstats -f S288c_2u_mt.fsa --type variation OM1A_sorted.bam > OM1A.txt 
Use the same reference genome files as the original alignment.
Output: .txt file containing the number of reads of each base at each position, for S288c and SK1 genomes. For example, “OM1A_S288c”.

Step 4. Isolate variant positions.
Input: .txt file from PySAMStats containing variant read numbers
To isolate variants, using RStudio, run the program ‘Run_isolator_octads.R’ or ‘Run_isolator_tetrads.R’ depending on whether you made octads or tetrads. Change filenames in the script to match your samples. The script depends on functions SNP_indel_isolator_S288c.R, SNP_indel_isolator_SK1.R, and a variant table. 
Variant tables are provided. ‘Variant_Table_NDT80_6’ is for NDT80-AR samples, and ‘Variant_Table_6’ can be used for any other sample.
Output: .txt file containing the number of reads of each base at each variant position, for SK1 and S288c genomes. For example, “OM1AOut_S288c.txt”.

Step 5. Combine members of an octad/tetrad and call the variants. 
Input: .txt file containing the number of reads of each base at each variant position
To call variants, run the R programs ‘Run_Caller_S288c’ and ‘Run_Caller_SK1’ on each sample. You will need to alter the filename inside the script. The scripts depend on the function scripts SNP_Indel_Caller_S288c and SNP_Indel_Caller_SK1.
Output: .txt files containing variant calls against the SK1 and S288c genomes. For example, “OM1_combined_S288c.txt”.

Step 6. Reconcile variant calls.
Input: .txt files containing variant calls against the SK1 and S288c genomes.
To reconcile the SK1 and S288c calls, run the script ‘Reconcile_SK1_S288c_calls.R’. You will need to change filenames inside the script. Also note, different files need to be input for msh2Δ tetrads with heteroduplex calls which you will need to alter if you are using this sample type, this process is explained inside the script.
Output: a “gr2” file that is compatible with the event viewer in R, and a “binary” file containing the final variant calls, compatible with the event caller. For example, “OM1_gr2”, “OM1_binary”.
Example binary files are in the supplementary files for this publication.

Step 7. Collect variant patterns on each chromosome.
Input: Binary files.
To collect variant changes into patterns, use the Python script ‘AnnotationScript’. It is necessary to edit the script to contain the correct filenames. Then run the script for each binary file. This script requires the file “ChrLen2011_Rep3” which contains the chromosome lengths.
The annot files output needs to be modified slightly in order to be compatible with the next step. This uses an R script called ‘annot files modifier.R’. Run this for each sample.
Output: annotated chromosome files e.g “NAnnot_OM1”

Step 8. Call events.
Input: Annotated chromosome files.
To call recombination events, run the Python script ‘EventsScript.py’. Edit the file to contain correct filenames. You can also adjust the minimum distance between events, e.g. if the limit is 1500, events occurring closer than 1500bp will be considered part of the same event.
This produces two files, “Events” and “UEvents”. The UEvents file is needed for events at the end of the chromosome. To combine the two files, use the R script ‘Run_event_sorterV3.R’. This also sorts the events in order of complexity which is useful for annotation.
Output: Sorted Events file, for example “SEvents_OM1_1500”

Step 9. Combine and image events.
Input: Sorted Events files, gr2 files, datatracks
To make combined master tables containing all the octad data, run the R script ‘Combining Events TablesV2.R’ to create master tables. These are input into the imager.
Create images of every event using the R script ‘Event_Imager_Unannotatedv11.R’. This requires a number of datatracks which are in provided. 
Output: Master table of events, images of each event.

Step 10. Manually annotate events.
Input: Event images and sorted event files
Manually examine each image to check the category is correct, look for errors and write notes. Edit the sorted events file with any changes, and rename it “AEvents” for annotated events. For example, ‘AEvents_OM1_1500.txt’.
Guidance for event categories can be found in this publication or in Marsolier Kergoat et al., 2017.
Also record the orientation of the DNA strands if possible, in the Orientation_table file.
Output: Annotated events tables, Orientation table.

Step 11. Make annotated images.
Input: Annotated events files, orientation table, master gr2 file, datatracks
To create a post-annotation master table, use the R script ‘Combining_Annotated_Events’. You don’t need to re-make the gr2 table. 
Create post-annotation images of every event using the R script “Event_Imager_Annotated”. This is similar to the first imager but includes some more information, such as strand orientation if you determined this.
Output: Annotated event table for analysis, annotated event images

Datatrack references

Mohibullah N, Keeney S (2017) Numerical and spatial patterning of yeast meiotic DNA breaks by Tel1. Genome Res 27: 278-288.

Panizza S, Mendoza MA, Berlinger M, Huang L, Nicolas A, Shirahige K et al. (2011) Spo11-accessory proteins link double-strand break sites to the chromosome axis in early meiotic recombination. Cell 146: 372-383.

Ito M, Kugou K, Fawcett JA, Mura S, Ikeda S, Innan H et al. (2014) Meiotic recombination cold spots in chromosomal cohesion sites. Genes Cells 19: 359-373.

Nagalakshmi U, Wang Z, Waern K, Shou C, Raha D, Gerstein M et al. (2008) The transcriptional landscape of the yeast genome defined by RNA sequencing. Science 320: 1344-1349.

Pan J, Sasaki M, Kniewel R, Murakami H, Blitzblau HG, Tischfield SE et al. (2011) A Hierarchical Combination of Factors Shapes the Genome-wide Topography of Yeast Meiotic Recombination Initiation. Cell 144: 719-731.


