#!/bin/bash
#Script to do joint-variant calling for PV AmpliSeq data for targeted regions for PV Vietnam design
#version 28/11/2022 E Kattenberg for github

#set file paths:
#1 to the folder with the demultiplexed fastq-files from the AMpliSeq sequencing reaction on MiSeq
Thepathfq=/home/user/fq/
#2 To the folder where the output will go
Thepath=/home/user/project_name/

#3. To the sample map with all samples and link to vcf to include in the joint genotyping
samplelist=$Thepath/Sample_list_PV_VTN_experiment1.txt
samplemap=$Thepath/alignments/database/samplemap_PV_experiment1_sort.txt
#see notes in previous script (2) how to make this samplemap

#set paths to programs on my computer
gatk=~/bin/gatk-4.1.2.0
TMPDIR=/home/user/tmp

#set path to reference genome  
#this is not the latest version, but the one I used for the Pv AmpliSeq Vietnam data publised in Kattenberg et al. 2022 
#You can update to latest version, download from PlasmoDB: https://plasmodb.org/plasmo/app/downloads/Current_Release/PvivaxP01/
ref=/home/user/PvP01/PlasmoDB-46/PlasmoDB-46_PvivaxP01_Genome.fasta


###########################################################################
#7. make database to combine genotypes specifying the target regions only #
###########################################################################

#Go to the folder with the samplemap
cd $Thepath/alignments/database

#create the database with all samples
#change the name of the database for your own study.
#change the intervals if you use a different AmpliSeq design, this one is for PV Vietnam

$gatk/gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport \
--genomicsdb-workspace-path PV_VTN_allsamples_Nov2022_database \
--sample-name-map $samplemap \
--tmp-dir=$TMPDIR \
-L PvP01_01_v1:121120-121446 \
-L PvP01_01_v1:164325-164650 \
-L PvP01_01_v1:442168-445683 \
-L PvP01_01_v1:671599-671928 \
-L PvP01_02_v1:153713-158934 \
-L PvP01_02_v1:197944-198259 \
-L PvP01_02_v1:373468-373796 \
-L PvP01_02_v1:594651-594975 \
-L PvP01_03_v1:127832-128135 \
-L PvP01_03_v1:220801-221130 \
-L PvP01_03_v1:334581-334905 \
-L PvP01_03_v1:552739-553924 \
-L PvP01_03_v1:782109-782432 \
-L PvP01_04_v1:420827-421155 \
-L PvP01_04_v1:514823-515147 \
-L PvP01_04_v1:885372-885697 \
-L PvP01_05_v1:192353-192680 \
-L PvP01_05_v1:451975-452280 \
-L PvP01_05_v1:693386-693708 \
-L PvP01_05_v1:701363-701690 \
-L PvP01_05_v1:1077249-1079282 \
-L PvP01_05_v1:1079319-1079539 \
-L PvP01_05_v1:1285062-1285380 \
-L PvP01_06_v1:45733-46061 \
-L PvP01_06_v1:278016-278338 \
-L PvP01_06_v1:646077-646399 \
-L PvP01_06_v1:944634-944957 \
-L PvP01_07_v1:260890-261212 \
-L PvP01_07_v1:595121-595445 \
-L PvP01_07_v1:754446-754775 \
-L PvP01_07_v1:1020359-1020683 \
-L PvP01_07_v1:1210982-1211301 \
-L PvP01_08_v1:44930-45091 \
-L PvP01_08_v1:53195-53514 \
-L PvP01_08_v1:442285-442612 \
-L PvP01_08_v1:878901-879225 \
-L PvP01_08_v1:1420967-1421291 \
-L PvP01_08_v1:1590968-1591291 \
-L PvP01_09_v1:303127-303447 \
-L PvP01_09_v1:448658-448978 \
-L PvP01_09_v1:539199-539524 \
-L PvP01_09_v1:1070335-1070662 \
-L PvP01_09_v1:1101045-1101353 \
-L PvP01_09_v1:1459281-1461078 \
-L PvP01_09_v1:1838486-1838812 \
-L PvP01_09_v1:1883812-1884132 \
-L PvP01_10_v1:351609-351936 \
-L PvP01_10_v1:478690-483160 \
-L PvP01_10_v1:483513-483836 \
-L PvP01_10_v1:523552-523876 \
-L PvP01_10_v1:826080-831521 \
-L PvP01_10_v1:1130890-1131215 \
-L PvP01_10_v1:1385634-1385960 \
-L PvP01_11_v1:159891-162171 \
-L PvP01_11_v1:245545-245867 \
-L PvP01_11_v1:383879-384206 \
-L PvP01_11_v1:503806-504133 \
-L PvP01_11_v1:720133-720439 \
-L PvP01_11_v1:1137154-1137482 \
-L PvP01_11_v1:1145107-1145430 \
-L PvP01_11_v1:1226851-1227143 \
-L PvP01_11_v1:1448738-1449062 \
-L PvP01_11_v1:1867852-1868180 \
-L PvP01_12_v1:94251-94573 \
-L PvP01_12_v1:484766-487079 \
-L PvP01_12_v1:844106-844430 \
-L PvP01_12_v1:1116221-1116544 \
-L PvP01_12_v1:1400006-1400335 \
-L PvP01_12_v1:1462978-1463300 \
-L PvP01_12_v1:1844773-1845096 \
-L PvP01_12_v1:1988181-1988504 \
-L PvP01_12_v1:2441294-2446380 \
-L PvP01_13_v1:65981-66304 \
-L PvP01_13_v1:162629-162913 \
-L PvP01_13_v1:659361-659691 \
-L PvP01_13_v1:1769965-1770289 \
-L PvP01_14_v1:344814-345140 \
-L PvP01_14_v1:743252-743582 \
-L PvP01_14_v1:1160566-1160887 \
-L PvP01_14_v1:1229330-1229655 \
-L PvP01_14_v1:1258494-1258818 \
-L PvP01_14_v1:1266125-1266440 \
-L PvP01_14_v1:1269755-1272327 \
-L PvP01_14_v1:1910971-1911297 \
-L PvP01_14_v1:2051946-2058196 \
-L PvP01_14_v1:3004151-3004474 



#############################################################
#8. Joint genotyping										#
#Conduct joint genotyping using all samples in the database #
#############################################################

#change database and output vcf name for your study

$gatk/gatk --java-options "-Xmx4g" GenotypeGVCFs \
-R $ref \
-V gendb://PV_VTN_allsamples_Nov2022_database \
-O VTN_Ampliseq_PV_Nov2022.vcf.gz \
--tmp-dir=$TMPDIR


#the final vcf file will be in the location: 
#$Thepath/alignments/database/VTN_Ampliseq_PV_Nov2022.vcf.gz
#There is an additional script to filter the variants retaining only high quality variants for analysis. 

##END##