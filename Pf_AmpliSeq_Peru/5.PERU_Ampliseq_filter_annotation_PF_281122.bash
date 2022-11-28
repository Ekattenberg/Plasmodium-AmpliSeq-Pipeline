#!/bin/bash
# Script to filter variants and do annotation for data from PF Peru Ampliseq
#version 28/11/2022 E Kattenberg for github

#set file paths:
#1 to the folder with the demultiplexed fastq-files from the AMpliSeq sequencing reaction on MiSeq
Thepathfq=/home/user/fq/
#2 To the folder where the output will go
Thepath=/home/user/project_name/

#3. To the list of samples to include in the variantcalling
samplelist=$Thepath/Sample_list_PF_Peru_experiment1.txt

#set paths to programs
sam=~/bin/samtools-1.9
gatk=~/bin/gatk-4.1.2.0
TMPDIR=/home/user/tmp

##set path to reference genome  
#this is not the latest version, but the one I used for the Pf AmpliSeq Peru data publised in Kattenberg et al. 2022 
#You can update to latest version, download from PlasmoDB: https://plasmodb.org/plasmo/app/downloads/Current_Release/Pfalciparum3D7/
ref=/home/user/3D7_genome/PlasmoDB-44_Pfalciparum3D7_Genome.fasta

######################################
#9. Filter vcf  					 #
######################################

#Go to the folder with the joint vcf file
cd $Thepath/alignments/database

#Subset to make sure the vcf includes only variants in the target region
$gatk/gatk SelectVariants \
-R $ref \
-V PERU_Ampliseq_PF_Nov2022.vcf.gz \
-O PERU_Ampliseq_PF_Nov2022_panel.vcf.gz \
-L Pf3D7_01_v3:179679-179985 \
-L Pf3D7_01_v3:180375-180467 \
-L Pf3D7_01_v3:190106-191218 \
-L Pf3D7_01_v3:191661-192230 \
-L Pf3D7_01_v3:192264-192925 \
-L Pf3D7_01_v3:192961-194344 \
-L Pf3D7_01_v3:194371-195570 \
-L Pf3D7_01_v3:195762-196800 \
-L Pf3D7_01_v3:196836-200312 \
-L Pf3D7_01_v3:200448-200978 \
-L Pf3D7_01_v3:204949-205193 \
-L Pf3D7_01_v3:339163-339466 \
-L Pf3D7_01_v3:464638-466535 \
-L Pf3D7_01_v3:466803-470193 \
-L Pf3D7_02_v3:519270-519567 \
-L Pf3D7_02_v3:694157-694456 \
-L Pf3D7_03_v3:361017-361320 \
-L Pf3D7_03_v3:849331-849578 \
-L Pf3D7_04_v3:532089-532395 \
-L Pf3D7_04_v3:691941-692036 \
-L Pf3D7_04_v3:747935-749935 \
-L Pf3D7_04_v3:770125-770439 \
-L Pf3D7_05_v3:921740-921999 \
-L Pf3D7_05_v3:957777-958361 \
-L Pf3D7_05_v3:958396-959096 \
-L Pf3D7_05_v3:959153-959988 \
-L Pf3D7_05_v3:960166-961735 \
-L Pf3D7_05_v3:962020-962193 \
-L Pf3D7_05_v3:1188236-1188508 \
-L Pf3D7_05_v3:1214304-1214615 \
-L Pf3D7_06_v3:148637-148858 \
-L Pf3D7_06_v3:635908-636217 \
-L Pf3D7_07_v3:403083-404022 \
-L Pf3D7_07_v3:404296-406440 \
-L Pf3D7_07_v3:455433-455668 \
-L Pf3D7_07_v3:782049-782219 \
-L Pf3D7_08_v3:500969-501114 \
-L Pf3D7_08_v3:548041-549952 \
-L Pf3D7_08_v3:549982-550332 \
-L Pf3D7_08_v3:803010-803329 \
-L Pf3D7_08_v3:1373986-1374244 \
-L Pf3D7_08_v3:1374249-1374474 \
-L Pf3D7_08_v3:1374486-1374705 \
-L Pf3D7_08_v3:1374711-1375419 \
-L Pf3D7_09_v3:230922-231210 \
-L Pf3D7_09_v3:1005112-1005419 \
-L Pf3D7_10_v3:340949-341230 \
-L Pf3D7_10_v3:1172615-1172823 \
-L Pf3D7_11_v3:416182-416488 \
-L Pf3D7_11_v3:874927-875071 \
-L Pf3D7_11_v3:1294408-1295115 \
-L Pf3D7_11_v3:1295124-1295390 \
-L Pf3D7_11_v3:1295423-1295672 \
-L Pf3D7_11_v3:1505334-1505635 \
-L Pf3D7_12_v3:717790-718460 \
-L Pf3D7_12_v3:718476-719820 \
-L Pf3D7_12_v3:1126866-1127122 \
-L Pf3D7_12_v3:1552012-1552292 \
-L Pf3D7_12_v3:1611167-1611448 \
-L Pf3D7_12_v3:2091976-2094267 \
-L Pf3D7_13_v3:1595960-1596123 \
-L Pf3D7_13_v3:1724599-1727164 \
-L Pf3D7_13_v3:1827453-1827763 \
-L Pf3D7_13_v3:2503123-2503359 \
-L Pf3D7_13_v3:2503390-2505580 \
-L Pf3D7_13_v3:2840639-2841793 \
-L Pf3D7_14_v3:293359-293982 \
-L Pf3D7_14_v3:294267-294832 \
-L Pf3D7_14_v3:832462-832742 \
-L Pf3D7_14_v3:1381895-1381990 \
-L Pf3D7_API_v3:24048-34146 \
-L Pf_M76611:3384-4634


#Hard Filtering variants
#MQ >50
$gatk/gatk VariantFiltration \
   -R $ref \
   -V PERU_Ampliseq_PF_Nov2022_panel.vcf.gz \
   -O PERU_Ampliseq_PF_Nov2022_panel_MQ.vcf.gz \
   --filter-name "MQ50" \
   --filter-expression "MQ0 < 50" 

#QUAl > 30
$gatk/gatk VariantFiltration \
   -R $ref \
   -V PERU_Ampliseq_PF_Nov2022_panel_MQ.vcf.gz \
   -O PERU_Ampliseq_PF_Nov2022_panel_MQQUAL.vcf.gz \
   --filter-name "QUAL" \
   --filter-expression "QUAL < 30.00"

#Overall DP > 500
$gatk/gatk VariantFiltration \
   -R $ref \
   -V PERU_Ampliseq_PF_Nov2022_panel_MQQUAL.vcf.gz \
   -O PERU_Ampliseq_PF_Nov2022_panel_MQQUALDP.vcf.gz \
   --filter-name "DP500" \
   --filter-expression "DP < 500" 

## QD SOR ReadPosRankSum
$gatk/gatk VariantFiltration \
    -V PERU_Ampliseq_PF_Nov2022_panel_MQQUALDP.vcf.gz \
    -filter "QD < 1.0" --filter-name "QD1" \
    -filter "SOR > 4.0" --filter-name "SOR3" \
    -filter "ReadPosRankSum < -10.0" --filter-name "ReadPosRankSum-10" \
    -O PERU_Ampliseq_PF_Nov2022_panel_MQQUALDPQDSORRPRS.vcf.gz

#remove the variants that are filtered
$gatk/gatk SelectVariants \
-R $ref \
-V PERU_Ampliseq_PF_Nov2022_panel_MQQUALDPQDSORRPRS.vcf.gz \
--exclude-filtered \
-O PERU_Ampliseq_PF_Nov2022_panel_filteredPASS.vcf.gz

#### Filter on GT DP of 5 (on top of Qual 30 and DP 100)
$gatk/gatk VariantFiltration \
-R $ref \
-V PERU_Ampliseq_PF_Nov2022_panel_filteredPASS.vcf.gz \
-O PERU_Ampliseq_PF_Nov2022_panel_filteredPASS_DP5.vcf.gz \
-G-filter-name "GTDP5" \
-G-filter "DP < 6" \
--set-filtered-genotype-to-no-call

#############################
#10. Annotation of variants	#
#############################

#annotate the variants using the program snpEff and dataset created with PlasmoDB 3D7 v44 (see manuals snpEff how to create a dataset for annotation with snpEff)
java -Xmx4g -jar ~/snpEff/snpEff.jar -c ~/snpEff/snpEff.config -v Pf3D7.44 PERU_Ampliseq_PF_Nov2022_panel_filteredPASS_DP5.vcf.gz > PERU_Ampliseq_PF_Nov2022_panel_filteredPASS_DP5_ANN.vcf

# Index the vcf and then compress the file
$gatk/gatk IndexFeatureFile \
-F PERU_Ampliseq_PF_Nov2022_panel_filteredPASS_DP5_ANN.vcf

bgzip -c PERU_Ampliseq_PF_Nov2022_panel_filteredPASS_DP5_ANN.vcf > PERU_Ampliseq_PF_Nov2022_panel_filteredPASS_DP5_ANN.vcf.gz
tabix -p vcf PERU_Ampliseq_PF_Nov2022_panel_filteredPASS_DP5_ANN.vcf.gz


##END##
