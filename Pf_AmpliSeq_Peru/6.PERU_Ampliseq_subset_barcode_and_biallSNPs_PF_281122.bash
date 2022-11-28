#!/bin/bash
# Script to subset barcode variants and biallelic SNPs Peru Ampliseq
#version 28/11/2022 E Kattenberg for github

#set file paths:
#To the folder where the output will go
Thepath=/home/user/project_name/

#set paths to programs
sam=~/bin/samtools-1.9
gatk=~/bin/gatk-4.1.2.0
TMPDIR=/home/malariology/FA4_Ampliseq/tmp

##set path to reference genome  
#this is not the latest version, but the one I used for the Pf AmpliSeq Peru data publised in Kattenberg et al. 2022 
#You can update to latest version, download from PlasmoDB: https://plasmodb.org/plasmo/app/downloads/Current_Release/Pfalciparum3D7/
ref=/home/user/3D7_genome/PlasmoDB-44_Pfalciparum3D7_Genome.fasta

#set path to folder with vcf
cd $Thepath/alignments/database


####extract GT table
$gatk/gatk VariantsToTable -V PERU_Ampliseq_PF_Nov2022_panel_filteredPASS_DP5_ANN.vcf.gz -F CHROM -F POS -F TYPE -GF GT -O PERU_Ampliseq_PF_Nov2022_panel_filteredPASS_DP5_ANN_GT.table


####Filter to SNP barcode regions only

$gatk/gatk SelectVariants \
-R $ref \
-V PERU_Ampliseq_PF_Nov2022_panel_filteredPASS_DP5_ANN.vcf.gz \
-O PERU_Ampliseq_PF_Nov2022_panel_filteredPASS_DP5_ANN_barcode.vcf.gz \
-L Pf3D7_01_v3:205066 \
-L Pf3D7_01_v3:339436 \
-L Pf3D7_02_v3:519457 \
-L Pf3D7_02_v3:694307 \
-L Pf3D7_03_v3:361199 \
-L Pf3D7_03_v3:849476 \
-L Pf3D7_04_v3:691961 \
-L Pf3D7_04_v3:770292 \
-L Pf3D7_05_v3:921893 \
-L Pf3D7_05_v3:1188394 \
-L Pf3D7_06_v3:148827 \
-L Pf3D7_06_v3:636044 \
-L Pf3D7_07_v3:455494 \
-L Pf3D7_07_v3:782111 \
-L Pf3D7_08_v3:501042 \
-L Pf3D7_08_v3:803172 \
-L Pf3D7_09_v3:231065 \
-L Pf3D7_09_v3:1005351 \
-L Pf3D7_10_v3:341106 \
-L Pf3D7_10_v3:1172712 \
-L Pf3D7_11_v3:874948 \
-L Pf3D7_11_v3:1505533 \
-L Pf3D7_12_v3:1127000 \
-L Pf3D7_12_v3:1552084 \
-L Pf3D7_13_v3:1595988 \
-L Pf3D7_13_v3:1827569 \
-L Pf3D7_14_v3:832594 \
-L Pf3D7_14_v3:1381943 

$gatk/gatk IndexFeatureFile \
-F PERU_Ampliseq_PF_Nov2022_panel_filteredPASS_DP5_ANN_barcode.vcf.gz

#extract GT per allele
$gatk/gatk VariantsToTable -V PERU_Ampliseq_PF_Nov2022_panel_filteredPASS_DP5_ANN_barcode.vcf.gz -F CHROM -F POS -F TYPE -GF GT -O PERU_Ampliseq_PF_Nov2022_panel_filteredPASS_DP5_ANN_barcode_GT.table

####Filter to biallic SNPs only and index file

$gatk/gatk SelectVariants \
-R $ref \
-V PERU_Ampliseq_PF_Nov2022_panel_filteredPASS_DP5_ANN.vcf.gz \
--select-type-to-include SNP \
--restrict-alleles-to BIALLELIC \
-O PERU_Ampliseq_PF_Nov2022_panel_filteredPASS_DP5_ANN_biSNPs.vcf.gz 

$gatk/gatk IndexFeatureFile \
-F PERU_Ampliseq_PF_Nov2022_panel_filteredPASS_DP5_ANN_biSNPs.vcf.gz

##END##
