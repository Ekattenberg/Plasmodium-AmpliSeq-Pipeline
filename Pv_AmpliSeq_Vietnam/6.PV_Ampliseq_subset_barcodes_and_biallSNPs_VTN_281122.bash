#!/bin/bash
# Script to subset barcode variants and biallelic SNPs Vietnam PV Ampliseq
#version 28/11/2022 E Kattenberg for github

#set file paths:
#To the folder where the output will go
Thepath=/home/user/project_name/

#set paths to programs
sam=~/bin/samtools-1.9
gatk=~/bin/gatk-4.1.2.0
TMPDIR=/home/malariology/FA4_Ampliseq/tmp

#set path to reference genome  
#this is not the latest version, but the one I used for the Pv AmpliSeq Vietnam data publised in Kattenberg et al. 2022 
#You can update to latest version, download from PlasmoDB: https://plasmodb.org/plasmo/app/downloads/Current_Release/PvivaxP01/
ref=/home/user/PvP01/PlasmoDB-46/PlasmoDB-46_PvivaxP01_Genome.fasta

#set path to folder with vcf
cd $Thepath/alignments/database


####extract GT table
$gatk/gatk VariantsToTable -V VTN_Ampliseq_PV_Nov2022_panel_filteredPASS_DP5_ANN.vcf.gz -F CHROM -F POS -F TYPE -GF GT -O VTN_Ampliseq_PV_Nov2022_panel_filteredPASS_DP5_ANN_GT.table

########################################
####Filter to SNP barcode regions only #
########################################
## for VivaxGEO barcode

$gatk/gatk SelectVariants \
-R $ref \
-V VTN_Ampliseq_PV_Nov2022_panel_filteredPASS_DP5_ANN.vcf.gz \
-O VTN_Ampliseq_PV_Nov2022_panel_filteredPASS_DP5_ANN_vivaxGEO_SNPS.vcf.gz \
-L PvP01_03_v1:220999 \
-L PvP01_05_v1:693449 \
-L PvP01_05_v1:1079461 \
-L PvP01_05_v1:1285366 \
-L PvP01_06_v1:646186 \
-L PvP01_07_v1:261070 \
-L PvP01_07_v1:595236 \
-L PvP01_08_v1:442363 \
-L PvP01_08_v1:879062 \
-L PvP01_08_v1:1420994 \
-L PvP01_09_v1:303291 \
-L PvP01_09_v1:448825 \
-L PvP01_09_v1:1101235 \
-L PvP01_09_v1:1884013 \
-L PvP01_10_v1:480601 \
-L PvP01_10_v1:483568 \
-L PvP01_10_v1:1131031 \
-L PvP01_11_v1:384043 \
-L PvP01_11_v1:503992 \
-L PvP01_11_v1:1137410 \
-L PvP01_11_v1:1226906 \
-L PvP01_11_v1:1448769 \
-L PvP01_11_v1:1868012 \
-L PvP01_12_v1:1116377 \
-L PvP01_12_v1:1400307 \
-L PvP01_12_v1:1463080 \
-L PvP01_12_v1:1988355 \
-L PvP01_13_v1:66121 \
-L PvP01_14_v1:344881 \
-L PvP01_14_v1:1160762 \
-L PvP01_14_v1:1229487 \
-L PvP01_14_v1:1258544 \
-L PvP01_14_v1:1266326 

$gatk/gatk IndexFeatureFile \
-F VTN_Ampliseq_PV_Nov2022_panel_filteredPASS_DP5_ANN_vivaxGEO_SNPS.vcf.gz

#extract GT per allele
$gatk/gatk VariantsToTable -V VTN_Ampliseq_PV_Nov2022_panel_filteredPASS_DP5_ANN_vivaxGEO_SNPS.vcf.gz -F CHROM -F POS -F TYPE -GF GT -O VTN_Ampliseq_PV_Nov2022_panel_filteredPASS_DP5_ANN_vivaxGEO_SNPS_GT.table

################################
## for Vietnam 42-SNP barcode
$gatk/gatk SelectVariants \
-R $ref \
-V VTN_Ampliseq_PV_Nov2022_panel_filteredPASS_DP5_ANN.vcf.gz \
--select-type-to-include SNP \
-O VTN_Ampliseq_PV_Nov2022_panel_filteredPASS_DP5_ANN_VTN_SNPs.vcf.gz \
-L PvP01_01_v1:121166 \
-L PvP01_01_v1:164620 \
-L PvP01_01_v1:671781 \
-L PvP01_02_v1:198112 \
-L PvP01_02_v1:373744 \
-L PvP01_02_v1:594798 \
-L PvP01_03_v1:127847 \
-L PvP01_03_v1:334738 \
-L PvP01_03_v1:782122 \
-L PvP01_04_v1:421012 \
-L PvP01_04_v1:514934 \
-L PvP01_04_v1:885624 \
-L PvP01_05_v1:192482 \
-L PvP01_05_v1:452277 \
-L PvP01_05_v1:701523 \
-L PvP01_06_v1:45794 \
-L PvP01_06_v1:278171 \
-L PvP01_06_v1:944771 \
-L PvP01_07_v1:754506 \
-L PvP01_07_v1:1020470 \
-L PvP01_07_v1:1211093 \
-L PvP01_08_v1:45083 \
-L PvP01_08_v1:53327 \
-L PvP01_08_v1:1591244 \
-L PvP01_09_v1:539410 \
-L PvP01_09_v1:1070402 \
-L PvP01_09_v1:1838632 \
-L PvP01_10_v1:351755 \
-L PvP01_10_v1:523718 \
-L PvP01_10_v1:1385673 \
-L PvP01_11_v1:245834 \
-L PvP01_11_v1:720439 \
-L PvP01_11_v1:1145363 \
-L PvP01_12_v1:94291 \
-L PvP01_12_v1:844166 \
-L PvP01_12_v1:1844936 \
-L PvP01_13_v1:162821 \
-L PvP01_13_v1:659592 \
-L PvP01_13_v1:1770129 \
-L PvP01_14_v1:743338 \
-L PvP01_14_v1:1911110 \
-L PvP01_14_v1:3004298 

#extract GT table
$gatk/gatk VariantsToTable -V VTN_Ampliseq_PV_Nov2022_panel_filteredPASS_DP5_ANN_VTN_SNPs.vcf.gz  -F CHROM -F POS -F TYPE -GF GT -O VTN_Ampliseq_PV_Nov2022_panel_filteredPASS_DP5_ANN_VTN_SNPs_GT.table


##########################################
#for VivaxGEO barcode AND Vietnam 42-SNP barcode
$gatk/gatk SelectVariants \
-R $ref \
-V VTN_Ampliseq_PV_Nov2022_panel_filteredPASS_DP5_ANN.vcf.gz \
--select-type-to-include SNP \
-O VTN_Ampliseq_PV_Nov2022_panel_filteredPASS_DP5_ANN_VivaxGEO_VTN_SNPs.vcf.gz \
-L PvP01_03_v1:220999 \
-L PvP01_05_v1:693449 \
-L PvP01_05_v1:1079461 \
-L PvP01_05_v1:1285366 \
-L PvP01_06_v1:646186 \
-L PvP01_07_v1:261070 \
-L PvP01_07_v1:595236 \
-L PvP01_08_v1:442363 \
-L PvP01_08_v1:879062 \
-L PvP01_08_v1:1420994 \
-L PvP01_09_v1:303291 \
-L PvP01_09_v1:448825 \
-L PvP01_09_v1:1101235 \
-L PvP01_09_v1:1884013 \
-L PvP01_10_v1:480601 \
-L PvP01_10_v1:483568 \
-L PvP01_10_v1:1131031 \
-L PvP01_11_v1:384043 \
-L PvP01_11_v1:503992 \
-L PvP01_11_v1:1137410 \
-L PvP01_11_v1:1226906 \
-L PvP01_11_v1:1448769 \
-L PvP01_11_v1:1868012 \
-L PvP01_12_v1:1116377 \
-L PvP01_12_v1:1400307 \
-L PvP01_12_v1:1463080 \
-L PvP01_12_v1:1988355 \
-L PvP01_13_v1:66121 \
-L PvP01_14_v1:344881 \
-L PvP01_14_v1:1160762 \
-L PvP01_14_v1:1229487 \
-L PvP01_14_v1:1258544 \
-L PvP01_14_v1:1266326 \
-L PvP01_01_v1:121166 \
-L PvP01_01_v1:164620 \
-L PvP01_01_v1:671781 \
-L PvP01_02_v1:198112 \
-L PvP01_02_v1:373744 \
-L PvP01_02_v1:594798 \
-L PvP01_03_v1:127847 \
-L PvP01_03_v1:334738 \
-L PvP01_03_v1:782122 \
-L PvP01_04_v1:421012 \
-L PvP01_04_v1:514934 \
-L PvP01_04_v1:885624 \
-L PvP01_05_v1:192482 \
-L PvP01_05_v1:452277 \
-L PvP01_05_v1:701523 \
-L PvP01_06_v1:45794 \
-L PvP01_06_v1:278171 \
-L PvP01_06_v1:944771 \
-L PvP01_07_v1:754506 \
-L PvP01_07_v1:1020470 \
-L PvP01_07_v1:1211093 \
-L PvP01_08_v1:45083 \
-L PvP01_08_v1:53327 \
-L PvP01_08_v1:1591244 \
-L PvP01_09_v1:539410 \
-L PvP01_09_v1:1070402 \
-L PvP01_09_v1:1838632 \
-L PvP01_10_v1:351755 \
-L PvP01_10_v1:523718 \
-L PvP01_10_v1:1385673 \
-L PvP01_11_v1:245834 \
-L PvP01_11_v1:720439 \
-L PvP01_11_v1:1145363 \
-L PvP01_12_v1:94291 \
-L PvP01_12_v1:844166 \
-L PvP01_12_v1:1844936 \
-L PvP01_13_v1:162821 \
-L PvP01_13_v1:659592 \
-L PvP01_13_v1:1770129 \
-L PvP01_14_v1:743338 \
-L PvP01_14_v1:1911110 \
-L PvP01_14_v1:3004298 

#extract GT 
$gatk/gatk VariantsToTable -V VTN_Ampliseq_PV_Nov2022_panel_filteredPASS_DP5_ANN_VivaxGEO_VTN_SNPs.vcf.gz -F CHROM -F POS -F TYPE -GF GT -O VTN_Ampliseq_PV_Nov2022_panel_filteredPASS_DP5_ANN_VivaxGEO_VTN_SNPs_GT.table



##################################
####Filter to biallic SNPs only and index file

$gatk/gatk SelectVariants \
-R $ref \
-V VTN_Ampliseq_PV_Nov2022_panel_filteredPASS_DP5_ANN.vcf.gz \
--select-type-to-include SNP \
--restrict-alleles-to BIALLELIC \
-O VTN_Ampliseq_PV_Nov2022_panel_filteredPASS_DP5_ANN_biSNPs.vcf.gz 

$gatk/gatk IndexFeatureFile \
-F VTN_Ampliseq_PV_Nov2022_panel_filteredPASS_DP5_ANN_biSNPs.vcf.gz

##END##
