#!/bin/bash
#Script to sort and index sam and bam files and do the variantcalling and generate the stats and samplemap of ampliseq data Peru AmpliSeq after alignment
#version 28/11/2022 E Kattenberg for github


#set file paths:
#1 to the folder with the demultiplexed fastq-files from the AMpliSeq sequencing reaction on MiSeq
Thepathfq=/home/user/fq/
#2 To the folder where the output will go
Thepath=/home/user/project_name/

#3. To the list of samples to include in the variantcalling
samplelist=$Thepath/Sample_list_PF_Peru_experiment1.txt

#set paths to programs on my computer
sam=~/bin/samtools-1.9
gatk=~/bin/gatk-4.1.2.0
TMPDIR=/home/user/tmp

#set path to reference genome  
#this is not the latest version, but the one I used for the Pf AmpliSeq Peru data publised in Kattenberg et al. 2022 
#You can update to latest version, download from PlasmoDB: https://plasmodb.org/plasmo/app/downloads/Current_Release/Pfalciparum3D7/
ref=/home/user/3D7_genome/PlasmoDB-44_Pfalciparum3D7_Genome.fasta


###################################################
#4. Alignment preprocessing						  #
#Convert to bam and sort, then index the bam file #
###################################################

#Go to the folder with the bwa alignments created with the 1st script
cd $Thepath/alignments

#1. Function to sort SAM file (named in the format *.PF_algn.bam) and convert to BAM using PICARD TOOLS SortSam.jar

sortsam(){
	java -Xmx5g -jar ~/bin/picard.jar SortSam INPUT=${f}.PF_algn.sam OUTPUT=$f.PF_sort.bam 		SORT_ORDER=coordinate TMP_DIR=$TMPDIR VALIDATION_STRINGENCY=SILENT
}

#loop to run the function for each sample in the samplelist
for f in $(cat $samplelist); do
	sortsam "$f" & 
done

wait


# 2. Set the function to index the bam files with samtoolsindex
index(){
	$sam/samtools index $f.PF_sort.bam
}

#loop to run the function for each sample in the samplelist, but for 6 samples at a time (1 sample/core)
N=6
(
for f in $(cat $samplelist); do
	index "$f" & 
done   
)
wait

#####################################
# 5. Haplotypecaller				#
# Variant calling   				#
#####################################

#Function to do the alignment for every sample individually using haplotypecaller in gvcf mode
variantsPF(){
	$gatk/gatk HaplotypeCaller \
     	-R $ref \
     	-I $f.PF_sort.bam \
     	--emit-ref-confidence GVCF \
	-O $f.PF_sort.GATK.g.vcf
}

#Go to the folder with the alignments in sorted bam files
cd $Thepath/alignments 

#loop to run the function for each sample in the samplelist, but for 6 samples at a time (1 sample/core)
N=6
(
for f in $(cat $samplelist); do
	((i=i%N)); ((i++==0))&& wait
	variantsPF "$f" & 
done
)

wait


##########################################################
# 6. Preparation of samplemap file for Joint genotyping  #
##########################################################

#make folder to store the samplemap file (and later the joint database)
cd $Thepath/alignments
mkdir database

#go to new folder
cd $Thepath/alignments/database

#Make samplemap:
#1. add names of PF alignment for samples from sample list
#2. add link to the vcf file (created in the previous step) for that sample 

for f in $(cat $samplelist);
do
	echo "$f	$Thepath/alignments/$f.PF_sort.GATK.g.vcf" >> $Thepath/alignments/database/samplemap_PF_experiment1.txt

done


#sort the samplemap
sort -k 1 samplemap_PF_experiment1.txt >> samplemap_PF_experiment1_sort.txt 

##Note: up until this point we have done everything for the samples individually, and I usually run the scripts 1,2 and 4 after each sequencing run.
# But for the joint genotyping I usually combine multiple RUNS, including all samples for a whole project together. Therefore, you will need to combine all the samples into one samplemap 


##END##