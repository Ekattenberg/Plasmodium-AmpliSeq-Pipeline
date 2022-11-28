#!/bin/bash
#Script to generate the sequencing stats of ampliseq data Peru AmpliSeq after alignment, and to do fastqc after trimming
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



#################################
#Generate sequencing statistics #
#################################

#go to the folder with alignments
cd $Thepath/alignments 

#function to collect sequencing metrics
 
statsPF(){
	java -Xmx2g -jar ~/bin/picard.jar CollectAlignmentSummaryMetrics \
	R=$ref \
	I=$f.PF_sort.bam \
	O=$f.PF_sort.stats.txt
}

#loop to run the function for each sample in the samplelist
for f in $(cat $samplelist); do
	statsPF "$f" & 
done   

wait

#extract headers from one of the samples, in this case sample 1_S1. Change filename to match one of your samples.
#head -n 7 1_S1.PF_sort.stats.txt | tail -1  > $Thepath/PF_experiment1_Stats.txt

#loop for all files in the folder alignments with extension ".PF_sort.stats.txt"
FILES1=$Thepath/alignments/*.PF_sort.stats.txt
for h in $FILES1
do
	echo "$h" >> $Thepath/stats_PF_names_PF_experiment1.txt
	head -n 10 $h | tail -1  >> $Thepath/PF_experiment1_Stats.txt
done

###############################################################
#Check quality of reads in trimmed fastq files with FastQC #
###############################################################

#this gives an idea of the effect of the trimming, not essential for the pipeline
#make folder for output
cd $Thepath/
mkdir fastqc_trim

#Go to folder with FastQC program
cd ~/bin/fastqc_v0.11.8/FastQC

#run fastqc for the fastq of the trimmed read 1 and read 2 for each sample
./fastqc $Thepath/fq_trim/*_paired_R1.fq.gz --outdir=$Thepath/fastqc_trim 
./fastqc $Thepath/fq_trim/*_paired_R2.fq.gz --outdir=$Thepath/fastqc_trim 

##END##