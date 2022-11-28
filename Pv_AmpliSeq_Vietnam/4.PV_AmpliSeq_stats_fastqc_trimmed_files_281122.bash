#!/bin/bash
#Script to generate the sequencing stats of ampliseq data Vietnam Pv AmpliSeq after alignment, and to do fastqc after trimming
#version 28/11/2022 E Kattenberg for github


#set file paths:
#1 to the folder with the demultiplexed fastq-files from the AMpliSeq sequencing reaction on MiSeq
Thepathfq=/home/user/fq/
#2 To the folder where the output will go
Thepath=/home/user/project_name/

#3. To the list of samples to include in the variantcalling
samplelist=$Thepath/Sample_list_PV_VTN_experiment1.txt

#set paths to programs on my computer
sam=~/bin/samtools-1.9
gatk=~/bin/gatk-4.1.2.0
TMPDIR=/home/user/tmp

#set path to reference genome  
#this is not the latest version, but the one I used for the Pv AmpliSeq Vietnam data publised in Kattenberg et al. 2022 
#You can update to latest version, download from PlasmoDB: https://plasmodb.org/plasmo/app/downloads/Current_Release/PvivaxP01/
ref=/home/user/PvP01/PlasmoDB-46/PlasmoDB-46_PvivaxP01_Genome.fasta


#################################
#Generate sequencing statistics #
#################################

#go to the folder with alignments
cd $Thepath/alignments 

#function to collect sequencing metrics
 
statsPV(){
	java -Xmx2g -jar ~/bin/picard.jar CollectAlignmentSummaryMetrics \
	R=$ref \
	I=$f.PV_sort.bam \
	O=$f.PV_sort.stats.txt
}

#loop to run the function for each sample in the samplelist
for f in $(cat $samplelist); do
	statsPV "$f" & 
done

wait

#extract headers from one of the samples, in this case sample 2018_GL-KRP-088_S50. Change filename to match one of your samples.
#head -n 7 2018_GL-KRP-088_S50.PV_sort.stats.txt | tail -1  > $Thepath/PV_VTN_experiment1_Stats.txt

#loop for all files in the folder alignments with extension ".PV_sort.stats.txt"
FILES1=$Thepath/alignments/*.PV_sort.stats.txt
for h in $FILES1
do
	echo "$h" >> $Thepath/stats_PV_names_VTN_experiment1.txt
	head -n 10 $h | tail -1  >> $Thepath/PV_VTN_experiment1_Stats.txt
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