#!/bin/bash
# Script for the data preparation of ampliseq data Peru AmpliSeq and alignment to the reference genome 
# version 28/11/2022 E Kattenberg for github

#set file paths: 
#1 to the folder with the demultiplexed fastq-files from the AMpliSeq sequencing reaction on MiSeq
Thepathfq=/home/user/fq/
#2 To the folder where the output will go
Thepath=/home/user/project_name/

#3. To the list of samples to include in the preprocessing and alignment
samplelist=$Thepath/Sample_list_PF_Peru_experiment1.txt
#I made the sample list in the terminal using the following code:
#cd $Thepathfq
#echo *_R1_001.fastq.gz >> $Thepath/Sample_list_PF_Peru_experiment1.txt
#adjust to remove all _R1_001.fastq.gz and save
#This way you have a nice record of the samples that were included in the analysis and can be used to make the metadata file


#set paths to programs on my device
trim=~/bin/Trimmomatic-0.39
bwa=~/bin/bwa-0.7.17
sam=~/bin/samtools-1.9
gatk=~/bin/gatk-4.1.2.0
TMPDIR=/home/user/tmp

#set path to reference genome  
#this is not the latest version, but the one I used for the Pf AmpliSeq Peru data publised in Kattenberg et al. 2022 
#You can update to latest version, download from PlasmoDB: https://plasmodb.org/plasmo/app/downloads/Current_Release/Pfalciparum3D7/
ref=/home/user/3D7_genome/PlasmoDB-44_Pfalciparum3D7_Genome.fasta

########################################################################
#1. check quality of reads in fastq files with FastQC and FastQ Screen #
########################################################################

#make folder to store the files
cd $Thepath/
mkdir fastqc_results
mkdir fastQscreen_res

#Go to folder with FastQC program
cd ~/bin/fastqc_v0.11.8/FastQC

#run fastqc for the fastq of read 1 and read 2 for each sample
./fastqc $Thepathfq/*_R1_001.fastq.gz --outdir=$Thepath/fastqc_results 
./fastqc $Thepathfq/*_R2_001.fastq.gz --outdir=$Thepath/fastqc_results 

#Go to folder with FastQ-Screen program
cd ~/bin/FastQ-Screen-0.14.1

#run fastq screen for the fastq of read 1 and read 2 for each sample. In fastqscreen we check for reads aligning to P. falciparuma and the human genome. We routinely also check for the other human infecting malaria species (Pv, Pm, Pk, Po)
./fastq_screen $Thepathfq/*_R1_001.fastq.gz --outdir=$Thepath/fastQscreen_res --aligner bwa
./fastq_screen $Thepathfq/*_R2_001.fastq.gz --outdir=$Thepath/fastQscreen_res --aligner bwa



#######################################
# 2. TRIMMOMATIC					  #
# trim adapters and low quality reads #
#######################################

#set path to file adapter sequences (see adpater sequences on illumina website)
adapters=$Thepath/adapters.fa
# this part isnt essential, as the MiSeq during the demultiplexing already removes adapter reads. 



#make folder to store the files
cd $Thepath/
mkdir fq_trim

#Go to the folder with the fastq files
cd $Thepathfq

#set the function for trimmomatic using the filenames in the folder
trimtask(){
	r1=$Thepathfq/${f}_L001_R1_001.fastq.gz
	r2=$Thepathfq/${f}_L001_R2_001.fastq.gz
	
	java -jar $trim/trimmomatic-0.39.jar PE \
		-phred33 \
		-trimlog $Thepath/fq_trim/${f}_trimreport.txt \
		$r1 \
		$r2 \
		$Thepath/fq_trim/${f}_paired_R1.fq.gz \
		$Thepath/fq_trim/${f}_unpaired_R1.fq.gz \
		$Thepath/fq_trim/${f}_paired_R2.fq.gz \
		$Thepath/fq_trim/${f}_unpaired_R2.fq.gz \
	 	ILLUMINACLIP:$adapters:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
}

#loop to run the function for each sample in the samplelist
for f in $(cat $samplelist); do
	trimtask "$f" 2>file_$f.err &
done

#wait till trimmomatic is done for all files before continuing
wait

#######################################
# 3. Burrows-Wheeler Alignment		  #
# alignment to 3D7 reference sequence #
#######################################

#make folder to store the files
cd $Thepath/
mkdir alignments

#Go to folder with Burrows-Wheeler aligner program
cd ~/bin/bwa-0.7.17/

#first time only: need to index for the reference sequence:
$bwa/bwa index $ref

#set function for Pf alignment task
task(){
	r1=$Thepath/fq_trim/${f}_paired_R1.fq.gz
	r2=$Thepath/fq_trim/${f}_paired_R2.fq.gz

	$bwa/bwa mem -M -R "@RG\tID:L001\tPL:ILLUMINA\tSM:${f}" -t 1\
	$ref \
	$r1 \
	$r2 \
	> $Thepath/alignments/${f}.PF_algn.sam
}

#loop to run the function for each sample in the samplelist, but for 6 samples at a time (1 sample/core) to make it run more smoothly on my computer
N=6
(
for f in $(cat $samplelist); do
	((i=i%N)); ((i++==0))&& wait
	task "$f" & 
done
)

wait

##END##