# Plasmodium-AmpliSeq-Pipeline
This github contains an in-house linux pipeline used for AmpliSeq alignment and variant calling for P. falciparum (Kattenberg et al. accepted for publication in Microbiology Spectrum) or P. vivax (Kattenberg et al. Front. Cell. Infect. Microbiol. 2022) AmpliSeq panels for molecular surveillance of malaria. The purpose of this repository is to be transparent in the data processing that was performed and to allow others to replicate the AmpliSeq data processing in the same way. 

For more information on the NGS assay designs for which we used these pipelines, please see:
1. Kattenberg JH, Fernandez-Miñope C, van Dijk NJ, Llacsahuanga Allcca L, Guetens P, Valdivia HO, Van Geertruyden JP, Rovira-Vallbona E, Monsieurs P, Delgado-Ratto C, Gamboa D, Rosanas-Urgell A. (2021) Malaria molecular surveillance in the Peruvian Amazon with novel highly multiplexed Plasmodium falciparum Ampliseq assay. medRxiv https://doi.org/10.1101/2021.11.12.21266245

2. Kattenberg JH, Nguyen HV, Nguyen HL, Sauve E, Nguyen NTH, Chopo-Pizarro A, Trimarsanto H, Monsieurs P, Guetens P, Nguyen XX, Esbroeck MV, Auburn S, Nguyen BTH and Rosanas-Urgell A (2022) Novel highly-multiplexed AmpliSeq targeted assay for Plasmodium vivax genetic surveillance use cases at multiple geographical scales. Front. Cell. Infect. Microbiol. 12:953187. https://doi.org/10.3389/fcimb.2022.953187

We have different versions of the AmpliSeq Custom targeted sequencing assay, for P. falciparum and P. vivax. Currently, the scripts included here were those used for the Pf AmpliSeq for Peru and for the Pv AmpliSeq for Vietnam, which are described in the publications above. The scripts can be found in the folders above 

In brief the different steps for the variant calling are given in 6 bash scripts:
- Script 1. Trimming poor quality reads and alignment, incl. quality and contamination analysis 
- Script 2. Variant calling with Haplotypecaller in gvcf mode
- Script 3. Joint genotyping
- Script 4. Alignment statistics 
- Script 5. Hard filtering variants and annotation
- Script 6. Subsetting variants

Several programs are required in order to run the scripts that should be downloaded elsewhere (see links) and installed on your computer: 
- Burrows-Wheeler aligner: https://bio-bwa.sourceforge.net/
- GATK 4: https://gatk.broadinstitute.org/
- Picard tools: https://github.com/broadinstitute/picard
- Samtools: http://www.htslib.org/
- SnpEff: http://pcingola.github.io/SnpEff/
- Trimmomatic: http://www.usadellab.org/cms/?page=trimmomatic
- FastQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
- FastQ Screen: https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/

Note: Using different versions of the programs than specified in the scripts provided here might result in differences in the code required to run the programs. Adjust the scripts accordingly.  
