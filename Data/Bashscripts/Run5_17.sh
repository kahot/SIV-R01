#!/bin/bash


mkdir Output/Run5_17
#0 adapter trimming
bbduk.sh in1=~/programs/BaseSpace/B670_run_5_and_Mouse_LAI_RT-215809595/FASTQ_Generation_2020-12-12_02_55_33Z-352361009/17/17_S17_L001_R1_001.fastq.gz in2=~/programs/Basespace/B670_run_5_and_Mouse_LAI_RT-215809595/FASTQ_Generation_2020-12-12_02_55_33Z-352361009/17/17_S17_L001_R2_001.fastq.gz  out=Output/Run5_17/Run5_17_adp.trimmed.fastq ref=/Users/kahotisthammer/programs/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 stats=Output/Run5_17/stats_0.txt

#1 Trim reads at both ends at Score<30
bbduk.sh in=Output/Run5_17/Run5_17_adp.trimmed.fastq out=Output/Run5_17/Run5_17_trimmed.q30.fastq qtrim=rl trimq=30 stats=Output/Run5_17/stats30_1.txt

#2. Kmer filtering
bbduk.sh in=Output/Run5_17/Run5_17_trimmed.q30.fastq out=Output/Run5_17/Run5_17_unmatched.q30.fq ref=/Users/kahotisthammer/programs/bbmap/resources/phix174_ill.ref.fa k=31 hdist=1 stats=Output/Run5_17/stats30_2.txt

#3.Remove reads with ave score <30
#bbduk.sh in=Output/Run5_17/Run5_17_unmatched.q30.fq out=Output/Run5_17/Run5_17_clean.q30.fq maq=30 

#4. Align the file using bwa to the reference 

bwa mem -t 8 -k 15 -a SIV2 Output/Run5_17/Run5_17_unmatched.q30.fq  > Output/Run5_17/Run5_17_BWAmapped.sam

#5. convert sam to bam
samtools view -S -b Output/Run5_17/Run5_17_BWAmapped.sam > Output/Run5_17/Run5_17_BWAmapped.bam


#6. sort the bam file
samtools sort  Output/Run5_17/Run5_17_BWAmapped.bam -o  Output/bam/Run5_17_BWA_sort.bam

#7. index the bam file
samtools index  Output/bam/Run5_17_BWA_sort.bam  Output/bam/Run5_17_BWA_sort.bam.bai

rm Output/Run5_17/Run5_17_trimmed.q30.fastq Output/Run5_17/Run5_17_unmatched.q30.fq Output/Run5_17/Run5_17_BWAmapped.sam

