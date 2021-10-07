#!/bin/bash


mkdir Output/Run5_4_
#0 adapter trimming
bbduk.sh in1=~/programs/BaseSpace/B670_run_5_and_Mouse_LAI_RT-215809595/FASTQ_Generation_2020-12-12_02_55_33Z-352361009/4_/4_S4_L001_R1_001.fastq.gz in2=~/programs/Basespace/B670_run_5_and_Mouse_LAI_RT-215809595/FASTQ_Generation_2020-12-12_02_55_33Z-352361009/4_/4_S4_L001_R2_001.fastq.gz  out=Output/Run5_4_/Run5_4__adp.trimmed.fastq ref=/Users/kahotisthammer/programs/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 stats=Output/Run5_4_/stats_0.txt

#1 Trim reads at both ends at Score<30
bbduk.sh in=Output/Run5_4_/Run5_4__adp.trimmed.fastq out=Output/Run5_4_/Run5_4__trimmed.q30.fastq qtrim=rl trimq=30 stats=Output/Run5_4_/stats30_1.txt

#2. Kmer filtering
bbduk.sh in=Output/Run5_4_/Run5_4__trimmed.q30.fastq out=Output/Run5_4_/Run5_4__unmatched.q30.fq ref=/Users/kahotisthammer/programs/bbmap/resources/phix174_ill.ref.fa k=31 hdist=1 stats=Output/Run5_4_/stats30_2.txt

#3.Remove reads with ave score <30
#bbduk.sh in=Output/Run5_4_/Run5_4__unmatched.q30.fq out=Output/Run5_4_/Run5_4__clean.q30.fq maq=30 

#4. Align the file using bwa to the reference 

bwa mem -t 8 -k 15 -a SIV2 Output/Run5_4_/Run5_4__unmatched.q30.fq  > Output/Run5_4_/Run5_4__BWAmapped.sam

#5. convert sam to bam
samtools view -S -b Output/Run5_4_/Run5_4__BWAmapped.sam > Output/Run5_4_/Run5_4__BWAmapped.bam


#6. sort the bam file
samtools sort  Output/Run5_4_/Run5_4__BWAmapped.bam -o  Output/bam/Run5_4__BWA_sort.bam

#7. index the bam file
samtools index  Output/bam/Run5_4__BWA_sort.bam  Output/bam/Run5_4__BWA_sort.bam.bai

rm Output/Run5_4_/Run5_4__trimmed.q30.fastq Output/Run5_4_/Run5_4__unmatched.q30.fq Output/Run5_4_/Run5_4__BWAmapped.sam

