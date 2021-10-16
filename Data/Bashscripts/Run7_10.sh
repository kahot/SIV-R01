#!/bin/bash


mkdir Output/Run7_10
#0 adapter trimming
bbduk.sh in1=~/programs/BaseSpace/B670_Run_7_Repeat_2-292249961/FASTQ_Generation_2021-09-07_14_35_26Z-457696525/10/10_S10_L001_R1_001.fastq.gz in2=~/programs/BaseSpace/B670_Run_7_Repeat_2-292249961/FASTQ_Generation_2021-09-07_14_35_26Z-457696525/10/10_S10_L001_R2_001.fastq.gz  out=Output/Run7_10/Run7_10_adp.trimmed.fastq ref=/Users/kahotisthammer/programs/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 stats=Output/Run7_10/stats_0.txt

#1 Trim reads at both ends at Score<30
bbduk.sh in=Output/Run7_10/Run7_10_adp.trimmed.fastq out=Output/Run7_10/Run7_10_trimmed.q30.fastq qtrim=rl trimq=30 stats=Output/Run7_10/stats30_1.txt

#2. Kmer filtering
bbduk.sh in=Output/Run7_10/Run7_10_trimmed.q30.fastq out=Output/Run7_10/Run7_10_unmatched.q30.fq ref=/Users/kahotisthammer/programs/bbmap/resources/phix174_ill.ref.fa k=31 hdist=1 stats=Output/Run7_10/stats30_2.txt

#3.Remove reads with ave score <30
#bbduk.sh in=Output/Run7_10/Run7_10_unmatched.q30.fq out=Output/Run7_10/Run7_10_clean.q30.fq maq=30 

#4. Align the file using bwa to the reference 

bwa mem -t 8 -k 15 -a SIV2 Output/Run7_10/Run7_10_unmatched.q30.fq  > Output/Run7_10/Run7_10_BWAmapped.sam

#5. convert sam to bam
samtools view -S -b Output/Run7_10/Run7_10_BWAmapped.sam > Output/Run7_10/Run7_10_BWAmapped.bam


#6. sort the bam file
samtools sort  Output/Run7_10/Run7_10_BWAmapped.bam -o  Output/bam/Run7_10_BWA_sort.bam

#7. index the bam file
samtools index  Output/bam/Run7_10_BWA_sort.bam  Output/bam/Run7_10_BWA_sort.bam.bai

rm Output/Run7_10/Run7_10_trimmed.q30.fastq Output/Run7_10/Run7_10_unmatched.q30.fq Output/Run7_10/Run7_10_BWAmapped.sam

