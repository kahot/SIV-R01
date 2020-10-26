#!/bin/bash

mkdir Run4_18

#mv fastq file
cp -r /Users/kahotisthammer/programs/BaseSpace/Ambrose1-41042018/FASTQ_Generation_2018-05-21_03_30_27Z-96699097/XX/*.gz Run4_18


#Run PID script
ruby TCS_SIV_run2.rb Run4_18

mkdir Run4_18/Out

#Rename the r1, r2 files and move to Output directory.
mv Run4_18/ENV1/consensus/r1.txt Run4_18/Out/R1.fasta
mv Run4_18/ENV1/consensus/r2.txt Run4_18/Out/R2.fasta

# Remove r1 and r2  from the sequence names

sed -i -2 's/_r1//g' Run4_18/Out/R1.fasta
sed -i -2 's/_r2//g' Run4_18/Out/R2.fasta


bwa mem -t 8 -k 15 -a SIV2 Run4_18/Out/R1.fasta Run4_18/Out/R2.fasta > Run4_18/Run4_18_BWAmapped.sam

samtools view -S -b Run4_18/Run4_18_BWAmapped.sam > Bam/Run4_18_BWAmapped.bam

samtools sort  Bam/Run4_18_BWAmapped.bam -o  Bam/Run4_18_BWAmapped_sort.bam

#7. index the bam file
samtools index  Bam/Run4_18_BWAmapped_sort.bam  Bam/Run4_18_BWAmapped_sort.bam.bai

rm Run4_18/Run4_18_BWAmapped.sam Bam/Run4_18_BWAmapped.bam
