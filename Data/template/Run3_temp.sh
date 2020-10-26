#!/bin/bash

mkdir Run4_18

#mv fastq file
cp -r ~/programs/BaseSpace/B670_Run_3-142079938/FASTQ_Generation_2019-10-03_13_19_45Z-193695834/XX/*.gz Run4_18


#Run PID script
ruby TCS_SIV.rb Run4_18

mkdir Run4_18/Out

#Rename the r1, r2 files and move to Output directory.
mv Run4_18/Reverse1/consensus/r1.txt Run4_18/Out/r1_1.fasta
mv Run4_18/Reverse2/consensus/r1.txt Run4_18/Out/r1_2.fasta
mv Run4_18/Reverse3/consensus/r1.txt Run4_18/Out/r1_3.fasta
mv Run4_18/Reverse1/consensus/r2.txt Run4_18/Out/r2_1.fasta
mv Run4_18/Reverse2/consensus/r2.txt Run4_18/Out/r2_2.fasta
mv Run4_18/Reverse3/consensus/r2.txt Run4_18/Out/r2_3.fasta

#combine the files

cat Run4_18/Out/r1_1.fasta Run4_18/Out/r1_2.fasta Run4_18/Out/r1_3.fasta > Run4_18/Out/R1.fasta
cat Run4_18/Out/r2_1.fasta Run4_18/Out/r2_2.fasta Run4_18/Out/r2_3.fasta > Run4_18/Out/R2.fasta

# Remove r1 and r2  from the sequence names

sed -i -2 's/_r1//g' Run4_18/Out/R1.fasta
sed -i -2 's/_r2//g' Run4_18/Out/R2.fasta


bwa mem -t 8 -k 15 -a SIV2 Run4_18/Out/R1.fasta Run4_18/Out/R2.fasta > Run4_18/Run4_18_BWAmapped.sam

samtools view -S -b Run4_18/Run4_18_BWAmapped.sam > Bam/Run4_18_BWAmapped.bam

samtools sort  Bam/Run4_18_BWAmapped.bam -o  Bam/Run4_18_BWAmapped_sort.bam

#7. index the bam file
samtools index  Bam/Run4_18_BWAmapped_sort.bam  Bam/Run4_18_BWAmapped_sort.bam.bai

rm Run4_18/Run4_18_BWAmapped.sam Bam/Run4_18_BWAmapped.bam

