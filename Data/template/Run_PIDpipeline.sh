#!/bin/bash

mkdir Run6_1_

#mv fastq file
cp -r /Users/kahotisthammer/programs/BaseSpace/B670_Run_6-268458190/FASTQ_Generation_2021-06-04_22_44_27Z-424462038/XX/*.gz Run6_1_


#Run PID script
ruby TCS_SIV.rb Run6_1_

mkdir Run6_1_/Out

cd Run6_1_

#Rename the r1, r2 files and move to Output directory.
mv Reverse1/consensus/r1.txt Out/r1_1.fasta
mv Reverse2/consensus/r1.txt Out/r1_2.fasta
mv Reverse3/consensus/r1.txt Out/r1_3.fasta
mv Reverse1/consensus/r2.txt Out/r2_1.fasta
mv Reverse2/consensus/r2.txt Out/r2_2.fasta
mv Reverse3/consensus/r2.txt Out/r2_3.fasta

cd ..


#combine the files

cat Run6_1_/Out/r1_1.fasta Run6_1_/Out/r1_2.fasta Run6_1_/Out/r1_3.fasta > Run6_1_/Out/R1.fasta
cat Run6_1_/Out/r2_1.fasta Run6_1_/Out/r2_2.fasta Run6_1_/Out/r2_3.fasta > Run6_1_/Out/R2.fasta

# Remove r1 and r2  from the sequence names

sed -i -2 's/_r1//g' Run6_1_/Out/R1.fasta
sed -i -2 's/_r2//g' Run6_1_/Out/R2.fasta


bwa mem -t 8 -B 1 -O 1 -a SIV2 Run6_1_/Out/R1.fasta Run6_1_/Out/R2.fasta > Run6_1_/Run6_1__BWAmapped.sam

#Sort the sam file
samtools sort  Run6_1_/Run6_1__BWAmapped.sam -o  Sam/Run6_1_.sort.sam

rm Run6_1_/Run6_1__BWAmapped.sam Run6_1_/Out/R1.fasta-2 Run6_1_/Out/R2.fasta-2 Run6_1_/*.fastq

#Hard clip the sam file for merging Forward and Rev in R

java -jar ~/programs/jvarkit/dist/biostar84452.jar Sam/Run6_1_.sort.sam > Sam/Run6_1_.clipped.sam

