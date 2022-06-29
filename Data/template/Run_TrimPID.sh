#!/bin/bash

mkdir Run6_1_


#Trim reads from ends (R) only for  Score>15
bbduk.sh in=dir/XX/10_S10_L001_R1_001.fastq out=Run6_1_/Run6_1_Q15.R1.fastq qtrim=r trimq=15
bbduk.sh in=dir/XX/10_S10_L001_R2_001.fastq out=Run6_1_/Run6_1_Q15.R2.fastq qtrim=r trimq=15


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


bwa mem -t 8 -B 2 -O 2 -a SIV2 Run6_1_/Out/R1.fasta Run6_1_/Out/R2.fasta > Run6_1_/Run6_1__BWAmapped.sam

#Sort the sam file
samtools sort  Run6_1_/Run6_1__BWAmapped.sam -o  Sam/Run6_1_.sort.sam

rm Run6_1_/Run6_1__BWAmapped.sam Run6_1_/Out/R1.fasta-2 Run6_1_/Out/R2.fasta-2 Run6_1_/Out/r1_1.fasta Run6_1_/Out/r1_2.fasta Run6_1_/Out/r1_3.fasta Run6_1_/Out/r2_1.fasta Run6_1_/Out/r2_2.fasta Run6_1_/Out/r2_3.fasta
rm Run6_1_/Run6_1_Q15.R1.fastq Run6_1_/Run6_1_Q15.R2.fastq   

#Hard clip the sam file for merging Forward and Rev in R

java -jar ~/programs/jvarkit/dist/biostar84452.jar Sam/Run6_1_.sort.sam > Sam/Run6_1_.clipped.sam

