#!/bin/bash


mkdir Output/M10
#0 adapter trimming
bbduk.sh in1=~/programs/Basespace/10_S10_L001_R1_001.fastq in2=~/programs/Basespace/10_S10_L001_R2_001.fastq  out=Output/M10/M10_adp.trimmed.fastq ref=/Users/kahotisthammer/programs/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 stats=Output/M10/stats_0.txt

#1 Trim reads at both ends at Score<35
bbduk.sh in=Output/M10/M10_adp.trimmed.fastq out=Output/M10/M10_trimmed.q35.fastq qtrim=rl trimq=35 stats=Output/M10/stats35_1.txt

#2. Kmer filtering
bbduk.sh in=Output/M10/M10_trimmed.q35.fastq out=Output/M10/M10_unmatched.q35.fq ref=/Users/kahotisthammer/programs/bbmap/resources/phix174_ill.ref.fa k=31 hdist=1 stats=Output/M10/stats35_2.txt

#3.Remove reads with ave score <30
#bbduk.sh in=Output/M10/M10_unmatched.q30.fq out=Output/M10/M10_clean.q30.fq maq=30 

#4. Align the file using bwa to the reference 

bwa mem -t 8 -k 15 -a SIV2 Output/M10/M10_unmatched.q35.fq  > Output/M10/M10_BWAmapped.sam

#5. convert sam to bam
samtools view -S -b Output/M10/M10_BWAmapped.sam > Output/M10/M10_BWAmapped.bam


#6. sort the bam file
samtools sort  Output/M10/M10_BWAmapped.bam -o  Output/bam/M10_BWA_sort.bam

#7. index the bam file
samtools index  Output/bam/M10_BWA_sort.bam  Output/bam/M10_BWA_sort.bam.bai

rm Output/M10/M10_trimmed.q35.fastq Output/M10/M10_unmatched.q35.fq Output/M10/M10_BWAmapped.sam

