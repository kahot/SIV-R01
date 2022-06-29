#!/bin/bash


mkdir Output/M10
#0 adapter trimming
bbduk.sh in1=dir1/dir2/XX/10_S10_L001_R1_001.fastq in2=dir1/dir2/XX/10_S10_L001_R2_001.fastq  out=Output/M10/M10_adp.trimmed.fastq ref=/Users/kahotisthammer/programs/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 stats=Output/M10/stats_0.txt

#Trim reads at both ends at Score>20
bbduk.sh in=Output/M10/M10_adp.trimmed.fastq out=Output/M10/M10_trimmed.q20.fastq qtrim=rl trimq=20 stats=Output/M10/stats20_1.txt
#Trim reads at both ends at Score>15
bbduk.sh in=Output/M10/M10_adp.trimmed.fastq out=Output/M10/M10_trimmed.q15.fastq qtrim=rl trimq=14 stats=Output/M10/stats14_1.txt


#Kmer filtering
bbduk.sh in=Output/M10/M10_trimmed.q20.fastq out=Output/M10/M10_unmatched.q20.fq ref=/Users/kahotisthammer/programs/bbmap/resources/phix174_ill.ref.fa k=31 hdist=1 stats=Output/M10/stats20_2.txt
bbduk.sh in=Output/M10/M10_trimmed.q15.fastq out=Output/M10/M10_unmatched.q15.fq ref=/Users/kahotisthammer/programs/bbmap/resources/phix174_ill.ref.fa k=31 hdist=1 stats=Output/M10/stats15_2.txt

#Align the file using bwa to the reference 
bwa mem -t 8 -k 15 -a SIV2 Output/M10/M10_unmatched.q20.fq  > Output/M10/M10_20mapped.sam
bwa mem -t 8 -k 15 -a SIV2 Output/M10/M10_unmatched.q15.fq  > Output/M10/M10_15mapped.sam

#5. convert sam to bam
samtools view -S -b Output/M10/M10_20mapped.sam > Output/M10/M10_20mapped.bam
samtools view -S -b Output/M10/M10_15mapped.sam > Output/M10/M10_15mapped.bam


#6. sort the bam file
samtools sort  Output/M10/M10_20mapped.bam -o  Output/bam/M10_20_sort.bam
samtools sort  Output/M10/M10_15mapped.bam -o  Output/bam/M10_15_sort.bam

#7. index the bam file
samtools index  Output/bam/M10_20_sort.bam  Output/bam/M10_20_sort.bam.bai
samtools index  Output/bam/M10_15_sort.bam  Output/bam/M10_15_sort.bam.bai

rm Output/M10/M10_trimmed.q20.fastq Output/M10/M10_trimmed.q15.fastq Output/M10/M10_unmatched.q20.fq Output/M10/M10_unmatched.q15.fq Output/M10/M10_20mapped.sam Output/M10/M10_15mapped.sam Output/M10/M10_20mapped.bam Output/M10/M10_15mapped.bam

