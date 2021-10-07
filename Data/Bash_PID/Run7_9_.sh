
#1. Map the consensus fasta files to ref 
bwa mem -t 8 -B 0 -O 0 -L 50 -w 200 SIVB670 Output/PID_Consensus/Run7_9_.consensus.fasta  > Output/bam_PID/Run7_9_.PID.BWA.sam

#bbmap.sh in=Output/PID_Consensus/Run7_9_.consensus.fasta out=Output/bam_PID/Run7_9_.PID.bbmap.sam

#2. convert sam to bam
samtools view -S -b Output/bam_PID/Run7_9_.PID.BWA.sam > Output/bam_PID/Run7_9_.PID.BWA.bam


#3. sort the bam file
samtools sort  Output/bam_PID/Run7_9_.PID.BWA.bam -o  Output/bam_PID/Run7_9_.PID.sort.bam

#4. index the bam file
samtools index  Output/bam_PID/Run7_9_.PID.sort.bam  Output/bam_PID/Run7_9_.PID.sort.bam.bai

rm Output/bam_PID/Run7_9_.PID.BWA.bam 

