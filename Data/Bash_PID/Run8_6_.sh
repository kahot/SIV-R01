
#1. Map the consensus fasta files to ref 
bwa mem -t 8 -B 0 -O 0 -L 50 -w 200 SIVB670 Output/PID_Consensus/Run8_6_.consensus.fasta  > Output/sam/Run8_6_.PID.BWA.sam

#bbmap.sh in=Output/PID_Consensus/Run8_6_.consensus.fasta out=Output/bam_PID/Run8_6_.PID.bbmap.sam

#2. convert sam to bam
samtools view -S -b Output/sam/Run8_6_.PID.BWA.sam > Output/bam_PID/Run8_6_.PID.BWA.bam


#3. sort the bam file
samtools sort  Output/bam_PID/Run8_6_.PID.BWA.bam -o  Output/bam_PID/Run8_6_.PID.sort.bam

#4. index the bam file
samtools index  Output/bam_PID/Run8_6_.PID.sort.bam  Output/bam_PID/Run8_6_.PID.sort.bam.bai

rm Output/bam_PID/Run8_6_.PID.BWA.bam 

