library(reshape)
library(msa)
library(ape)
library(seqinr)
library(bios2mds)
library(DataCombine)
library(GenomicAlignments)
source("Rscripts/baseRscript2.R")

bams<-list.files("Output/bam_PID/",pattern="^Run4.+.bam$")
#parameter setting for stackStringsFromBam
gr <- GRanges(seqnames="AY032751env",  IRanges(215, 691))

#Create an alignment file from bam files pos 215:691
for (i in 1:length(bams)){
    index<-paste0("Output/bam_PID/",bams[i],'.bai')
    bamfile <- BamFile(paste0("Output/bam_PID/",bams[i]), index=index)

    fname<-substr(bams[i], 1,7)
    seq<-stackStringsFromBam(bamfile, param=gr)
    
    #write the alignment fasta file
    writeXStringSet(seq,paste0("Output/PID_Con_Alignment/", fname, ".Alignment.fasta"),width =500)
}
    
#Replace '+' with '-' for the fasta files saved above in terminal
#     find PID_Con_Alignment -type f -exec sed -i "" 's/+/-/g' {} \;

