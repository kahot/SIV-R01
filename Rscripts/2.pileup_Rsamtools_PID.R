library(Rsamtools)
library(stringr)

source("Rscripts/pileupFreq.R")

#number of samples to process
bamfiles<-list.files("Output/bam_PID/",pattern="bam$")
bamfiles<-list.files("Output/bam_PID/",pattern="^Run8_.*bam$")


#dir.create("Output/CSV_PIDcon/")
#for PID consensus mapping
for (i in 1:length(bamfiles)){
    bam<-bamfiles[i]
    index<-paste0("Output/bam_PID/",bam,'.bai')
    bf<-BamFile(paste0("Output/bam_PID/",bam), index=index)
    
    file.name<-paste(bam)
    file.name<-substr(file.name,start=1,stop=7 )
    p_param <- PileupParam(max_depth=50000,include_insertions=TRUE,include_deletions=TRUE)
    result<-pileup(bf, pileupParam = p_param, distinguish_strands=FALSE,ignore_query_Ns=FALSE)
    summary<-pileupFreq(result)
    
    summary$TotalReads<-rowSums(summary[3:6])
    maxr<-max(summary$TotalReads)
    print(file.name)
    cat("The maximum number of read depth is ", maxr)
    cat("\n")
    write.csv(summary, file=paste0("Output/CSV_PID/",file.name,".PID.csv",collapse=""))
}


