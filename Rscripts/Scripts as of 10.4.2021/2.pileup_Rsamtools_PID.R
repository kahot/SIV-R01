library(Rsamtools)
library(stringr)

source("Rscripts/pileupFreq.R")

#number of samples to process
bamfiles<-list.files("Output/bam_PID/",pattern="bam$")
bamfiles<-list.files("Output/bam_PID/",pattern="^Run7_.*bam$")


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
    write.csv(summary, file=paste0("Output/CSV_PIDcon/",file.name,".PID.csv",collapse=""))
    
}



### mapped PID output R1 and R2 to the reference without merging ###
#bamfiles<-list.files("~/programs/PID-master/Bam/",pattern="bam$")
##dir.create("Output/CSV_PID/")
#
#for (i in 1:length(bamfiles)){
#        bam<-bamfiles[i]
#        index<-paste0("~/programs/PID-master/Bam/",bam,'.bai')
#        bf<-BamFile(paste0("~/programs/PID-master/Bam/",bam), index=index)
#
#        file.name<-paste(bam)
#        file.name<-substr(file.name,start=1,stop=7 )
#        p_param <- PileupParam(max_depth=50000,include_insertions=TRUE,include_deletions=TRUE)
#        result<-pileup(bf, pileupParam = p_param, distinguish_strands=FALSE,ignore_query_Ns=FALSE)
#        summary<-pileupFreq(result)
#
#        summary$TotalReads<-rowSums(summary[3:6])
#        maxr<-max(summary$TotalReads)
#        print(file.name)
#        cat("The maximum number of read depth is ", maxr)
#        cat("\n")
#        write.csv(summary, file=paste0("Output/CSV_PID/",file.name,".csv",collapse=""))
#        
#  }
