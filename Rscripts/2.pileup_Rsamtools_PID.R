library(Rsamtools)
library(stringr)

source("Rscripts/pileupFreq.R")

#number of sampels to process
bamfiles<-list.files("~/programs/PID-master/Bam/",pattern="bam$")


#dir.create("Output/CSV/")

#select 1, 4, 14, 9, 10,
#s2<-c(20,32)
#bamfiles<-bamfiles[s2]

for (i in 2:length(bamfiles)){
        bam<-bamfiles[i]
        index<-paste0("~/programs/PID-master/Bam/",bam,'.bai')
        bf<-BamFile(paste0("~/programs/PID-master/Bam/",bam), index=index)

        file.name<-paste(bam)
        file.name<-substr(file.name,start=1,stop=7 )
        p_param <- PileupParam(max_depth=50000,include_insertions=TRUE)
        result<-pileup(bf, pileupParam = p_param, distinguish_strands=FALSE,ignore_query_Ns=FALSE)
        summary<-pileupFreq(result)

        summary$TotalReads<-rowSums(summary[3:6])
        maxr<-max(summary$TotalReads)
        print(file.name)
        cat("The maximum number of read depth is ", maxr)
        cat("\n")
        write.csv(summary, file=paste0("Output/CSV_PID/",file.name,".csv",collapse=""))
        
  }
