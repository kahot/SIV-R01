library(Rsamtools)
library(stringr)

source("Rscripts/pileupFreq.R")

#number of sampels to process
bamfiles<-list.files("Output/bam/",pattern="bam$")
dir.create("Output/CSV/")

for (i in 1:length(bamfiles)){
        bam<-bamfiles[i]
        index<-paste0("Output/bam/",bam,'.bai')
        bf<-BamFile(paste0("Output/bam/",bam), index=index)

        file.name<-paste(bam)
        file.name<-substr(file.name,start=1,stop=7 )
        p_param <- PileupParam(max_depth=70000)
        result<-pileup(bf, pileupParam = p_param, distinguish_strands=FALSE,ignore_query_Ns=FALSE)
        summary<-pileupFreq(result)

        summary$TotalReads<-rowSums(summary[3:6])
        maxr<-max(summary$TotalReads)
        print(file.name)
        cat("The maximum number of read depth is ", maxr)
        cat("\n")
        write.csv(summary, file=paste0("Output/CSV/",file.name,".csv",collapse=""))
        #assign(file.name, summary)
        
  }
