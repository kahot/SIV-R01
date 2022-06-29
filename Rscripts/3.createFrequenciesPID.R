library(tidyverse)
source("Rscripts/BaseRscript.R")

#dir.create("Output/SeqData/")

###########################

SIVFiles<-list.files("Output/CSV_PID/",pattern="csv")
#SIVFiles<-list.files("Output/CSV_PID/",pattern="^Run0|^Run2|^Run3")
SIVFiles<-list.files("Output/CSV_PID/",pattern="Run8")

coding.start<-190
coding.end<-681
no<-data.frame("pos"=c(coding.start:coding.end))
for (i in 1:length(SIVFiles)){
        print(i)
        id<-substr(paste(SIVFiles[i]),start=1,stop=7)
        print(id)
        SeqData<-read.csv(paste0("Output/CSV_PID/",SIVFiles[i]), row.names = 1, stringsAsFactors = F)
        
        SeqData<-SeqData[,-1]
        colnames(SeqData)[1]<-"pos"
        colnames(SeqData)[8:9]<-c("deletion","insertion") #deletion here is more like N
        colnames(SeqData)[2:5]<-c("a","c","g","t")

        #determine the majority nucleotide base at each site
        SeqData$MajNt<-apply(SeqData[,2:5],1,function(x) c("a","c","g","t")[which.max(x)])
        
        #read the refrence sequence:
        SeqData<-merge(no,SeqData,by="pos",all.x=T)
        reference<-read.dna("Data/AY032751.fasta", format = "fasta",as.character=TRUE)
        ref.code<-reference[coding.start:coding.end]
        SeqData$ref<-ref.code
        
         
        SeqData$transition.maj<-NA
        SeqData$transition.ref<-NA
        for (j in 1:nrow(SeqData)) SeqData$transition.maj[j]<-transition(SeqData$MajNt[j])        
        for (j in 1:nrow(SeqData)) SeqData$transition.ref[j]<-transition(SeqData$ref[j])
        
        #rearrange the columns
        SeqData<-SeqData[,c("a","c","g","t","deletion","insertion","N","pos","TotalReads","MajNt","ref","transition.maj","transition.ref")]
     
        #determine Transition mutation freq of every site.
        for (k in 1:nrow(SeqData)){
                if (is.na(SeqData$MajNt[k])) {
                        SeqData$freq.Ts[k]<-NA #transition mutations
                        SeqData$freq.Ts.ref[k]<-NA

                        SeqData$freq.transv[k]<-NA #transversion mutations
                        SeqData$freq.transv.ref[k]<-NA
                        SeqData$freq.transv1[k]<-NA
                        SeqData$freq.transv2[k]<-NA
                        SeqData$freq.transv1.ref[k]<-NA
                        SeqData$freq.transv2.ref[k]<-NA

                        SeqData$freq.mutations.ref[k]<-NA #all mutations
                        SeqData$freq.mutations[k]<-NA
                        
                        }
                else {
                      MajNum <- SeqData [k,which(c("a","c","g","t")==SeqData$MajNt[k])]
                      MutNum1<- SeqData [k,which(c("a","c","g","t")==SeqData$transition.maj[k])]
                      WTNum <- SeqData [k,which(c("a","c","g","t")==SeqData$ref[k])]
                      MutNum2<- SeqData [k,which(c("a","c","g","t")==SeqData$transition.ref[k])]
                      
                      SeqData$freq.Ts[k]<-MutNum1/SeqData$TotalReads[k]
                      SeqData$freq.Ts.ref[k]<-MutNum2/SeqData$TotalReads[k]
                      
 
                      #Frequenceis for specific transversion mutations (1 & 2)
                      Tvs1Num<-SeqData[k,which(c("a","c","g","t")==(transv1(SeqData$MajNt[k])))]
                      Tvs2Num<-SeqData[k,which(c("a","c","g","t")==(transv2(SeqData$MajNt[k])))]
                      SeqData$freq.transv1[k]<-Tvs1Num/SeqData$TotalReads[k]
                      SeqData$freq.transv2[k]<-Tvs2Num/SeqData$TotalReads[k]
                      Tvs1rNum<-SeqData[k,which(c("a","c","g","t")==(transv1(SeqData$ref[k])))]
                      Tvs2rNum<-SeqData[k,which(c("a","c","g","t")==(transv2(SeqData$ref[k])))]
                      SeqData$freq.transv1.ref[k]<-Tvs1rNum/SeqData$TotalReads[k]
                      SeqData$freq.transv2.ref[k]<-Tvs2rNum/SeqData$TotalReads[k]
                      
                      #mutation frequencies of all transversion mutations
                      SeqData$freq.transv[k]<-SeqData$freq.transv1[k]+SeqData$freq.transv2[k]
                      SeqData$freq.transv.ref[k]<-SeqData$freq.transv1.ref[k]+ SeqData$freq.transv2.ref[k]
                      
                      #Frequencies of all SNPs (no indels)
                      AllMutNum<-SeqData$TotalReads[k]-MajNum
                      AllMutNum2<-SeqData$TotalReads[k]-WTNum
                      
                      SeqData$freq.mutations[k]<-AllMutNum/SeqData$TotalReads[k]
                      SeqData$freq.mutations.ref[k]<-AllMutNum2/SeqData$TotalReads[k]
                      }
        }
                
        write.csv(SeqData,paste0("Output/SeqData_PID/SeqData_",id,".con.csv"))

}


## Max and Average read depth ###

csvs<-list.files("Output/CSV_PID/")
depth<-data.frame(files=csvs)
for (i in 1:length(csvs)){
    id<-substr(paste(csvs[i]),start=1,stop=7)
    
    df<-read.csv(paste0("Output/CSV_PID/",csvs[i]), row.names = 1, stringsAsFactors = F)
    ave<-mean(df$TotalReads, na.rm=T)
    mx<-max(df$TotalReads, na.rm=T)
    nb<-substr(id, start=6,stop=7)
    nb<-gsub("_","",nb)
    depth$File.name[i]<-id
    depth$run[i]<-substr(id, start=4,stop=4)
    depth$number[i]<-nb
    depth$Average[i]<-round(ave, digits = 1)
    depth$Max[i]<-mx
}

depth<-depth[,-1]
samples<-read.csv("Data/SamplesNoduplicates.csv",stringsAsFactors = F)
samples<-samples[,c("File.name","Tissue")]
depth2<-merge(depth, samples, by="File.name")

depth2$Tissue2<-ifelse(depth2$Tissue=="Plasma"|depth2$Tissue=="Stock", "Plasma", "Tissue")
aggregate(depth2["Average"], by=list(depth2$run,depth2$Tissue2), FUN=mean)


write.csv(depth2, "Output/ReadDeapth_all.csv")
