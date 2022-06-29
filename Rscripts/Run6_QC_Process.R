
#### Check mapped raw files (processed at Q15) 

library(Rsamtools)
library(stringr)
source("Rscripts/baseRscript.R")
source("Rscripts/pileupFreq.R")

bamfiles<-list.files("Output/bam/",pattern="^Run6_.*15_sort.bam$")
bamfiles<-list.files("Output/bam/",pattern="^Run6_.*15mapped.bam$")

for (i in 1:length(bamfiles)){
    bam<-bamfiles[i]
    index<-paste0("Output/bam/",bam,'.bai')
    bf<-BamFile(paste0("Output/bam/",bam), index=index)
    
    file.name<-paste(bam)
    file.name<-substr(file.name,start=1,stop=7 )
    p_param <- PileupParam(max_depth=300000,include_insertions=TRUE,include_deletions=TRUE)
    result<-pileup(bf, pileupParam = p_param, distinguish_strands=FALSE,ignore_query_Ns=FALSE)
    summary<-pileupFreq(result)
    
    summary$TotalReads<-rowSums(summary[3:6])
    maxr<-max(summary$TotalReads)
    print(file.name)
    cat("The maximum number of read depth is ", maxr)
    cat("\n")
    write.csv(summary, file=paste0("Output/QC/",file.name,".QC.csv",collapse=""))
}

#Run6_6_, Run6_23, Run6_18, Run6_19, 

###########################
library(tidyverse)
SIVFiles<-list.files("Output/QC/",pattern=".QC.csv")

coding.start<-190
coding.end<-681
no<-data.frame("pos"=c(coding.start:coding.end))
Seq<-list()
for (i in 1:length(SIVFiles)){
    print(i)
    id<-substr(paste(SIVFiles[i]),start=1,stop=7)
    print(id)
    #SeqData<-read.csv(paste0("Output/CSV_PID/",SIVFiles[i]), row.names = 1, stringsAsFactors = F)
    SeqData<-read.csv(paste0("Output/QC/",SIVFiles[i]), row.names = 1, stringsAsFactors = F)
    
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
            
            #mutation frequencies of all transversion mutataions
            SeqData$freq.transv[k]<-SeqData$freq.transv1[k]+SeqData$freq.transv2[k]
            SeqData$freq.transv.ref[k]<-SeqData$freq.transv1.ref[k]+ SeqData$freq.transv2.ref[k]

            #Frequencies of all SNPs (no indels)
            AllMutNum<-SeqData$TotalReads[k]-MajNum
            AllMutNum2<-SeqData$TotalReads[k]-WTNum
            
            SeqData$freq.mutations[k]<-AllMutNum/SeqData$TotalReads[k]
            SeqData$freq.mutations.ref[k]<-AllMutNum2/SeqData$TotalReads[k]
        }
    }
    Seq[[i]]<-SeqData
    names(Seq)[i]<-id
}
############################### 
#Create overview files
library(dplyr)
Overview<-list()
for (i in 1:length(Seq)){   
    id<-names(Seq)[i]
    print(id)
    OverviewDF<-Seq[[i]]
    
    TypeOfSite<-c() 
    TypeOfSite.tv1<-c()
    TypeOfSite.tv2<-c()
    TypeOfSite2<-c()
    TypeOfSite2.tv1<-c()
    TypeOfSite2.tv2<-c()
    
    for (codon in 1:(nrow(OverviewDF)/3)) { #modify based on reading frame
        positions <- c(codon*3-2,codon*3-1, codon*3)  #starting with codon 1
        WTcodon <- OverviewDF$MajNt[positions]  
        Refcodon<-OverviewDF$ref[positions] 
        if (is.na(WTcodon[1])|is.na(WTcodon[2])|is.na(WTcodon[3])){ 
            WTcodon<-c('n','n','n')
            mutant1codon<-c('n','n','n')
            mutant2codon<-c('n','n','n')
            mutant3codon<-c('n','n','n')
            
            mutant1codon.tv1 <- c('n','n','n')
            mutant2codon.tv1 <- c('n','n','n')
            mutant3codon.tv1 <-c('n','n','n')
            
            mutant1codon.tv2 <- c('n','n','n')
            mutant2codon.tv2 <- c('n','n','n')
            mutant3codon.tv2 <- c('n','n','n')
        }
        
        else {                        
            mutant1codon <- c(transition(WTcodon[1]), WTcodon[2:3])  
            mutant2codon <- c(WTcodon[1],transition(WTcodon[2]), WTcodon[3])
            mutant3codon <- c(WTcodon[1:2], transition(WTcodon[3]))
            
            #transversion mutation to 'a' or 'c'
            mutant1codon.tv1 <- c(transv1(WTcodon[1]), WTcodon[2:3]) 
            mutant2codon.tv1 <- c(WTcodon[1],transv1(WTcodon[2]), WTcodon[3])
            mutant3codon.tv1 <- c(WTcodon[1:2], transv1(WTcodon[3]))
            #transversion mutation to 'g' or 't'
            mutant1codon.tv2 <- c(transv2(WTcodon[1]), WTcodon[2:3])  
            mutant2codon.tv2 <- c(WTcodon[1],transv2(WTcodon[2]), WTcodon[3])
            mutant3codon.tv2 <- c(WTcodon[1:2], transv2(WTcodon[3]))
        }
        #compare to the ref seq
        mutant1codon2 <- c(transition(Refcodon[1]), Refcodon[2:3])  #If the first position has transistion mutation, it's labeld as mutatnt1codon.
        mutant2codon2 <- c(Refcodon[1],transition(Refcodon[2]), Refcodon[3])
        mutant3codon2 <- c(Refcodon[1:2], transition(Refcodon[3]))
        
        #transversion mutation to 'a' or 'c'
        mutant1codon2.tv1 <- c(transv1(Refcodon[1]), Refcodon[2:3]) 
        mutant2codon2.tv1 <- c(Refcodon[1],transv1(Refcodon[2]), Refcodon[3])
        mutant3codon2.tv1 <- c(Refcodon[1:2], transv1(Refcodon[3]))
        #transversion mutation to 'g' or 't'
        mutant1codon2.tv2 <- c(transv2(Refcodon[1]), Refcodon[2:3])  
        mutant2codon2.tv2 <- c(Refcodon[1],transv2(Refcodon[2]), Refcodon[3])
        mutant3codon2.tv2 <- c(Refcodon[1:2], transv2(Refcodon[3]))
        
        
        
        TypeOfSite<-c(TypeOfSite,typeofsitefunction(WTcodon,mutant1codon))
        TypeOfSite<-c(TypeOfSite,typeofsitefunction(WTcodon,mutant2codon))
        TypeOfSite<-c(TypeOfSite,typeofsitefunction(WTcodon,mutant3codon))
        
        TypeOfSite.tv1<-c(TypeOfSite.tv1,typeofsitefunction(WTcodon,mutant1codon.tv1))
        TypeOfSite.tv1<-c(TypeOfSite.tv1,typeofsitefunction(WTcodon,mutant2codon.tv1))
        TypeOfSite.tv1<-c(TypeOfSite.tv1,typeofsitefunction(WTcodon,mutant3codon.tv1))
        
        TypeOfSite.tv2<-c(TypeOfSite.tv2,typeofsitefunction(WTcodon,mutant1codon.tv2))
        TypeOfSite.tv2<-c(TypeOfSite.tv2,typeofsitefunction(WTcodon,mutant2codon.tv2))
        TypeOfSite.tv2<-c(TypeOfSite.tv2,typeofsitefunction(WTcodon,mutant3codon.tv2))
        
        
        TypeOfSite2<-c(TypeOfSite2,typeofsitefunction(Refcodon,mutant1codon2))
        TypeOfSite2<-c(TypeOfSite2,typeofsitefunction(Refcodon,mutant2codon2))
        TypeOfSite2<-c(TypeOfSite2,typeofsitefunction(Refcodon,mutant3codon2))
        
        TypeOfSite2.tv1<-c(TypeOfSite2.tv1,typeofsitefunction(Refcodon,mutant1codon2.tv1))
        TypeOfSite2.tv1<-c(TypeOfSite2.tv1,typeofsitefunction(Refcodon,mutant2codon2.tv1))
        TypeOfSite2.tv1<-c(TypeOfSite2.tv1,typeofsitefunction(Refcodon,mutant3codon2.tv1))
        
        TypeOfSite2.tv2<-c(TypeOfSite2.tv2,typeofsitefunction(Refcodon,mutant1codon2.tv2))
        TypeOfSite2.tv2<-c(TypeOfSite2.tv2,typeofsitefunction(Refcodon,mutant2codon2.tv2))
        TypeOfSite2.tv2<-c(TypeOfSite2.tv2,typeofsitefunction(Refcodon,mutant3codon2.tv2))     
    } 
    
    OverviewDF$Type<-TypeOfSite[1:length(OverviewDF$pos)]
    OverviewDF$Type.tv1<-TypeOfSite.tv1[1:length(OverviewDF$pos)]
    OverviewDF$Type.tv2<-TypeOfSite.tv2[1:length(OverviewDF$pos)]
    
    OverviewDF$Type.r<-TypeOfSite2[1:length(OverviewDF$pos)]
    OverviewDF$Type.tv1.r<-TypeOfSite2.tv1[1:length(OverviewDF$pos)]
    OverviewDF$Type.tv2.r<-TypeOfSite2.tv2[1:length(OverviewDF$pos)]
    
    Overview[[i]]<-OverviewDF[,-c(1:7)]
    
    names(Overview)[i]<-id   
}
###############################
Overview_sum<-list()
for (i in 1:length(Overview)){
    
    OverviewDF<-Overview[[i]]
    id<-names(Overview)[i]
    
    OverviewDF$WTAA<-""
    OverviewDF$MUTAA<-""
    OverviewDF$TVS1_AA<-""
    OverviewDF$TVS2_AA<-""
    OverviewDF$TVS2_AA<-""
    for (k in 1:nrow(OverviewDF)){
        
        if (k%%3==1){
            if (is.na(OverviewDF$MajNt[k])|is.na(OverviewDF$MajNt[k+1])|is.na(OverviewDF$MajNt[k+2])) { OverviewDF$WTAA[k]<-"NA"
            OverviewDF$MUTAA[k]<-"NA"
            OverviewDF$TVS1_AA[k]<-"NA"
            OverviewDF$TVS2_AA[k]<-"NA"}
            else {  OverviewDF$WTAA[k] = seqinr::translate(OverviewDF$MajNt[c(k,k+1,k+2)])
            OverviewDF$MUTAA[k] = seqinr::translate(c(transition(OverviewDF$MajNt[k]),OverviewDF$MajNt[c(k+1,k+2)]))
            OverviewDF$TVS1_AA[k] = seqinr::translate(c(transv1(OverviewDF$MajNt[k]),OverviewDF$MajNt[c(k+1,k+2)]))
            OverviewDF$TVS2_AA[k] = seqinr::translate(c(transv2(OverviewDF$MajNt[k]),OverviewDF$MajNt[c(k+1,k+2)]))}
        } 
        if (k%%3==2){
            if (is.na(OverviewDF$MajNt[k-1])|is.na(OverviewDF$MajNt[k])|is.na(OverviewDF$MajNt[k+1]))  {OverviewDF$WTAA[k]<-"NA"
            OverviewDF$MUTAA[k]<-"NA"
            OverviewDF$TVS1_AA[k]<-"NA"
            OverviewDF$TVS2_AA[k]<-"NA"}
            else {  OverviewDF$WTAA[k] = seqinr::translate(OverviewDF$MajNt[c(k-1,k,k+1)])
            OverviewDF$MUTAA[k] = seqinr::translate(c(OverviewDF$MajNt[c(k-1)],transition(OverviewDF$MajNt[k]),OverviewDF$MajNt[c(k+1)]))
            OverviewDF$TVS1_AA[k] = seqinr::translate(c(OverviewDF$MajNt[c(k-1)],transv1(OverviewDF$MajNt[k]),OverviewDF$MajNt[c(k+1)]))
            OverviewDF$TVS2_AA[k] = seqinr::translate(c(OverviewDF$MajNt[c(k-1)],transv2(OverviewDF$MajNt[k]),OverviewDF$MajNt[c(k+1)]))}
        }
        if (k%%3==0){
            if (is.na(OverviewDF$MajNt[k-2])|is.na(OverviewDF$MajNt[k-1])|is.na(OverviewDF$MajNt[k]))  {  OverviewDF$WTAA[k]<-"NA"
            OverviewDF$MUTAA[k]<-"NA"
            OverviewDF$TVS1_AA[k]<-"NA"
            OverviewDF$TVS2_AA[k]<-"NA"}
            else {  OverviewDF$WTAA[k] = seqinr::translate(OverviewDF$MajNt[c(k-2,k-1,k)])
            OverviewDF$MUTAA[k] = seqinr::translate(c(OverviewDF$MajNt[c(k-2,k-1)],transition(OverviewDF$MajNt[k])))
            OverviewDF$TVS1_AA[k] = seqinr::translate(c(OverviewDF$MajNt[c(k-2,k-1)],transv1(OverviewDF$MajNt[k])))
            OverviewDF$TVS2_AA[k] = seqinr::translate(c(OverviewDF$MajNt[c(k-2,k-1)],transv2(OverviewDF$MajNt[k])))}
        }
        
    }
    #Add whether AA change is drastic & makes CpG
    OverviewDF$bigAAChange<-0
    OverviewDF$bigAAChange.tv1<-0
    OverviewDF$bigAAChange.tv2<-0
    OverviewDF$makesCpG <- 0
    OverviewDF$makesCpG.tvs <- 0
    OverviewDF$makesCpG.tv1 <- 0
    OverviewDF$makesCpG.tv2 <- 0
    
    for(j in 1:nrow(OverviewDF)){
        WT <- amCat(OverviewDF[j,'WTAA'])
        MUT <- amCat(OverviewDF[j,'MUTAA'])
        MUT1<-amCat(OverviewDF[j,'TVS1_AA'])
        MUT2<-amCat(OverviewDF[j,'TVS2_AA'])
        
        if (WT != MUT) OverviewDF$bigAAChange[j] <- 1
        if (WT != MUT1) OverviewDF$bigAAChange.tv1[j] <- 1
        if (WT != MUT2) OverviewDF$bigAAChange.tv2[j] <- 1
        
        trip <- OverviewDF$MajNt[c(j, j+1,j+2)]
        if (is.na(trip[1])|is.na(trip[2])|is.na(trip[3])) 
            next
        else{
            if (trip[1] == "c" & trip[2] == "a" ) OverviewDF$makesCpG[j] <- 1 
            if (trip[2] == "t" & trip[3] == "g")  OverviewDF$makesCpG[j] <- 1
            if (trip[1] == "c" & (trip[2]=="c"|trip[2]=='t')) OverviewDF$makesCpG.tvs[j] <- 1
            if (trip[3] == "g" & (trip[2]=="a"|trip[2]=="g")) OverviewDF$makesCpG.tvs[j] <- 1
            
            if (trip[1] == "c" & (trip[2]=="c"|trip[2]=='t')) OverviewDF$makesCpG.tv2[j] <- 1                                
            if (trip[3] == "g" & (trip[2]=="a"|trip[2]=="g")) OverviewDF$makesCpG.tv1[j] <- 1
            
        }
    } 
    
    
    write.csv(OverviewDF,paste0("Output/QC/",id,"_QCoverview.csv"))
    
    #filter the sites with reads<10000
    remove<-which(OverviewDF$TotalReads<1000)
    OverviewDF[remove, 7:16]<-NA
    write.csv(OverviewDF,paste0("Output/QC/",id,"_QCfiltered.overview.con.csv"))
}

#Pull out the depth and freq at the six sites of intesest
### Sites of interest
diffpos<-c(475,494,500,515,518,578)
aapos<-ceiling(diffpos/3)
summary<-read.csv("Output/MF_PID/HighFreq/Summary_highFreq_sites.csv", row.names = 1)
summary2<-summary[summary$pos %in% diffpos,]
for (i in 1:nrow(summary2)){
    if (summary2$Type.1[i]=="Ts") summary2$mut[i]<-transition(summary2$ref[i])
    if (summary2$Type.1[i]=="Tv1") summary2$mut[i]<-transv1 (summary2$ref[i])
    if (summary2$Type.1[i]=="Tv2") summary2$mut[i]<-transv2 (summary2$ref[i])
}
summary2$mut[summary2$Type.1=="Ts"]<-transition(summary2$ref)
MutNames<-paste0(summary2$ref, summary2$pos,summary2$mut)

files<-list.files("Output/QC/", pattern="_QCoverview.csv")
DF<-list()
for (i in 1: length(files)){
    df<-read.csv(paste0("Output/QC/",files[i]), stringsAsFactors = F, row.names = 1)
    DF[[i]]<-df[df$pos %in% diffpos,]
    names(DF)[i]<-substr(files[i], 1,7)
}

#Find the read depth and freq at the six sites from Q15 raw data
rawMF<-data.frame()
rawDepth<-data.frame()
for (i in 1:nrow(summary2)){
    muttype<-summary2$Type.1[i]
    if (muttype=="Tv1") muttype="transv1"
    if (muttype=="Tv2") muttype="transv2"
    
    df<-data.frame(File.name=names(DF))
    df[,"MF"]<-sapply(DF, function(x) x[i,paste0("freq.",muttype,".ref")])
    df$Mutation<-MutNames[i]
    rawMF<-rbind(rawMF, df)
    df2<-data.frame(File.name=names(DF))
    df2[,"Depth"]<-sapply(DF, function(x) x[i,"TotalReads"])
    df2$Mutation<-MutNames[i]
    rawDepth<-rbind(rawDepth, df2)

}

rawMF$Method<-"Raw Q15"
rawDepth$Method<-"Raw Q15"


####################################
### extract mf and depth for the six sites from PID processed data
IVfiles<-list.files("Output/Overview_PIDcon/", pattern="Run6")

Ov<-list()
for (i in 1:length(SIVfiles)){ 
    df<-read.csv(paste0("Output/Overview_PIDcon/",SIVfiles[i]),stringsAsFactors=F, row.names = 1)
    Ov[[i]]<-df[df$pos %in% diffpos,]
    names(Ov)[i]<-substr(paste(SIVfiles[i]),start=1,stop=7)
}

#Find the read depth and freq at the six sites from Q15 raw data

PIDmf<-data.frame()
PIDdepth<-data.frame()
for (i in 1:nrow(summary2)){
    muttype<-summary2$Type.1[i]
    if (muttype=="Tv1") muttype="transv1"
    if (muttype=="Tv2") muttype="transv2"
    
    df<-data.frame(File.name=names(Ov))
    df[,"MF"]<-sapply(Ov, function(x) x[i,paste0("freq.",muttype,".ref")])
    df$Mutation<-MutNames[i]
    PIDmf<-rbind(PIDmf, df)
    df2<-data.frame(File.name=names(Ov))
    df2[,"Depth"]<-sapply(Ov, function(x) x[i,"TotalReads"])
    df2$Mutation<-MutNames[i]
    PIDdepth<-rbind(PIDdepth, df2)
    
}
PIDmf$Method<-"PID"
PIDdepth$Method<-"PID"


MF<-rbind(rawMF,PIDmf)
Depth<-rbind(rawDepth,PIDdepth)

ggplot()+
    ylab("Mutation frequency")+xlab("Sample")+
    facet_wrap(~ Mutation, ncol=3)+
    geom_point(data=MF, aes(x=File.name, y=MF, color=Method),alpha=0.8,position=position_dodge(width = 0.5),size=2)+
    theme_bw()+
    scale_color_manual(values=paste0(cols2[c(7,4)],"CC"))+
    theme(axis.text.x = element_text(angle=90),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("Output/QC/MF_comparison_Run6.pdf", width = 10, height = 6 )

ggplot()+
    ylab("log(Read depth)")+xlab("Sample")+
    facet_wrap(~ Mutation, ncol=3)+
    geom_bar(data=Depth, aes(x=File.name,y=log10(Depth), fill=Method), alpha=0.7, position=position_dodge(width = 0.5), stat="identity")+
    theme_bw()+
    scale_fill_manual(values=paste0(cols2[c(7,4)]))+
    theme(axis.text.x = element_text(angle=90),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("Output/QC/Depth_comparison_Run6.pdf", width = 10, height = 6 )


write.csv(rawMF, "Output/QC/RawMF.6sites.csv")
write.csv(rawDepth, "Output/QC/RawDepth.6sites.csv")


MF2<-cbind(MF, Depth[,2])
colnames(MF2)[5]<-"Depth"
ggplot(MF2, aes(x=File.name))+
    facet_wrap(~ Mutation, ncol=6)+
    geom_point(aes(y=MF*4, color=Method),alpha=0.8,position=position_dodge(width = 0.5),size=2)+
    geom_bar(aes(y=log10(Depth), fill=Method), alpha=0.3, position=position_dodge(width = 0.5), stat="identity")+
    geom_point(aes(y=MF*4, color=Method),alpha=0.8,position=position_dodge(width = 0.5),size=2)+
    scale_y_continuous(
            name="Mutation frequency", breaks = c(0,2,4),labels = c(0,0.5,1),
            sec.axis=sec_axis(~., name="log(Read depth)"))+
    xlab("Sample")+
    theme_bw()+
    scale_color_manual(values=cols2[c(7,4)])+
    scale_fill_manual(values=cols2[c(7,4)])+
    theme(axis.text.x = element_text(angle=90),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("Output/QC/MF_Depth_comparison_Run6.pdf", width = 25, height = 6 )


write.csv(MF2,"Output/QC/MF_Depth_comparison_Run6.csv")

#####
rawMF<-data.frame()
for (i in 1:nrow(summary2)){
    df<-data.frame(File.name=names(DF))
    df[,"MF"]<-sapply(DF, function(x) x[i,paste0("freq.mutations.ref")])
    df$Mutation<-MutNames[i]
    rawMF<-rbind(rawMF, df)
    
}

ggplot()+
    ylab("Mutation frequency")+xlab("Sample")+
    facet_wrap(~ Mutation, ncol=3)+
    geom_point(data=rawMF, aes(x=File.name, y=MF),alpha=0.8,position=position_dodge(width = 0.5),size=2)+
    theme_bw()+
    scale_color_manual(values=paste0(cols2[c(7,4)],"CC"))+
    theme(axis.text.x = element_text(angle=90),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
