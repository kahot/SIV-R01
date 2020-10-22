library(ggplot2)
library(reshape)
library(ggpubr)
library(plotrix)
library(grid)
library(purrr)
library(zoo)
library(tidyverse)
library(gplots)
library(gridExtra)
library(DataCombine)
source("Rscripts/baseRscript2.R")

# read the files saved in Overview_output:
SIVFiles_overview<-list.files("Output/OverviewF/",pattern="overview.csv")
Seqfiles<-list.files("Output/SeqData/", pattern=".csv")
Overview_summary<-list()
for (i in 1:length(SIVFiles_overview)){ 
    ov<-read.csv(paste0("Output/OverviewF/",SIVFiles_overview[i]),stringsAsFactors=FALSE, row.names = 1)
    SeqDF<-read.csv(paste0("Output/SeqData/",Seqfiles[i]),stringsAsFactors=FALSE, row.names = 1)
    ov<-merge(SeqDF[,c("pos","a","c","g","t")], ov, by="pos")
    Overview_summary[[i]]<-ov
    names(Overview_summary)[i]<-substr(paste(SIVFiles_overview[i]),start=1,stop=7)
}
dir.create("Output/Overview_D/")
######
Overview.d<-list()
pi<-data.frame(SampleID=names(Overview_summary))
for (i in 1:length(Overview_summary)){
    dat<-Overview_summary[[i]]
    filename<-names(Overview_summary)[i]
    
    dat$D<-NA
    for (k in 1:nrow(dat)){
        if (is.na(dat$MajNt[k])|is.na(dat$freq.Ts[k]))  next 
        else{
            ac<-dat[k, "a"]*dat[k,"c"]
            ag<-dat[k, "a"]*dat[k,"g"]
            at<-dat[k, "a"]*dat[k,"t"]
            cg<-dat[k, "c"]*dat[k,"g"]
            ct<-dat[k, "c"]*dat[k,"t"]
            gt<-dat[k, "g"]*dat[k,"t"]
            m<-(dat[k,"TotalReads"]^2-dat[k,"TotalReads"])/2
            dat$D[k]<-(ac+ag+at+cg+ct+gt)/m
        }
    }
    Overview.d[[i]]<-dat
    names(Overview.d)[i]<-filename
    
    n<-nrow(dat[!is.na(dat$D),])
    pi$SampleID[i]<-filename
    pi$Pi[i]<-(sum(dat$D, na.rm=T))/n
    write.csv(dat,paste0("Output/Overview_D/",filename,"_overviewD.csv"))
}


## Plot the results:   
## Sample ID information ###
samples2<-read.csv("Data/SamplesNoduplicates.csv",stringsAsFactors = F)
ids<-c("A21918","A22317","A22517", "A22617","A22217","A23918","A22117")
list.animal<-split(samples2, samples2$Monkey)
monkeyList<-list()
k=1
for (i in 1:length(list.animal)){
    if (nrow(list.animal[[i]])>1){
        monkeyList[[k]]<-list.animal[[i]]
        names(monkeyList)[k]<-names(list.animal)[i]
        k=k+1
    }
}
monkeys<-names(monkeyList)
monkeys2<-monkeyList[ids]

tb<-c(17.5,20,20,20,17, 17.5, 19)
tbs<-data.frame(cbind(ids,tb))
tbs$tb<-as.numeric(as.character(tbs$tb))


stock<-pi$Pi[pi$SampleID=="Run0_17"]

Plot<-list()
Pi<-list()
for (i in 1:length(monkeys2)){
    print(ids[i])
    sample<-monkeys2[[i]]
    monkey<-names(monkeys2)[i]
    nucd<-pi[pi$SampleID %in%as.vector(sample$File.name), ]
    nucd<-merge(nucd, sample[,c("File.name","Week")], by.x="SampleID", by.y="File.name")
    nucd$SampleID<-as.character(nucd$SampleID)
    nucd = nucd[order(nucd[,'Week']),]
    nucd<-InsertRow(nucd, c("Stock",stock,0), RowNum=1)
    nucd$Week<-as.integer(nucd$Week)
    nucd$Pi<-as.numeric(nucd$Pi)
    Pi[[i]]<-nucd
    names(Pi)[i]<-monkey
    
    tbweek<-tbs$tb[tbs$ids==monkey]
    
    Plot[[i]]<-ggplot(data=nucd, aes(x=Week, y=Pi))+
        geom_point(color="blue")+ylab("Nucleotide diversity")+
        theme_bw()+geom_vline(xintercept=tbweek, col="deeppink")+ylim(0,0.043)+
        ggtitle(paste0(monkey))+theme(plot.title = element_text(size=12))
    
    
}


pdf("Output/NucDiversity/Nuc.diversity_overweek.pdf", width = 8, height = 7)
do.call(grid.arrange, c(Plot, ncol=2))
dev.off()

