library(ggplot2)
library(reshape)
library(ggpubr)
library(ggthemes)
library(plotrix)
library(grid)
library(colorspace)
source("Rscripts/baseRscript.R")
source("Rscripts/label_scientific.R")

cols2<-qualitative_hcl(6, palette="Dark3")


###############################       

samps<-read.csv("Data/SamplesNoduplicates.csv",stringsAsFactors = F)
samp<-samps[samps$Tissue=="Plasma",]

#siv.vl<-read.csv("Data/ViralLoads.csv", row.names = 1)
siv.reads<-read.csv("Output/ReadDeapth_all.csv", row.names = 1)
siv.reads<-siv.reads[siv.reads$Tissue=="Plasma",]

#siv.vl$id<-paste0(siv.vl$Monkey,".",siv.vl$Week)
#samp$id<-paste0(gsub("A",'',samp$Monkey),".",samp$Week)
#siv.vl<-merge(siv.vl, samp[,c("File.name", "id")], by="id")

#siv<-merge(siv.vl, siv.reads[c("File.name","Average")], by="File.name")

colnames(siv.reads)[4]<-"depth.PID"
siv<-siv.reads[,1:4]
#missing RUn2_10/Run2_5_ remove the file
siv<-siv[siv$File.name!="Run2_10",]
siv<-siv[siv$File.name!="Run2_5_",]
#remove Run5/7 files
siv<-siv[grep("Run0|Run2|Run3|Run4", siv$File.name),]
nonPIDfnames<-paste0(siv$File.name,"_overview.csv")
nonPID<-data.frame(File.name=siv$File.name)



#read depth without primerID (47 samples total)
for (i in 1:length(nonPIDfnames)){ 
    df<-read.csv(paste0("Output/Overview/",nonPIDfnames[i]), row.names = 1)
    m<-mean(df$TotalReads, na.rm=T)
    nonPID$depth.nonPID[i]<-m
}

siv<-merge(siv, nonPID, by="File.name")

p<-mean(siv$depth.PID)
np<-mean(siv$depth.nonPID)

p/np*100 #11.14656%
#starting RNA  10,000 RNA copies 

#remove Run3 as PID read depths exceeds 10,000 ????
siv2<-siv[siv$run!=3,]

## of unique sequence recovered 
p2<-mean(siv2$depth.PID)
p2 # 2685.825

#length of env sequenced ~500
#HCV sequenced was about !8500 bp
#8500/500 =17
#divide the sequence depth with 17 ro get genome wide coverage

mean(siv2$depth.nonPID)/17/10000*100
# 40.4% 
mean(siv2$depth.PID)/17/10000*100
#1.58%

#factoring into inefficiency of non-amplicon sequencing
mean(siv2$depth.PID)/50/10000*100
#0.537165 %
mean(siv2$depth.PID)/100/10000*100
#0.2685825


#average viral load of HCV 1A (per 500ul used)
tmplte<-20201990/2

0.016*tmplte
#161615
0.00537165*tmplte
#54259.01
