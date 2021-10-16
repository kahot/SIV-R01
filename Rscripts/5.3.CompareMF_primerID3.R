library(reshape)
library(ggpubr)
library(ggthemes)
library(plotrix)
library(grid)
library(colorspace)
source("Rscripts/baseRscript.R")
source("Rscripts/label_scientific.R")

cols2<-qualitative_hcl(6, palette="Dark3")


# read the files saved in Overview (unfiltered):

#File processed without primerID
ov1<-read.csv("Output/Overview/Run3_9__overview.csv", row.names = 1, stringsAsFactors = F)
colnames(ov1)[1]<-"posF"
ov1$pos<-ov1$posF-54

#dir.create("Output/MF_PID/Comparison")

#PID-Master processed
ov2<-read.csv("Output/Overview_PIDcon/Run3_9__overview.csv", row.names = 1, stringsAsFactors = F)

mf1<-ov1[,c("pos","freq.Ts")]
mf2<-ov2[,c("pos","freq.Ts")]
colnames(mf2)[2]<-"PID_processed"
Ts<-merge(mf1,mf2,by="pos")
colnames(Ts)[2]<-"noPrimerID"
write.csv(Ts, "Output/MF_PID/Comparison/Run3_9_Ts.new_comparison.csv")

Ts.m<-melt(Ts, id.vars = "pos")
    
ggplot(data=Ts.m, aes(x=pos, y=value, color=variable))+
        geom_point(size=0.7)+theme_bw()+ylab("Mutation frequency")+
        xlab("ENV position")+
        scale_y_continuous(trans = 'log10', labels=label_scientific)+
            theme(legend.title=element_blank())
ggsave("Output/MF_PID/Comparison/MFplot_log.pdf", width = 8, height = 4)
    
    
ggplot(data=Ts.m, aes(x=pos, y=value, color=variable))+
        geom_point(size=0.7)+theme_bw()+ theme(legend.title=element_blank())
ggsave("Output/MF_PID/Comparison/MFplot.pdf", width = 8, height = 4)
    
  
    
    
#plot the difference
Ts$diff<-Ts$noPrimerID-Ts$PID_processed
ggplot(data=Ts, aes(x=pos, y=abs(diff)))+
    geom_point(size=0.5, color="blue")+theme_bw()+ylab("Difference in Mutation frequency")+
    xlab("ENV position")+
    scale_y_continuous(trans = 'log10', labels=label_scientific)+
    theme(legend.title=element_blank())
ggsave("Output/MF_PID/Comparison/MF_diff.log.noIDvsPID.pdf", width = 7.5, height = 4)    
mean(abs(Ts$diff), na.rm=T) #0.01858248

ggplot(data=Ts, aes(x=pos, y=diff))+
    geom_point(size=0.5, color="blue")+theme_bw()+ylab("Difference in Mutation frequency")+
    xlab("ENV position")+
    theme(legend.title=element_blank())
ggsave("Output/MF_PID/Comparison/MF_diff.noIDvsPID.pdf", width = 7.5, height = 4)    



#### 



# # of sites freq(No primerID) > freq(primerID)
nrow(Ts[Ts$diff>0,])# 347
# of sites freq(No primerID) < freq(primerID)
nrow(Ts[Ts$diff<0,]) #91

OV<-list()
OV[[1]]<-ov1
OV[[2]]<-ov2

summaryP<-data.frame(Method=c("noPrimerID", "PID"))
summaryP$ave.depth<-unlist(lapply(OV, function(x) mean(x$TotalReads, na.rm=T)))
summaryP$ave.Ts.freq <-unlist(lapply(OV, function(x) mean(x$freq.Ts, na.rm=T)))
summaryP$ave.Tvs.freq<-unlist(lapply(OV, function(x) mean(x$freq.transv, na.rm=T)))
summaryP$ave.All.freq<-unlist(lapply(OV, function(x) mean(x$freq.mutations, na.rm=T)))

#      Method ave.depth ave.Ts.freq ave.Tvs.freq ave.All.freq
#1 noPrimerID  65378.33 0.010836888  0.003062286   0.01389917
#2        PID  21020.60 0.009896337  0.003565465   0.01346180

summaryP$Method<-factor(summaryP$Method, levels=c("noPrimerID", "PID"))

sumPm<-melt(summaryP, id.vars = "Method")
sumPm<-sumPm[-c(1:3),]
ggplot(data=summaryP, aes(x=Method, y=ave.depth))+
    ylab('Average depth')+
    geom_bar( stat="identity", color="lightblue", fill=cols[1])+theme_bw()
#ggsave("Output/primID_comparison/Compare_depth.pdf", height = 3, width = 4)

percent.increase2<-apply(summaryP[,3:5], 2, function(x) x[1]/x[2])
#ave.Ts.freq ave.Tvs.freq ave.All.freq 
#1.0950402    0.8588742    1.0324898 
