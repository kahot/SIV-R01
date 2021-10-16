library(ggplot2)
library(reshape)
library(ggpubr)
library(ggthemes)
library(plotrix)
library(grid)
library(gridExtra)
library(colorspace)
library(cowplot)
source("Rscripts/baseRscript2.R")
cols<-qualitative_hcl(6, palette="Dark3")

#Run5_18 is the control file
df<-read.csv(paste0("Output/OverviewF/Run5_18_filtered.overview.csv"),stringsAsFactors=F, row.names = 1)
df2<-read.csv(paste0("Output/Overview/Run5_18_overview.csv"),stringsAsFactors=F, row.names = 1)

ave<-round(mean(df$freq.mutations.ref, na.rm=T)*100,4)
ggplot(df, aes(x=pos, y=freq.mutations.ref*100))+
    ylab("Mutation frequency %")+xlab("")+
    geom_point(size=0.7, color="blue")+theme_bw()+
    ggtitle("Control (raw data)")+theme(plot.title = element_text(size=12))+
    annotate(geom='text', x=590, y=3.1, label=paste0("Ave. ", ave, "%"),color ='gray20', size=3.5,hjust=0)
ggsave("Output/Control_MF.pdf", width = 7, height = 2)

df$pos[df$freq.mutations.ref>=0.02] #197
df2<-df[df$pos!=197,]
mean(df$freq.mutations.ref, na.rm=T)   #0.001657977
mean(df2$freq.mutations.ref, na.rm=T)  #0.00156201
median(df$freq.mutations.ref, na.rm = T) #0.001478378
median(df2$freq.mutations.ref, na.rm = T) #0.001477813

# Non-filtered data (not removing sites <1000 reads)
df2<-read.csv(paste0("Output/Overview/Run5_18_overview.csv"),stringsAsFactors=F, row.names = 1)

ggplot(df2, aes(x=pos, y=freq.mutations.ref))+
    ylab("Mutation frequency")+xlab("")+
    geom_point(size=0.7)+theme_bw()+
    ggtitle("Control")+theme(plot.title = element_text(size=12))
#ggsave("Output/Control_MF_nonfiltered.pdf", width = 7, height = 2)

mean(df2$freq.mutations.ref, na.rm=T)  #0.001600204
median(df2$freq.mutations.ref, na.rm = T) #0.001355751

#PID processed data
#should be the same 
df3<-read.csv(paste0("Output/Overview_PID/Run5_18_overview.csv"),stringsAsFactors=F, row.names = 1)
df4<-read.csv(paste0("Output/OverviewF_PID/Run5_18_filtered.overview.csv"),stringsAsFactors=F, row.names = 1)

#non filtered PID 
ave<-round(mean(df3$freq.mutations.ref, na.rm=T)*100,4)
ggplot(df3, aes(x=pos, y=freq.mutations*100))+
    ylab("Mutation frequency %")+xlab("")+ylim(0,1)+
    geom_point(size=0.7, alpha=0.7, color='Blue')+theme_bw()+
    ggtitle("Control (PID processed)")+theme(plot.title = element_text(size=12))+
    annotate(geom='text', x=600, y=0.95, label=paste0("Ave. ", ave, "%"),color ='gray30', size=3,hjust=0)
ggsave("Output/Control_PID.pdf", width = 7, height = 2)


ave<-round(mean(df4$freq.mutations.ref, na.rm=T)*100,4)
ggplot(df4, aes(x=pos, y=freq.mutations.ref*100))+
    ylab("Mutation frequency %")+xlab("")+
    geom_point(size=0.7, color="Blue")+theme_bw()+
    ggtitle("Control")+theme(plot.title = element_text(size=12))+
annotate(geom='text', x=590, y=0.92, label=paste0("Ave. ", ave, "%"),color ='gray30', size=3.7,hjust=0)
ggsave("Output/Control_PID.pdf", width = 7, height = 2)

df4<-df4[!is.na(df4$TotalReads),]
nrow(df4[df4$freq.mutations>0.0027,]) #only 2 points above 0.27% mutation freq.
