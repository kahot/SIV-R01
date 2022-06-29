library(reshape)
library(ggpubr)
library(ggthemes)
library(plotrix)
library(grid)
library(gridExtra)
library(colorspace)
library(cowplot)
source("Rscripts/baseRscript.R")
cols<-qualitative_hcl(6, palette="Dark3")

#Run5_18 is the control file

#non-PID processed file
df<- read.csv("Output/OverviewF/Run5_18_filtered.overview.csv",stringsAsFactors=F, row.names = 1)
#df2<-read.csv("Output/Overview/Run5_18_overview.csv",stringsAsFactors=F, row.names = 1)

ave<-round(mean(df$freq.mutations.ref, na.rm=T)*100,4)
ggplot(df, aes(x=pos, y=freq.mutations.ref*100))+
    ylab("Mutation frequency %")+xlab("")+
    geom_point(size=0.7, color="blue")+theme_bw()+
    ggtitle("Control (raw data)")+theme(plot.title = element_text(size=12))+
    annotate(geom='text', x=590, y=3.1, label=paste0("Ave. ", ave, "%"),color ='gray20', size=3.5,hjust=0)
ggsave("Output/Control_MF.pdf", width = 7, height = 2)

#which one has the high freq?
df$pos[df$freq.mutations.ref>=0.02] #position 197

#Removing pos 197?
df2<-df[df$pos!=197,]
mean(df$freq.mutations.ref, na.rm=T)   #0.001657977
mean(df2$freq.mutations.ref, na.rm=T)  #0.00156201
median(df$freq.mutations.ref, na.rm = T) #0.001478378
median(df2$freq.mutations.ref, na.rm = T) #0.001477813
#Does not change the results much


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
#df4<-read.csv(paste0("Output/OverviewF_PID/Run5_18_filtered.overview.con.csv"),stringsAsFactors=F, row.names = 1)
#no reads above 1000 for Run5_18

#Filter at reads >100
df4<-df3
df4[which(df4$TotalReads<100), 7:16]<-NA
    
#non filtered PID 
ave<-round(mean(df3$freq.mutations.ref, na.rm=T)*100,4)
ggplot(df3, aes(x=pos, y=freq.mutations*100))+
    ylab("Mutation frequency %")+xlab("")+ylim(0,1)+
    geom_point(size=0.7, alpha=0.7, color='Blue')+theme_bw()+
    ggtitle("Control (PID processed)")+theme(plot.title = element_text(size=12))+
    annotate(geom='text', x=600, y=0.95, label=paste0("Ave. ", ave, "%"),color ='gray30', size=3,hjust=0)
ggsave("Output/Control_PID_nonFiltered.pdf", width = 7, height = 2)


ave<-round(mean(df4$freq.mutations.ref, na.rm=T)*100,4)
ggplot(df4, aes(x=pos, y=freq.mutations.ref*100))+
    ylab("Mutation frequency %")+xlab("")+ylim(0,1)+
    geom_point(size=0.7, color="Blue")+theme_bw()+
    ggtitle("Control")+theme(plot.title = element_text(size=12))+
annotate(geom='text', x=590, y=0.92, label=paste0("Ave. ", ave, "%"),color ='gray30', size=3.7,hjust=0)
ggsave("Output/Control_PID_filtered100.pdf", width = 7, height = 2)

df4<-df4[!is.na(df4$freq.Ts),]
nrow(df4[df4$freq.mutations>0.0027,]) #only 2 points above 0.27% mutation freq.


######## STOCK #########

######   #Plot the stock virus ######
# Stock is Run0_17
stock<-read.csv("Output/OverviewF_PID/Run0_17_filtered.overview.con.csv",stringsAsFactors=F, row.names = 1)

Ps<-ggplot(data=stock, aes(x=pos, y=freq.mutations*100))+
    ylab("% Divesity")+xlab("")+ylim(0,100)+
    geom_point(size=0.7, color=cols2[4])+theme_bw()+
    ggtitle("Stock")+theme(plot.title = element_text(size=12))
#ggsave(plot=Ps, "Output/Figures/Stock_mf_PID.pdf", width = 7, height = 2)

#How many sites are over 1, 5, 10%?
s2<-stock[stock$freq.mutations.ref>=0.01,]
s2=s2[!is.na(s2$pos),]  #68 positions over 1%
s2<-stock[stock$freq.mutations.ref>=0.005,]
s2=s2[!is.na(s2$pos),]  #95 positions over 0.5%
s2<-stock[stock$freq.mutations.ref>=0.05,]
s2=s2[!is.na(s2$pos),]  #22 positions over 5%
s2<-stock[stock$freq.mutations.ref>=0.1,]
s2=s2[!is.na(s2$pos),]  #16 positions over 10%

df<-stock
#calculate total ns freq at each position
df$ns1<-as.numeric(apply(df[,c("Type.r","freq.Ts.ref")],1, function(x) if (x["Type.r"]=="nonsyn") x["freq.Ts.ref"] else 0))
df$ns2<-as.numeric(apply(df[,c("Type.tv1.r","freq.transv1.ref")],1, function(x) if (x["Type.tv1.r"]=="nonsyn") x=x["freq.transv1.ref"] else 0))
df$ns3<-as.numeric(apply(df[,c("Type.tv2.r","freq.transv2.ref")],1, function(x) if (x["Type.tv2.r"]=="nonsyn") x=x["freq.transv2.ref"] else 0))
df$ns<-df$ns1+df$ns2+df$ns3
#calculate total syn freq at each position
df$syn1<-as.numeric(apply(df[,c("Type.r","freq.Ts.ref")],1, function(x) if (x["Type.r"]=="syn") x["freq.Ts.ref"] else 0))
df$syn2<-as.numeric(apply(df[,c("Type.tv1.r","freq.transv1.ref")],1, function(x) if (x["Type.tv1.r"]=="syn") x=x["freq.transv1.ref"] else 0))
df$syn3<-as.numeric(apply(df[,c("Type.tv2.r","freq.transv2.ref")],1, function(x) if (x["Type.tv2.r"]=="syn") x=x["freq.transv2.ref"] else 0))
df$syn<-df$syn1+df$syn2+df$syn3
#calculate total stop freq at each position
df$stop1<-as.numeric(apply(df[,c("Type.r","freq.Ts.ref")],1, function(x) if (x["Type.r"]=="stop") x["freq.Ts.ref"] else 0))
df$stop2<-as.numeric(apply(df[,c("Type.tv1.r","freq.transv1.ref")],1, function(x) if (x["Type.tv1.r"]=="stop") x=x["freq.transv1.ref"] else 0))
df$stop3<-as.numeric(apply(df[,c("Type.tv2.r","freq.transv2.ref")],1, function(x) if (x["Type.tv2.r"]=="stop") x=x["freq.transv2.ref"] else 0))
df$stop<-df$stop1+df$stop2+df$stop3

df2<-df[,c("pos","ns")]
colnames(df2)[2]<-"freq"
df2$Type<-"nonsyn"
df2<-df2[!is.na(df2$freq),]

df3<-df[,c("pos","syn")]
colnames(df3)[2]<-"freq"
df3$Type<-"syn"
df3<-df3[!is.na(df3$freq),]

#combine syn and nonsyn dataframes
mut<-rbind(df2,df3)

cutoff=0.005
#Plot only the sites with freq > cutoff to avoid too many points
mut$freq[mut$freq<cutoff]<-NA

#select colors
MFcolors<-c("#CA0020","#0571B0")


#nonsyn mean
nonsyn.ave<-format(round(mean(df2$freq,na.rm = T)*100,2),nsmall = 2) 
# 0.9980099
#Syn mean
syn.ave<-format(round(mean(df3$freq,na.rm = T)*100,2),nsmall = 2)
# 0.3761383


ave<-format(round(mean(stock$freq.mutations.ref, na.rm=T)*100,2),nsmall = 2)


ggplot(data=mut, aes(x=pos, y=freq*100, color=Type, fill=Type))+
    geom_bar(stat = "identity", width=0.1)+
    scale_color_manual(values=MFcolors,guide = 'none')+
    scale_fill_manual(values=MFcolors, guide = 'none')+
    scale_y_continuous(limits=c(0,50))+
    theme_bw()+
    theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank())+
    ylab("% Diversity")+xlab('Genome position')+
    theme(legend.title = element_blank())+
    geom_text(data=mut[mut$freq>0.1,], aes(label=pos), hjust=0.3, vjust=-.8, size=2,show.legend = FALSE)+
    annotate(geom="text", x=610, y=49, hjust=0,label=paste0("Nonsyn: ave ",nonsyn.ave,"%"),color =MFcolors[1], size=2.5)+
    annotate(geom="text", x=610,  y=46, hjust=0, label=paste0("Syn: ave ",syn.ave,"%"),color =MFcolors[2], size=2.5)+
    annotate("segment", x = 580, xend = 600, y = 49, yend = 49, colour = MFcolors[1]) +
    annotate("segment", x = 580, xend = 600, y = 46, yend = 46, colour = MFcolors[2])+
    annotate(geom='text', x=610, y=43, label=paste0("Total diversity: ", ave,"%"),color ='gray20', size=2.5,hjust=0)
ggsave("Output/Figures/Stock_Diversity.pdf", width =6 ,height = 3.5)
ggsave("Output/Figures/Stock_Diversity.png", width =6 ,height = 3.5, unit="in", dpi = 300)

ggplot(data=stock, aes(x=pos, y=freq.mutations*100))+
    ylab("% Divesity")+xlab("")+ylim(0,100)+
    geom_point(size=0.7, color=cols2[4])+theme_bw()+
    ggtitle("Stock")+theme(plot.title = element_text(size=12))+
    annotate(geom='text', x=600, y=96, label=paste0("Diversity: ", ave,"%"),color ='gray30', size=3,hjust=0)
ggsave("Output/Figures/Stock_mf_PID.pdf", width = 7, height = 2)


