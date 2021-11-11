library(reshape2)
library(gridExtra)
library(colorspace)
library(cowplot)
library(DataCombine)
library(dplyr)
source("Rscripts/baseRscript.R")
cols2<-qualitative_hcl(6, palette="Dark3")


######   #Plot the stock virus ######
# Stock is Run0_17
stock<-read.csv("Output/OverviewF_PIDcon/Run0_17_filtered.overview.con.csv",stringsAsFactors=F, row.names = 1)

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
MFcolors<-c("#fb8072","#FF9300","#9437FF")

ggplot(data=mut, aes(x=pos, y=freq*100, color=Type, fill=Type))+
    geom_bar(stat = "identity", width=0.1)+
    scale_color_manual(values=MFcolors[2:3],guide = 'none')+
    scale_fill_manual(values=MFcolors[2:3], guide = 'none')+
    scale_y_continuous(limits=c(0,50))+
    theme_bw()+
    theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank())+
    ylab("% Diversity")+xlab('Genome position')+
    theme(legend.title = element_blank())+
    geom_text(data=mut[mut$freq>0.1,], aes(label=pos), hjust=0.3, vjust=-.8, size=2,show.legend = FALSE)+
    annotate(geom="text", x=640, y=49, hjust=0,label="nonsyn",color =MFcolors[2], size=2.5)+
    annotate(geom="text", x=640,  y=46, hjust=0, label="syn",color =MFcolors[3], size=2.5)+
    annotate("segment", x = 610, xend = 630, y = 49, yend = 49, colour = MFcolors[2]) +
    annotate("segment", x = 610, xend = 630, y = 46, yend = 46, colour = MFcolors[3]) 
ggsave("Output/Figures/Stock_Diversity.pdf", width =6 ,height = 3.5)


#nonsyn mean
mean(df2$freq,na.rm = T) 
# 0.009980099
#Syn mean
mean(df3$freq,na.rm = T)
# 0.003761383

ave<-format(round(mean(stock$freq.mutations.ref, na.rm=T)*100,2),nsmall = 2)
#1.44%


Ps<-ggplot(data=stock, aes(x=pos, y=freq.mutations*100))+
    ylab("% Divesity")+xlab("")+ylim(0,100)+
    geom_point(size=0.7, color=cols2[4])+theme_bw()+
    ggtitle("Stock")+theme(plot.title = element_text(size=12))+
    annotate(geom='text', x=600, y=96, label=paste0("Diversity: ", ave,"%"),color ='gray30', size=3,hjust=0)
ggsave(plot=Ps, "Output/Figures/Stock_mf_PID.pdf", width = 7, height = 2)


