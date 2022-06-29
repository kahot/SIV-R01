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
df<- read.csv("Output/OverviewF/Control_filtered.overview.csv",stringsAsFactors=F, row.names = 1)

ave<-round(mean(df$freq.mutations.ref, na.rm=T)*100,4)
ggplot(df, aes(x=pos, y=freq.mutations.ref*100))+
    ylab("Mutation frequency %")+xlab("")+
    geom_point(size=0.7, color="blue")+theme_bw()+
    ggtitle("Control (raw data)")+theme(plot.title = element_text(size=12))+
    annotate(geom='text', x=590, y=3.1, label=paste0("Ave. ", ave, "%"),color ='gray20', size=3.5,hjust=0)
ggsave("Output/Control/Control_Q35.pdf", width = 7, height = 2)

#which one has the high freq?
df$pos[df$freq.mutations.ref>=0.02] #position 197

#Removing pos 197?
df2<-df[df$pos!=197,]
mean(df2$freq.mutations.ref, na.rm=T)  #0.00156201
median(df2$freq.mutations.ref, na.rm = T) #0.001477813

