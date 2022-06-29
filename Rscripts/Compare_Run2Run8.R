library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)
library(colorspace)
library(DataCombine)
library(reshape2)
#source("Rscripts/baseRscript2.R")
cols8<-qualitative_hcl(8, palette="Dark3")

# read all Overview files (procssed without filtering:
#Use 'Filtered' overviews for Run1-Run4 (except for Run0_14 (500 reads) and some Run5 files with high read depths) 
files2<-list.files("Output/OverviewF_PID/",pattern="Run2")
files8<-list.files("Output/OverviewF_PID/",pattern="Run8")


ov2<-list()
for (i in 1:length(files2)){ 
        overviews<-read.csv(paste0("Output/OverviewF_PID/",files2[i]),stringsAsFactors=F, row.names = 1)
        ov2[[i]]<-overviews
        names(ov2)[i]<-substr(paste(files2[i]),start=1,stop=7)
}

ov8<-list()
for (i in 1:length(files8)){ 
    overviews<-read.csv(paste0("Output/OverviewF_PID/",files8[i]),stringsAsFactors=F, row.names = 1)
    ov8[[i]]<-overviews
    names(ov8)[i]<-substr(paste(files8[i]),start=1,stop=7)
}

files2
files8

samples<-read.csv("Data/SamplesplusRun8.csv")
samples<-samples[samples$Run==2|samples$Run==8,]

#Plot the old version and new version (trimmed at Q15 and filtered >100)
for (i in 1:10){
    sample<-samples$File.name[i]
    sample2<-samples$File.name[samples$Monkey==samples$Monkey[i]&samples$Week==samples$Week[i]&samples$Run==8]
    
    df1<-ov2[[sample]]
    df2<-ov8[[sample2]]
    
    
    ave1<-format(round(mean(df1$freq.mutations.ref, na.rm=T)*100,2),nsmall = 2)
    ave2<-format(round(mean(df2$freq.mutations.ref, na.rm=T)*100,2),nsmall = 2)
    
    p1<-ggplot(data=df1, aes(x=pos, y=freq.mutations.ref*100))+
        ylab("% Mutation")+xlab("")+ylim(0,100)+
        geom_point(size=0.8, color=cols8[5])+theme_bw()+
        ggtitle(paste(sample,samples$Monkey[i], samples$Week[i]))+theme(plot.title = element_text(size=12))+
        annotate(geom='text', x=600, y=96, label=paste0("Mean: ", ave1,"%"),color ='gray20', size=3.3,hjust=0)
    p2<-ggplot(data=df2, aes(x=pos, y=freq.mutations.ref*100))+
        ylab("% Mutation")+xlab("")+ylim(0,100)+
        geom_point(size=0.8, color=cols8[8])+theme_bw()+
        ggtitle(paste(sample2,samples$Monkey[i], samples$Week[i]))+theme(plot.title = element_text(size=12))+
        annotate(geom='text', x=600, y=96, label=paste0("Mean: ", ave2,"%"),color ='gray20', size=3.3,hjust=0)
    
    df1<-df1[,c("pos","freq.mutations.ref")]
    colnames(df1)[2]<-"Run2"
    df2<-df2[,c("pos","freq.mutations.ref")]
    colnames(df2)[2]<-"Run8"
    
    df<-merge(df1,df2, by="pos")
    dfm<-melt(df, id.vars="pos")
    
    p3<-ggplot(data=dfm, aes(x=pos, y=value*100, color=variable))+
        ylab("% Mutation")+xlab("")+ylim(0,100)+
        geom_point(size=0.9)+theme_bw()+
        scale_color_manual(values=c(paste0(cols8[5],"B3"),paste0(cols8[8],"66")))+
        theme(legend.position = "bottom", legend.title=element_blank())
    
    pdf(paste0("Output/QC3/Run2.vs.Run8_", sample2, ".pdf"), height = 7, width = 6)
    grid.arrange(p1,p2,p3, ncol=1)
    dev.off()
}


