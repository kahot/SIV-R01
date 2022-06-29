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
filesF<-list.files("Output/OverviewF_PIDold/",pattern=".csv")

files<-list.files("Output/Overview_PIDold/", pattern="Run5|Run6|Run7")
files<-files[files!="Run5_12_overview.csv"]
files<-files[files!="Run5_13_overview.csv"]
files<-files[files!="Run5_15_overview.csv"]
files<-c(files,"Run0_14_overview.csv")

Old<-list()
for (i in 1:length(filesF)){ 
        overviews<-read.csv(paste0("Output/OverviewF_PIDold/",filesF[i]),stringsAsFactors=F, row.names = 1)
        Old[[i]]<-overviews
        names(Old)[i]<-substr(paste(filesF[i]),start=1,stop=7)
}
n<-length(filesF)

for (i in (n+1):(length(files)+n)){
    df<-read.csv(paste0("Output/Overview_PIDold/", files[i-n]),stringsAsFactors=F, row.names = 1)
    Old[[i]]<-df
    names(Old)[i]<-substr(paste(files[i-n]),start=1,stop=7)
}

####
#Read the filtered and PID processed version
nonf<-list.files("Output/Overview_PID/",pattern=".csv")

non<-list()
fil10<-list()
for (i in 1:length(nonf)){ 
    df<-read.csv(paste0("Output/Overview_PID/",nonf[i]),stringsAsFactors=F, row.names = 1)
    non[[i]]<-df
    names(non)[i]<-substr(paste(nonf[i]),start=1,stop=7)
    remove<-which(df$TotalReads<10)
    df[remove, 7:16]<-NA
    fil10[[i]]<-df
    names(fil10)[i]<-substr(paste(nonf[i]),start=1,stop=7)
    
}

filf<-list.files("Output/OverviewF_PID/",pattern=".csv")
fil100<-list()
for (i in 1:length(filf)){ 
    df<-read.csv(paste0("Output/OverviewF_PID/",filf[i]),stringsAsFactors=F, row.names = 1)
    fil100[[i]]<-df
    names(fil100)[i]<-substr(paste(filf[i]),start=1,stop=7)
}



#Plot the old version and new version (trimmed at Q15 and filtered >100)
for (i in 1: length(fil100)){
    id<-names(fil100)[i]
    newDf<-fil100[[id]]
    oldDf<-Old[[id]]
    ave1<-format(round(mean(newDf$freq.mutations.ref, na.rm=T)*100,2),nsmall = 2)
    ave2<-format(round(mean(oldDf$freq.mutations.ref, na.rm=T)*100,2),nsmall = 2)
    
    p1<-ggplot(data=newDf, aes(x=pos, y=freq.mutations.ref*100))+
        ylab("% Mutation")+xlab("")+ylim(0,100)+
        geom_point(size=0.8, color=cols8[5])+theme_bw()+
        ggtitle("Trimmed at Q15")+theme(plot.title = element_text(size=12))+
        annotate(geom='text', x=600, y=96, label=paste0("Mean: ", ave1,"%"),color ='gray20', size=3.3,hjust=0)
    p2<-ggplot(data=oldDf, aes(x=pos, y=freq.mutations.ref*100))+
        ylab("% Mutation")+xlab("")+ylim(0,100)+
        geom_point(size=0.8, color=cols8[8])+theme_bw()+
        ggtitle("Non-trimmed")+theme(plot.title = element_text(size=12))+
        annotate(geom='text', x=600, y=96, label=paste0("Mean: ", ave2,"%"),color ='gray20', size=3.3,hjust=0)
    
    newDf<-newDf[,c("pos","freq.mutations.ref")]
    colnames(newDf)[2]<-"Trimmed"
    oldDf<-oldDf[,c("pos","freq.mutations.ref")]
    colnames(oldDf)[2]<-"Non-trimmed"
    
    df<-merge(newDf,oldDf, by="pos")
    dfm<-melt(df, id.vars="pos")
    
    p3<-ggplot(data=dfm, aes(x=pos, y=value*100, color=variable))+
        ylab("% Mutation")+xlab("")+ylim(0,100)+
        geom_point(size=0.9)+theme_bw()+
        scale_color_manual(values=c(paste0(cols8[5],"B3"),paste0(cols8[8],"66")))+
        theme(legend.position = "bottom", legend.title=element_blank())
    
    pdf(paste0("Output/QC2/Compare_", id, ".pdf"), height = 7, width = 6)
    grid.arrange(p1,p2,p3, ncol=1)
    dev.off()
}


# Filtered at 10 vs. 100

for (i in 1: length(filf)){
    for (i in 11:20){
    id<-names(fil100)[i]
    nonDF<-non[[id]]
    DF10<-fil10[[id]]
    filDF<-fil100[[id]]
    ave1<-format(round(mean(nonDF$freq.mutations.ref, na.rm=T)*100,2),nsmall = 2)
    ave2<-format(round(mean(DF10$freq.mutations.ref, na.rm=T)*100,2),nsmall = 2)
    ave3<-format(round(mean(filDF$freq.mutations.ref, na.rm=T)*100,2),nsmall = 2)
    
    p1<-ggplot(data=nonDF, aes(x=pos, y=freq.mutations.ref*100))+
        ylab("% Mutation")+xlab("")+ylim(0,100)+
        geom_point(size=0.8, color=cols8[7])+theme_bw()+
        ggtitle(paste0(id, " Not filtered"))+theme(plot.title = element_text(size=12))+
        annotate(geom='text', x=600, y=96, label=paste0("Mean: ", ave1,"%"),color ='gray20', size=3.3,hjust=0)
    
    p2<-ggplot(data=DF10, aes(x=pos, y=freq.mutations.ref*100))+
        ylab("% Mutation")+xlab("")+ylim(0,100)+
        geom_point(size=0.8, color=cols8[6])+theme_bw()+
        ggtitle("Filtered for sites with reads >10")+theme(plot.title = element_text(size=12))+
        annotate(geom='text', x=600, y=96, label=paste0("Mean: ", ave2,"%"),color ='gray20', size=3.3,hjust=0)
    p3<-ggplot(data=filDF, aes(x=pos, y=freq.mutations.ref*100))+
        ylab("% Mutation")+xlab("")+ylim(0,100)+
        geom_point(size=0.8, color=cols8[5])+theme_bw()+
        ggtitle("Filtered for sites with reads >100")+theme(plot.title = element_text(size=12))+
        annotate(geom='text', x=600, y=96, label=paste0("Mean: ", ave3,"%"),color ='gray20', size=3.3,hjust=0)
    
    pdf(paste0("Output/QC2/over100/", id, ".pdf"), height = 4.7, width = 6)
    grid.arrange(p1,p2,p3, ncol=1)
    dev.off()
}



#Read depth comparison

ori<-read.csv("Output/ReadDeapth_all_old.csv",row.names = 1)
new<-read.csv("Output/ReadDeapth_all.csv",row.names = 1)

depth<-ori[, c("File.name","Average","Max","Tissue")]
colnames(depth)<-c("File.name","Average.ori","Max.ori","Tissue")
depth<-merge(depth, new[,c("File.name","Average","Max")], by="File.name")

ave<-depth[,c(1,2,5)]
colnames(ave)[2:3]<-c("Non-trimmed","Trimmed")
avem<-melt(ave)
ggplot(avem, aes(x=File.name, y=value, color=variable, fill=variable))+
    geom_bar(position=position_dodge(width = 0.8), stat = "identity", alpha=0.6)+
    theme_classic()+xlab('')+ylab('Read depth')+
    theme(axis.text.x = element_text(angle=90, size=5), legend.title = element_blank())
ggsave("Output/QC2/Read_depth_comparison_old.vs.new.pdf", width = 14, height = 4)




stock<-non[["Run0_17"]]
ggplot(stock, aes(x=pos, y=TotalReads))+
    geom_bar(stat="identity")+ylab("Total reads")+xlab('')
ggsave("Output/QC2/Stock_depth_notFiltered.pdf", width = 7, height = 3.5)


for (i in 1:length(nonf)){
    df<-non[[i]]
    ggplot(df, aes(x=pos, y=TotalReads))+
        geom_bar(stat="identity")+ylab("Total reads")+xlab('')+ggtitle(names(non)[i])
    ggsave(paste0("Output/QC2/Reads/", names(non)[i],"_depth_notFiltered.pdf"), width = 7, height = 3.5)
}



#Cutoff 10 vs. 100 for filtering 

diversity<-data.frame(ID=names(fil100))
for (i in 1: length(filf)){
    id<-names(fil100)[i]
    nonDF<-non[[id]]
    DF10<-fil10[[id]]
    filDF<-fil100[[id]]
    diversity$No.filter[i]       <-mean(nonDF$freq.mutations.ref, na.rm=T)*100
    diversity$Filter.at.10[i] <-mean(DF10$freq.mutations.ref, na.rm=T)*100
    diversity$Filter.at.100[i]<-mean(filDF$freq.mutations.ref, na.rm=T)*100
}
    
div.m<-melt(diversity, id.vars="ID")
div.m$Run<-substr(div.m$ID, 4,4)

ggplot(div.m,aes(x=ID, y=value, color=variable))+
    geom_point()

runs<-unique(div.m$Run)
P<-list()
for (i in 1:7){
    run<-runs[i]
    df<-div.m[div.m$Run==run,]
    ave.diff<-format(round(mean(df$value[df$variable=="Filter.at.10"]-df$value[df$variable=="Filter.at.100"], na.rm = T),3),nsmall = 3)
    ymax<-max(df$value, na.rm=T)
    
    
    P[[i]]<-ggplot(df,aes(x=ID, y=value, color=variable))+
        geom_point(position=position_dodge(width=0.7))+
        ylab('Diversity')+xlab('')+ggtitle(paste0("Run ",run))+
        annotate(geom='text', x=0.5, y=ymax, label=paste0("Ave diff: ", ave.diff,"%"),
                 color ='gray20',size=3,hjust=0)+
        theme(axis.text.x = element_text(angle=90), legend.position = "none")
    
}    

p<-ggplot(div.m[div.m$Run==0,],aes(x=ID, y=value, color=variable))+
    geom_point(position=position_dodge(width=0.7))+
    ylab('Diversity')+xlab('')+ggtitle(paste0("Run ",run))+
    theme(axis.text.x = element_text(angle=90), legend.title = element_blank())
P[[8]]<- cowplot::get_legend(p)

pdf(paste0("Output/QC2/over100/FilteringResults.compareByRun.pdf"), height =10, width =8)
do.call(grid.arrange, c(P, ncol=2))
dev.off()


##  Conclusion: Filter at 10 is fine except for Run2




SIVFiles<-list.files("Output/Overview_PID/")
Overview<-list()
for (i in 1:length(SIVFiles)){ 
    overviews<-read.csv(paste0("Output/Overview_PID/",SIVFiles[i]),stringsAsFactors=F, row.names = 1)
    Overview[[i]]<-overviews
    names(Overview)[i]<-substr(paste(SIVFiles[i]),start=1,stop=7)
}


Plot<-list()
Plot2<-list()

#Remove the low read sites from Run2 samples and re-save.

gaps<-data.frame(ID=names(Overview))
for (i in 1:length(Overview)){
    df<-Overview[[i]]
    id<-names(Overview)[i]
    remove<-which(df$TotalReads<10)
    df[remove, 7:16]<-NA
    #write.csv(df,paste0("Output/OverviewF_PID/",id,"_filtered.overview.csv"))
    
    df<-df[df$pos>=215,]
    gaps$gap.size[i]<-nrow(df[is.na(df$freq.mutations.ref),])
    
}
    
gaps$Run<-substr(gaps$ID, 4,4)

runcolors<-sapply(gaps$Run, function(x) {
                                    if (x=="2"|x=="4"|x=="6") "black"
                                    else "blue"
})

ggplot(data=gaps, aes(x=ID, y=gap.size))+
    ylab("Gap size")+xlab("")+
    geom_point(size=0.7)+theme_classic()+
    theme(axis.text.x = element_text(angle=90, color=runcolors), legend.title = element_blank())

ggplot(data=gaps, aes(x=ID, y=gap.size))+
    ylab("Gap size")+xlab("")+
    geom_point(size=0.7)+theme_classic()+
    theme(axis.text.x = element_text(angle=90, color=runcolors, size=5), legend.title = element_blank())
ggsave("Output/QC2/Gapsize_filterat10.pdf",width = 8, height = 4)
  