library(ggplot2)
library(gridExtra)
library(colorspace)
source("Rscripts/baseRscript2.R")
cols<-qualitative_hcl(6, palette="Dark3")

# read the files saved in Overview_output:
SIVFiles_overview<-list.files("Output/OverviewF_PIDcon/",pattern="Run2")
Overview<-list()
for (i in 1:length(SIVFiles_overview)){ 
        overviews<-read.csv(paste0("Output/OverviewF_PIDcon/",SIVFiles_overview[i]),stringsAsFactors=F, row.names = 1)
        Overview[[i]]<-overviews
        names(Overview)[i]<-substr(paste(SIVFiles_overview[i]),start=1,stop=7)
}

# Plot mutation freq across the genome for each week by monkeys

#Plot mutation freq (Divergence)
Plot<-list()
Plot2<-list()
for (i in 1:length(Overview)){
    df<-Overview[[i]]
    #find the total read counts
    m<-df$TotalReads[100]
    #remove the sites that have 1000 or less total reads than 'm'
    lowreads<-which(df$TotalReads<=m-1000)
    df[lowreads,c(7:16)]<-NA
    
    #Plot[[i]]<-ggplot(data=df, aes(x=pos, y=freq.mutations.ref))+
    #    ylab("Mutation frequency")+xlab("")+ylim(0,1)+
    #    geom_point(size=0.7)+theme_bw()
    #Plot2[[i]]<-ggplot(data=df, aes(x=pos, y=TotalReads))+
    #    ylab("Read count")+xlab("")+
    #    geom_point(size=0.7)+theme_bw()

 write.csv(df,paste0("Output/OverviewF_PIDcon/",SIVFiles_overview[i]))
}

pdf(paste0("Output/OverviewF_PID/Center_removed_run2.pdf"), width = 7, height = 15)
do.call(grid.arrange, c(Plot, ncol=2))
dev.off()

pdf(paste0("Output/OverviewF_PID/TotalReads_run2.pdf"), width = 7, height = 15)
do.call(grid.arrange, c(Plot2, ncol=2))
dev.off()

