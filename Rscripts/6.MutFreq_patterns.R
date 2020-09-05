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
cols2<-qualitative_hcl(6, palette="Dark3")

# read the files saved in Overview_output:
SIVFiles_overview<-list.files("Output/OverviewF/",pattern="filtered.overview.csv")

Overview<-list()
for (i in 1:length(SIVFiles_overview)){ 
        overviews<-read.csv(paste0("Output/OverviewF/",SIVFiles_overview[i]),stringsAsFactors=F, row.names = 1)
        Overview[[i]]<-overviews
        names(Overview)[i]<-substr(paste(SIVFiles_overview[i]),start=1,stop=7)
}


# Plot mutation freq across the genome for each week by monkeys

samples<-read.csv("Data/Samples.csv",stringsAsFactors = F)

list.animal<-split(samples, samples$Monkey)

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

for (i in 1:length(monkeyList)){
    sample<-monkeyList[[i]]
    sample<-sample[order(sample$Week),]
    ovDF<-Overview[as.vector(sample$File.name)]
    monkey<-names(monkeyList)[i]

    Plot<-list()
    for (j in 1:length(ovDF)){
        df<-ovDF[[j]]
        week<-sample$Week[j]
        if (week<=17) col=cols[2]
        else {col=cols2[1]}
        Plot[[j]]<-ggplot(data=df, aes(x=pos, y=freq.Ts.ref))+
            ylab("Mutation frequency")+xlab("")+
            geom_point(size=0.7, color=col)+theme_bw()+
            ggtitle(paste0(monkey," Week ",week))+theme(plot.title = element_text(size=12))
    }
    pdf(paste0("Output/MF/",monkey,".pdf"), width = 7, height = length(ovDF)*2)
    do.call(grid.arrange, c(Plot, ncol=1))
    dev.off()
}




for (i in 1:length(monkeyList)){
    sample<-monkeyList[[i]]
    sample<-sample[order(sample$Week),]
    ovDF<-Overview[as.vector(sample$File.name)]
    monkey<-names(monkeyList)[i]
    
    Plot<-list()
    for (j in 1:length(ovDF)){
        df<-ovDF[[j]]
        week<-sample$Week[j]
        if (week<=17) col=cols[2]
        else {col=cols2[1]}
        Plot[[j]]<-ggplot(data=df, aes(x=pos, y=freq.Ts))+
            ylab("Mutation frequency")+xlab("")+ylim(0,0.5)+
            geom_point(size=0.7, color=col)+theme_bw()+
            ggtitle(paste0(monkey," Week ",week))+theme(plot.title = element_text(size=12))
    }
    pdf(paste0("Output/MF/",monkey,".maj.pdf"), width = 7, height = length(ovDF)*2)
    do.call(grid.arrange, c(Plot, ncol=1))
    dev.off()
}


for (i in 1:length(monkeyList)){
    sample<-monkeyList[[i]]
    sample<-sample[order(sample$Week),]
    ovDF<-Overview[as.vector(sample$File.name)]
    monkey<-names(monkeyList)[i]
    
    Plot<-list()
    for (j in 1:length(ovDF)){
        df<-ovDF[[j]]
        week<-sample$Week[j]
        if (week<=17) col=cols[2]
        else {col=cols2[1]}
        Plot[[j]]<-ggplot(data=df, aes(x=pos, y=freq.mutations))+
            ylab("Mutation frequency")+xlab("")+ylim(0,0.5)+
            geom_point(size=0.7, color=col)+theme_bw()+
            ggtitle(paste0(monkey," Week ",week))+theme(plot.title = element_text(size=12))
    }
    pdf(paste0("Output/MF/",monkey,".all.pdf"), width = 7, height = length(ovDF)*2)
    do.call(grid.arrange, c(Plot, ncol=1))
    dev.off()
}

