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

## second run (Ambrose1 5-18-18 has much better quality) A22617 week8 (Run1_12, rerun at Q30, cutoff=1000) i=6
#Plot mut freq relative to ref seq
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

######   #Plot the stock virus ######
df<-Overview[[9]]
ggplot(data=df, aes(x=pos, y=freq.Ts.ref))+
    ylab("Mutation frequency")+xlab("")+
    geom_point(size=0.7, color=col)+theme_bw()+
    ggtitle("Stock")+theme(plot.title = element_text(size=12))
ggsave("Output/Stock_mf.pdf", width = 7, height = 2)




#Plot minor variant freq
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

# Plot all mutation freq
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
        Plot[[j]]<-ggplot(data=df, aes(x=pos, y=freq.mutations.ref))+
            ylab("Mutation frequency")+xlab("")+ylim(0,1)+
            geom_point(size=0.7, color=col)+theme_bw()+
            ggtitle(paste0(monkey," Week ",week))+theme(plot.title = element_text(size=12))
    }
    pdf(paste0("Output/MF/",monkey,".all.pdf"), width = 7, height = length(ovDF)*2)
    do.call(grid.arrange, c(Plot, ncol=1))
    dev.off()
}



###### 
## Add nucleotide to the position info
i==3,5,6
for (i in 1:length(monkeyList)){
    sample<-monkeyList[[i]]
    sample = sample[order(sample[,'Week'],-sample[,'Run']),]
    sample = sample[!duplicated(sample$Week),]
    sample<-sample[order(sample$Week),]
    ovDF<-Overview[as.vector(sample$File.name)]
    monkey<-names(monkeyList)[i]
   
    sites<-lapply(ovDF, function(x) x$pos[x$freq.Ts.ref>0.1])
    sites<-unique(unlist(sites))
    sites<-sites[!is.na(sites)]
    #select only the positions in sites
    mfDF<-lapply(ovDF, function(x) x<-x[x$pos %in% sites,])
    #select the transition mutation freqeuncies of all files
    Ts<-unname(sapply(mfDF, `[[`, "freq.Ts.ref"))
    Ts<-data.frame(t(Ts))
    colnames(Ts)<-paste0("Pos_",sites)
    Ts$Week<-sample$Week
    Tsm<-melt(Ts, id.vars = "Week")
    colnames(Tsm)[2:3]<-c("Position", "MF")
    rown<-ceiling(length(sites)/5)
    ggplot(data=Tsm, aes(x=Week, y=MF))+
        ylab("Mutation frequency")+xlab("Weeks")+
        facet_wrap(~ Position, nrow=rown, ncol=5)+
        geom_point(size=2, color=col)+theme_bw()+ggtitle(paste(monkey))+
        geom_vline(xintercept=17, col="blue")+
        ggsave(paste0("Output/Timeseries/", monkey, ".pdf"), heigh=rown*2, width =10)
}




## Select th ones that are interesting

#Don't have post infecion files for A22117

ids<-c("A21918","A22317","A22517", "A22617","A22217","A23918","A22117")

monkeys2<-monkeyList[ids]
Pos<-list()
for (i in 1:length(monkeys2)){
    sample<-monkeys2[[i]]
    sample = sample[order(sample[,'Week'],-sample[,'Run']),]
    sample = sample[!duplicated(sample$Week),]
    #sample<-sample[order(sample$Week),]
    ovDF<-Overview[as.vector(sample$File.name)]
    monkey<-names(monkeys2)[i]
    
    sites<-lapply(ovDF, function(x) x$pos[x$freq.Ts.ref>0.1])
    sites<-unique(unlist(sites))
    sites<-sites[!is.na(sites)]
    sites<-sites[order(sites)]
    #select only the positions in sites
    mfDF<-lapply(ovDF, function(x) x<-x[x$pos %in% sites,])
    nt<-mfDF[[1]][,"ref"]
    Pos[[i]]<-sites
    names(Pos)[[i]]<-monkey
    Ts<-unname(sapply(mfDF, `[[`, "freq.Ts.ref"))
    Ts<-data.frame(t(Ts))
    colnames(Ts)<-paste0("Pos_",sites, ' ',nt)
    Ts$Week<-sample$Week
    
    #nt<-unname(sapply(mfDF, `[[`, "MajNt"))
    #nt<-    unname(sapply(mfDF, `[[`, "ref"))
    #nt<-data.frame(t(nt))
    #colnames(nt)<-paste0("Pos_",sites)
    #nt$Week<-sample$Week
    
    Tsm<-melt(Ts, id.vars = "Week")
    colnames(Tsm)[2:3]<-c("Position", "MF")
    rown<-ceiling(length(sites)/5)
    
    #Plot<-list()
    #for (j in 1:length(ovDF)){
    #    df<-ovDF[[j]]
    #    week<-sample$Week[j]
    #    if (week<=17) col=cols[2]
    #    else {col=cols2[1]}
    #    Plot[[j]]<-ggplot(data=df, aes(x=pos, y=freq.Ts.ref))+
    #        ylab("Mutation frequency")+xlab("")+ylim(0,1)+
    #        geom_point(size=0.7, color=col)+theme_bw()+
    #        ggtitle(paste0(monkey," Week ",week))+theme(plot.title = element_text(size=12))
    #}
    #pdf(paste0("Output/MF/",monkey,".nodup.pdf"), width = 7, height = length(ovDF)*2)
    #do.call(grid.arrange, c(Plot, ncol=1))
    #dev.off()
    #
    #ggplot(data=Tsm, aes(x=Week, y=MF))+
    #    ylab("Mutation frequency")+xlab("Weeks")+
    #    facet_wrap(~ Position, nrow=rown, ncol=5)+
    #    geom_point(size=2, color=col)+theme_bw()+ggtitle(paste(monkey))+
    #    geom_vline(xintercept=17, col="blue")+
    #    ggsave(paste0("Output/Timeseries2/", monkey, ".pdf"), heigh=rown*2, width =10)
}

#Find the over lapping positions
positions<-unique(unlist(unname(Pos)))
positions<-positions[order(positions)]
cpos<-data.frame(Pos=positions)
for (i in 1:length(monkeys2)){
    vec<-Pos[[i]]

    cpos[,names(monkeys2)[i]]<-apply(cpos['Pos'],2, function(x) ifelse(x %in% vec, "Y", "N"))
}

write.csv(cpos, "Output/HighMutfreq_overtime.csv")
cpos$codon<-apply(cpos["Pos"], 1, function(x) if(x%%3==0) x=3 else x=x%%3)

df<-OverviewDF[OverviewDF$pos %in% positions,]
df<-df[,c("pos","ref","Type.r","Type","WTAA","MUTAA","bigAAChange","makesCpG")]
colnames(cpos)[1]<-"pos"
cpos<-merge(cpos, df, by="pos")

write.csv(cpos, "Output/HighMutfreq_sites.csv")

