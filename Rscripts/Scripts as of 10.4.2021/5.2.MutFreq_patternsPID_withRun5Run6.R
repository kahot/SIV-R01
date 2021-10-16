library(ggplot2)
library(gridExtra)
library(colorspace)
source("Rscripts/baseRscript2.R")
cols<-qualitative_hcl(6, palette="Dark3")

# read the files saved in Overview_output:
SIVFiles_overview<-list.files("Output/OverviewF_PIDcon/",pattern=".csv")
SIVfiles<-list.files("Output/Overview_PIDcon/", pattern="Run5|Run6|Run7")
SIVfiles<-SIVfiles[SIVfiles!="Run5_12_overview.csv"]
SIVfiles<-SIVfiles[SIVfiles!="Run5_13_overview.csv"]
SIVfiles<-SIVfiles[SIVfiles!="Run5_15_overview.csv"]
SIVfiles<-c(SIVfiles,"Run0_14_overview.csv")
Overview<-list()
for (i in 1:length(SIVFiles_overview)){ 
        overviews<-read.csv(paste0("Output/OverviewF_PIDcon/",SIVFiles_overview[i]),stringsAsFactors=F, row.names = 1)
        Overview[[i]]<-overviews
        names(Overview)[i]<-substr(paste(SIVFiles_overview[i]),start=1,stop=7)
}
n<-length(SIVFiles_overview)

for (i in (n+1):(length(SIVfiles)+n)){
    df<-read.csv(paste0("Output/Overview_PIDcon/", SIVfiles[i-n]),stringsAsFactors=F, row.names = 1)
    Overview[[i]]<-df
    names(Overview)[i]<-substr(paste(SIVfiles[i-n]),start=1,stop=7)
}


# Plot mutation freq across the genome for each week by monkeys

samples<-read.csv("Data/SamplesNoduplicates.csv",stringsAsFactors = F)
samples$Week<-as.integer(samples$Week)
samples$Tissue<-factor(samples$Tissue, levels=c("Stock","Plasma","LN","HLN","Lung","Unknown"))

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

#TB infection week
tbs<-read.csv("Data/TBinfection.week.csv",stringsAsFactors = F)

monkeys2<-monkeyList[tbs$ids]

#####   #Plot the stock virus ######
stock<-Overview[[which(names(Overview)=="Run0_17")]]
ave<-format(round(mean(stock$freq.mutations.ref, na.rm=T)*100,2),nsmall = 2)

Ps2<-ggplot(data=stock, aes(x=pos, y=freq.Ts.ref))+
    ylab("Mutation frequency")+xlab("")+ylim(0,1)+
    geom_point(size=0.7, color=cols[4])+theme_bw()+
    ggtitle("Stock")+theme(plot.title = element_text(size=12))+
    annotate(geom='text', x=600, y=0.96, label=paste0("MF: ", ave,"%"),color ='gray30', size=3,hjust=0)


#Plot mutation freq (Divergence)
for (i in 10:length(monkeys2)){
    sample<-monkeys2[[i]]
    sample<-sample[order(sample$Week),]
    ovDF<-Overview[as.vector(sample$File.name)]
    monkey<-names(monkeys2)[i]
    tbweek<-tbs$tb[tbs$ids==monkey]
    Plot<-list()
    Plot[[1]]<-Ps2
    for (j in 1:length(ovDF)){
        df<-ovDF[[j]]
        week<-sample$Week[j]
        
        if (week<=tbweek|is.na(tbweek)) col=cols[3]
        else {col=cols[1]}
        ts<-ifelse(sample$Tissue2[j]=="Plasma", "",sample$Tissue2[j])
        ave<-round(mean(df$freq.mutations.ref, na.rm=T)*100,digit= 2)
        div<-round(mean(df$freq.mutations, na.rm=T)*100,2)
        
        Plot[[j+1]]<-ggplot(data=df, aes(x=pos, y=freq.mutations.ref))+
            ylab("Mutation frequency")+xlab("")+ylim(0,1)+
            geom_point(size=0.7, color=col)+theme_bw()+
            ggtitle(paste0(monkey," Week ",week," ",ts))+
            theme(plot.title = element_text(size=11))+
            annotate(geom='text', x=600, y=0.96, label=paste0("MF: ", ave,"%"),color ='gray30', size=3,hjust=0)+
            annotate(geom='text', x=600, y=0.82, label=paste0("Diversity: ", div,"%"),color ='gray30', size=3,hjust=0)
        
    }
    pdf(paste0("Output/MF_PID/withRun7/MF.",monkey,".pdf"), width = 7, height = length(ovDF)*2)
    do.call(grid.arrange, c(Plot, ncol=1))
    dev.off()
}


#Plot minor variant freq (diversity)
for (i in 1:length(monkeys2)){
    sample<-monkeys2[[i]]
    sample<-sample[order(sample$Week),]
    ovDF<-Overview[as.vector(sample$File.name)]
    monkey<-names(monkeys2)[i]
    tbweek<-tbs$tb[tbs$ids==monkey]
    Plot<-list()
    Plot[[1]]<-Ps2
    for (j in 1:length(ovDF)){
        df<-ovDF[[j]]
        week<-sample$Week[j]
        if (week<=tbweek) col=cols[3]
        else {col=cols[1]}
        ts<-ifelse(sample$Tissue2[j]=="Plasma", "",sample$Tissue2[j])
        ave<-round(mean(df$freq.mutations.ref, na.rm=T)*100,digit= 2)
        div<-round(mean(df$freq.mutations, na.rm=T)*100,2)
        
        Plot[[j+1]]<-ggplot(data=df, aes(x=pos, y=freq.mutations*100))+
            ylab("Diversity")+xlab("")+ylim(0,50)+
            geom_point(size=0.7, color=col)+theme_bw()+
            ggtitle(paste0(monkey," Week ",week," ",ts))+
            theme(plot.title = element_text(size=11))+
            annotate(geom='text', x=600, y=47, label=paste0("MF: ", ave,"%"),color ='gray30', size=3,hjust=0)+
            annotate(geom='text', x=600, y=40, label=paste0("Diversity: ", div,"%"),color ='gray30', size=3,hjust=0)
        
    }
    pdf(paste0("Output/MF_PID/withRun7/MVF.",monkey,".pdf"), width = 7, height = length(ovDF)*2)
    do.call(grid.arrange, c(Plot, ncol=1))
    dev.off()
}



### LOOK at the read depth vs. trasnv freq (for low quality samples (Run2))#### 
## Run2 samples are messy in the middle (pos 485 - 550ish) ##
df<-Overview[["Run2_10"]]
mean(df$freq.Ts.ref, na.rm=T)
#0.01476899
mean(df$freq.transv.ref, na.rm=T)
#0.02486384

plot(df$freq.transv.ref)
plot(df$TotalReads)

pt<-list()
pt[[1]]<-ggplot(data=df, aes(x=pos, y=freq.mutations*100))+
    ylab("Diversity")+xlab("")+ylim(0,50)+
    geom_point(size=0.7, color=col)+theme_bw()
pt[[2]]<-ggplot(data=df, aes(x=pos, y=freq.Ts*100))+
    ylab("Diversity Ts")+xlab("")+ylim(0,50)+
    geom_point(size=0.7, color=col)+theme_bw()
pt[[3]]<-ggplot(data=df, aes(x=pos, y=freq.transv*100))+
    ylab("Diversity Tvs")+xlab("")+ylim(0,50)+
    geom_point(size=0.7, color=col)+theme_bw()

pt[[4]]<-ggplot(data=df, aes(x=pos, y=TotalReads))+
    ylab("Total Reads")+xlab("")+
    geom_point(size=0.7, color="blue")+theme_bw()

pdf(paste0("Output/MF_PID/withRun5.6/MVF.vs.Reads.A22117.week16.pdf"), width = 7, height = 8)
do.call(grid.arrange, c(pt, ncol=1))
dev.off()

pt<-list()
for (i in 1: length(Overview)){
    df<-Overview[[i]]
    filename<-names(Overview)[i]
    mk<-samples$Monkey[samples$File.name==filename]
    wk<-samples$Week[samples$File.name==filename]
    pt[[i]]<-ggplot(data=df, aes(x=pos, y=TotalReads))+
        ylab("Total Reads")+xlab("")+
        geom_point(size=0.7, color="blue")+theme_bw()+
        ggtitle(paste(filename, mk, "Week", wk))+
        theme(plot.title = element_text(size=11))
}
pdf(paste0("Output/MF_PID/withRun5.6/MVF.vs.Reads.allfiles.pdf"), width = 15, height = 97)
do.call(grid.arrange, c(pt, ncol=3))
dev.off()



# Plot Ts mutation freq 
for (i in 1:length(monkeys2)){
    sample<-monkeys2[[i]]
    sample<-sample[order(sample$Week),]
    ovDF<-Overview[as.vector(sample$File.name)]
    monkey<-names(monkeys2)[i]
    tbweek<-tbs$tb[tbs$ids==monkey]
    Plot<-list()
    Plot[[1]]<-Ps2
    for (j in 1:length(ovDF)){
        df<-ovDF[[j]]
        week<-sample$Week[j]
        if (week<=tbweek) col=cols[3]
        else {col=cols[1]}
        ts<-ifelse(sample$Tissue2[j]=="Plasma", "",sample$Tissue2[j])
        ave<-round(mean(df$freq.mutations.ref, na.rm=T)*100,digit= 2)
        div<-round(mean(df$freq.mutations, na.rm=T)*100,2)
        
        Plot[[j+1]]<-ggplot(data=df, aes(x=pos, y=freq.Ts.ref))+
            ylab("Transition MF")+xlab("")+ylim(0,1)+
            geom_point(size=0.7, color=col)+theme_bw()+
            ggtitle(paste0(monkey," Week ",week," ",ts))+
            theme(plot.title = element_text(size=11))+
            annotate(geom='text', x=600, y=0.96, label=paste0("MF: ", ave,"%"),color ='gray30', size=3,hjust=0)+
            annotate(geom='text', x=600, y=0.82, label=paste0("Diversity: ", div,"%"),color ='gray30', size=3,hjust=0)
        
    }
    pdf(paste0("Output/MF_PID/withRun5.6/TsMF.",monkey,".pdf"), width = 7, height = length(ovDF)*2)
    do.call(grid.arrange, c(Plot, ncol=1))
    dev.off()
}

############################
## Diversity Shift over time ##

monkeys2<-monkeyList[tbs$ids]
morder<-c("A22517","A22617","A22117","A22317","A22217","A21918","A23918","A34019","A34119","A34219")
plots<-list()
for (i in 1:length(morder)){
    sample<-monkeys2[[morder[i]]]
    sample<-sample[order(sample$Week),]
    ovDF<-Overview[as.vector(sample$File.name)]
    monkey<-morder[i]
    tbweek<-tbs$tb[tbs$ids==monkey]
    diversity<-sample[,c(2,3,5,7)]
    monkey<-gsub("A",'',monkey)
    for (j in 1:length(ovDF)){
        df<-ovDF[[j]]
        diversity$Diversity[j]<-mean(df$freq.mutations, na.rm=T)*100
        diversity$MF[j]<-mean(df$freq.mutations.ref, na.rm=T)*100
    }
    diversity<-InsertRow(diversity,c("Stock","Stock",0,"Stock",mean(stock$freq.mutations.ref, na.rm=T)*100,mean(stock$freq.mutations.ref, na.rm=T)*100 ),1)
    diversity$Week<-as.integer(diversity$Week)
    diversity$Diversity<-as.numeric(diversity$Diversity)
    diversity$MF<-as.numeric(diversity$MF)
    #diversity$Tissue<-factor(diversity$Tissue, levels=c("Stock","Plasma","LN","HLN","Unknown"))
    plots[[i]]<-ggplot()+
            ylab("Diversity")+xlab("")+
            geom_point(data=diversity, aes(x=Week, y=Diversity, color=Tissue), size=2)+theme_bw()+
            ggtitle(paste0(monkey))+ylim(0.4,2.2)+
            theme(plot.title = element_text(size=11))+
            scale_color_manual(values=cols[c(1,6,4,2,3)])+
            theme(panel.grid.major.x  = element_blank(),panel.grid.minor.x = element_blank())+
            geom_vline(xintercept=tbweek, col="deeppink")+
            scale_x_continuous(breaks=seq(0,32,2), limits = c(0,32))+
            geom_line(data=diversity[diversity$Tissue=="Plasma"|diversity$Tissue=="Stock",],aes(x=Week, y=Diversity), color=cols[6])
    
    
    
    
}
pdf("Output/MF_PID/withRun7/Diversity_overtime_all.pdf", width = 8.5, height = 11)
do.call(grid.arrange, c(plots, ncol=2))
dev.off()






