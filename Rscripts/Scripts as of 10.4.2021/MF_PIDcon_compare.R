library(ggplot2)
library(gridExtra)
library(colorspace)
source("Rscripts/baseRscript2.R")
cols<-qualitative_hcl(6, palette="Dark3")

# read the files saved in Overview_output:
files<-list.files("Output/OverviewF_PIDcon/",pattern=".csv")
#files2<-list.files("Output/OverviewF_PIDcon/",pattern=".csv")

PID<-list()
PIDcon<-list()
for (i in 1:length(files)){ 
    fname<-substr(paste(files[i]),start=1,stop=7)
    PIDcon[[i]]<-read.csv(paste0("Output/OverviewF_PIDcon/",files[i]),stringsAsFactors=F, row.names = 1)
    names(PIDcon)[i]<-fname
    PID[[i]]<-read.csv(paste0("Output/OverviewF_PID/",fname,"_filtered.overview.csv"),stringsAsFactors=F, row.names = 1)
    names(PID)[i]<-fname
}


for (i in 1:length(PID)){
    fname<-names(PID)[i]
    df1<-PID[[i]]
    df1$Process<-"PID"
    df2<-PIDcon[[i]]
    df2$Process<-"PID consensus"
    DF<-rbind(df1,df2)
    
    ggplot(data=DF, aes(x=pos, y=freq.mutations))+
        ylab("MVF")+xlab("")+ylim(0,0.5)+
        geom_point(size=1, color=paste0(cols[4],"66"))+theme_bw()+
        facet_wrap(~ Process, nrow=2, ncol=1)+
        ggsave(paste0("Output/MF_PID/PIDvsPIDcon/",fname,".pdf"), width = 7,height = 4)
        
    ggplot(data=DF, aes(x=pos, y=freq.mutations, color=Process))+
            ylab("MVF")+xlab("")+ylim(0,0.5)+
            geom_point(size=1.5, shape=21)+theme_bw()+
            scale_color_manual(values=paste0(cols[c(1,5)],"B3"))+
        ggsave(paste0("Output/MF_PID/PIDvsPIDcon/",fname,".overlap.pdf"), width = 7,height = 2.2)

    }








DF<-Overview[[1]]
DF$Group<-"PIDconsensus"
DF2<-Overview[[2]]
DF2$Group<-"PID"

DF<-rbind(DF,DF2)
#Plot minor variant freq (Ts)
ggplot(data=DF, aes(x=pos, y=freq.Ts, color=Group))+
    ylab("MVF")+xlab("")+ylim(0,0.5)+
    geom_point(size=0.8)+theme_bw()+
    scale_color_manual(values=paste0(cols[c(1,5)],"B3"))


ggplot(data=DF, aes(x=pos, y=freq.mutations, color=Group))+
    ylab("MVF")+xlab("")+ylim(0,0.5)+
    geom_point(size=0.8)+theme_bw()+
    scale_color_manual(values=paste0(cols[c(1,5)],"B3"))

mean(DF$freq.mutations[DF$Group=="PIDconsensus"], na.rm = T)
# 0.03337619
mean(DF$freq.mutations[DF$Group=="PID"], na.rm = T)
# 0.03340586

mean(DF$freq.transv[DF$Group=="PIDconsensus"], na.rm = T)
# 0.02234762
mean(DF$freq.transv[DF$Group=="PID"], na.rm = T)
# 0.0223701
    
Plot<-list()
for (j in 1:length(Overview)){
        df<-Overview[[j]]
        Plot[[j]]<-ggplot(data=df, aes(x=pos, y=freq.Ts))+
            ylab("MVF")+xlab("")+ylim(0,0.5)+
            geom_point(size=0.7, color=cols[5])+theme_bw()+
            ggtitle(names(Overview)[j])+theme(plot.title = element_text(size=11))
}
pdf(paste0("Output/MF_PID/Run2_10_comapre.pdf"), width = 7, height = 5)
do.call(grid.arrange, c(Plot, ncol=1))
dev.off()

for (j in 1:length(Overview)){
    df<-Overview[[j]]
    Plot[[j]]<-ggplot(data=df, aes(x=pos, y=freq.Ts))+
        ylab("MVF")+xlab("")+ylim(0,0.5)+
        geom_point(size=0.7, color=cols[5])+theme_bw()+
        ggtitle(names(Overview)[j])+theme(plot.title = element_text(size=11))
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
    pdf(paste0("Output/MF_PID/",monkey,".all.pdf"), width = 7, height = length(ovDF)*2)
    do.call(grid.arrange, c(Plot, ncol=1))
    dev.off()
}


