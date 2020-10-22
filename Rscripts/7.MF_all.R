library(ggplot2)
library(reshape)
library(ggpubr)
library(ggthemes)
library(plotrix)
library(grid)
library(gridExtra)
library(colorspace)
library(cowplot)
library(DataCombine)
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



######   #Plot the stock virus ######
# Stock is Run0_17
stock<-Overview[[which(names(Overview)=="Run0_17")]]
Ps<-ggplot(data=stock, aes(x=pos, y=freq.Ts.ref))+
    ylab("Mutation frequency")+xlab("")+ylim(0,1)+
    geom_point(size=0.7, color=cols2[4])+theme_bw()+
    ggtitle("Stock")+theme(plot.title = element_text(size=12))
ggsave("Output/Stock_mf.pdf", width = 7, height = 2)

s2<-stock[stock$freq.mutations.ref>=0.1,]



## Select th ones that are interesting

ids<-c("A21918","A22317","A22517", "A22617","A22217","A23918","A22117")
#Tb infection week
tb<-c(17.5,20,20,20,17, 17.5, 19)
tbs<-data.frame(cbind(ids,tb))
tbs$tb<-as.numeric(as.character(tbs$tb))


samples<-read.csv("Data/SamplesNoduplicates.csv",stringsAsFactors = F)
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
monkeys2<-monkeyList[ids]

#1. Transition only
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
    Ts<-InsertRow(Ts, c(stock$freq.Ts.ref[stock$pos %in% sites],0), RowNum=1)
    
    Tsm<-melt(Ts, id.vars = "Week")
    colnames(Tsm)[2:3]<-c("Position", "MF")
    rown<-ceiling(length(sites)/5)
    
    Plot<-list()
    Plot[[1]]<-Ps
    for (j in 1:length(ovDF)){
        df<-ovDF[[j]]
        week<-sample$Week[j]
        if (week<=17) col=cols2[3]
        else {col=cols2[1]}
        Plot[[(j+1)]]<-ggplot(data=df, aes(x=pos, y=freq.Ts.ref))+
            ylab("Mutation frequency")+xlab("")+ylim(0,1)+
            geom_point(size=0.7, color=col)+theme_bw()+
            ggtitle(paste0(monkey," Week ",week))+theme(plot.title = element_text(size=12))
    }
    #pdf(paste0("Output/MF/",monkey,".stockp.pdf"), width = 7, height = length(ovDF)*2)
    #do.call(grid.arrange, c(Plot, ncol=1))
    #dev.off()
    
    tbweek<-tbs$tb[tbs$ids==monkey]
    ggplot(data=Tsm, aes(x=Week, y=MF))+
        ylab("Mutation frequency")+xlab("Weeks")+
        facet_wrap(~ Position, nrow=rown, ncol=5)+
        geom_point(size=2, color=cols2[1])+theme_bw()+ggtitle(paste(monkey))+
        geom_vline(xintercept=tbweek, col="blue")
    #+
    #    ggsave(paste0("Output/Timeseries2/", monkey, ".pdf"), heigh=rown*2, width =10)
}

#Find the overlapping positions
positions<-unique(unlist(unname(Pos)))
positions<-positions[order(positions)]
cpos<-data.frame(Pos=positions)
for (i in 1:length(monkeys2)){
    vec<-Pos[[i]]

    cpos[,names(monkeys2)[i]]<-apply(cpos['Pos'],2, function(x) ifelse(x %in% vec, "Y", "N"))
}

#write.csv(cpos, "Output/HighMutfreq_overtime.csv")
cpos$codon<-apply(cpos["Pos"], 1, function(x) if(x%%3==0) x=3 else x=x%%3)


df<-read.csv("Output/Overview.ref.csv", stringsAsFactors = F, row.names = 1)
df<-df[df$pos %in% positions,]
colnames(cpos)[1]<-"pos"
cpos<-merge(cpos, df, by="pos")

write.csv(cpos, "Output/HighMutfreq_sites.csv")


######################
### include transversion

Pos<-list()
PosType<-list()
for (i in 1:length(monkeys2)){
    print(ids[i])
    sample<-monkeys2[[i]]
    sample = sample[order(sample[,'Week'],-sample[,'Run']),]
    sample = sample[!duplicated(sample$Week),]
    ovDF<-Overview[as.vector(sample$File.name)]
    monkey<-names(monkeys2)[i]
    
    sites<-lapply(ovDF, function(x) x$pos[x$freq.Ts.ref>0.1])
    sites<-unique(unlist(sites))
    sites<-sites[!is.na(sites)]
    sites<-sites[order(sites)]
    df1<-data.frame(pos=sites, type="Ts")
    
    sites1<-lapply(ovDF, function(x) x$pos[x$freq.transv1.ref>0.1])
    sites1<-unique(unlist(sites1))
    sites1<-sites1[!is.na(sites1)]
    sites1<- sites1[order(sites1)]
    df2<-data.frame(pos=sites1, type="Tv1")
    
    sites2<-lapply(ovDF, function(x) x$pos[x$freq.transv2.ref>0.1])
    sites2<-unique(unlist(sites2))
    sites2<-sites2[!is.na(sites2)]
    sites2<-sites2[order(sites2)]
    df3<-data.frame(pos=sites2, type="Tv2")
    
    mutSites<-rbind(df1,rbind(df2,df3))
    msites<-unique(c(sites,sites1,sites2))
    msites<-msites[order(msites)]
    
    
    #select only the positions in msites
    mfDF1<-lapply(ovDF, function(x) x<-x[x$pos %in% sites,])
    mfDF2<-lapply(ovDF, function(x) x<-x[x$pos %in% sites1,])
    mfDF3<-lapply(ovDF, function(x) x<-x[x$pos %in% sites2,])
    
    #Transition
    Ts<-unname(sapply(mfDF1, `[[`, "freq.Ts.ref"))
    Ts<-data.frame(t(Ts))
    nt1<-mfDF1[[1]][,"ref"]
    colnames(Ts)<-paste0("Pos_",sites, ' ',nt1, ' Ts')
    Ts<-InsertRow(Ts, c(stock$freq.Ts.ref[stock$pos %in% sites]), RowNum=1)
    
    #Tranv1
    Tv1<-unname(sapply(mfDF2, `[[`, "freq.transv1.ref"))
    if (nrow(mfDF2[[1]])==1) Tv1<-data.frame(Tv1)
    else { 
    Tv1<-data.frame(t(Tv1)) 
    }
    nt2<-mfDF2[[1]][,"ref"]
    colnames(Tv1)<-paste0("Pos_",sites1, ' ',nt2, ' Tv1')
    Tv1<-InsertRow(Tv1, c(stock$freq.transv1.ref[stock$pos %in% sites1]), RowNum=1)
    
    #Tranv2
    Tv2<-unname(sapply(mfDF3, `[[`, "freq.transv2.ref"))
    if (nrow(mfDF3[[1]])==1) Tv2<-data.frame(Tv2)
    else { 
        Tv2<-data.frame(t(Tv2))
        }
    
    nt3<-mfDF3[[1]][,"ref"]
    colnames(Tv2)<-paste0("Pos_",sites2, ' ',nt3, ' Tv2')
    Tv2$Week<-sample$Week
    Tv2<-InsertRow(Tv2, c(stock$freq.transv2.ref[stock$pos %in% sites2],0), RowNum=1)
    
    M<-cbind(Ts, Tv1,Tv2)
    
    colorder<-order(colnames(M))
    M<-M[,colorder]
    
    PosType[[i]]<-mutSites
    names(PosType)[i]<-monkey
    Pos[[i]]<-msites
    names(Pos)[i]<-monkey
               
       
    
    Mm<-melt(M, id.vars = "Week")
    colnames(Mm)[2:3]<-c("Position", "MF")
    rown<-ceiling((ncol(M)-1)/5)
    
    Plot<-list()
    Plot[[1]]<-Ps
    for (j in 1:length(ovDF)){
        df<-ovDF[[j]]
        week<-sample$Week[j]
        if (week<=17) col=cols2[3]
        else {col=cols2[1]}
        Plot[[(j+1)]]<-ggplot(data=df, aes(x=pos, y=freq.Ts.ref))+
            ylab("Mutation frequency")+xlab("")+ylim(0,1)+
            geom_point(size=0.7, color=col)+theme_bw()+
            ggtitle(paste0(monkey," Week ",week))+theme(plot.title = element_text(size=12))
    
    }
    pdf(paste0("Output/MF/All.",monkey,".stock.pdf"), width = 7, height = length(ovDF)*2)
    do.call(grid.arrange, c(Plot, ncol=1))
    dev.off()
    
    tbweek<-tbs$tb[tbs$ids==monkey]
    ggplot(data=Mm, aes(x=Week, y=MF))+
        ylab("Mutation frequency")+xlab("Weeks")+
        facet_wrap(~ Position, nrow=rown, ncol=5)+
        geom_point(size=2, color=cols2[1])+theme_bw()+ggtitle(paste(monkey))+
        geom_vline(xintercept=tbweek, col="blue")+
        ggsave(paste0("Output/Timeseries2/", monkey, ".all.pdf"), heigh=rown*2, width =10)
}

#Find the over lapping positions
positions<-unique(unlist(unname(Pos)))
positions<-positions[order(positions)]
cpos<-data.frame(Pos=positions)
for (i in 1:length(monkeys2)){
    vec<-Pos[[i]]
    
    cpos[,names(monkeys2)[i]]<-apply(cpos['Pos'],2, function(x) ifelse(x %in% vec, "Y", "N"))
}

#write.csv(cpos, "Output/HighMutfreq_overtime.csv")
cpos$codon<-apply(cpos["Pos"], 1, function(x) if(x%%3==0) x=3 else x=x%%3)


df<-read.csv("Output/Overview.ref.csv", stringsAsFactors = F, row.names = 1)
df<-df[df$pos %in% positions,]
colnames(cpos)[1]<-"pos"
cpos<-merge(cpos, df, by="pos")

write.csv(cpos, "Output/HighMutfreq_sites_all.csv")



