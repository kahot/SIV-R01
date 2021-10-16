library(ggplot2)
library(reshape2)
#library(ggpubr)
#library(ggthemes)
#library(plotrix)
#library(grid)
library(gridExtra)
library(colorspace)
library(cowplot)
library(DataCombine)
library(dplyr)
source("Rscripts/baseRscript2.R")
cols2<-qualitative_hcl(6, palette="Dark3")

# read the files saved in Overview_output:
SIVFiles_overview<-list.files("Output/OverviewF_PID/",pattern="filtered.overview.csv")

Overview<-list()
for (i in 1:length(SIVFiles_overview)){ 
        overviews<-read.csv(paste0("Output/OverviewF_PID/",SIVFiles_overview[i]),stringsAsFactors=F, row.names = 1)
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
ggsave("Output/Stock_mf_PID.pdf", width = 7, height = 2)

s2<-stock[stock$freq.mutations.ref>=0.01,]
s2=s2[!is.na(s2$pos),]  #76 positions over 1%
s2<-stock[stock$freq.mutations.ref>=0.005,]
s2=s2[!is.na(s2$pos),]  #99 positions over 0.5%
s2<-stock[stock$freq.mutations.ref>=0.05,]
s2=s2[!is.na(s2$pos),]  #24 positions over 5%

s2<-stock[stock$freq.mutations.ref>=0.1,]
s2=s2[!is.na(s2$pos),]  #19 positions over 10%

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
# 0.01144362
#Syn mean
mean(df3$freq,na.rm = T)
# 0.003735265

ave<-format(round(mean(stock$freq.mutations.ref, na.rm=T)*100,2),nsmall = 2)

Ps2<-ggplot(data=stock, aes(x=pos, y=freq.Ts.ref))+
    ylab("Mutation frequency")+xlab("")+ylim(0,1)+
    geom_point(size=0.7, color=cols2[4])+theme_bw()+
    ggtitle("Stock")+theme(plot.title = element_text(size=12))+
    annotate(geom='text', x=600, y=0.96, label=paste0("MF: ", ave,"%"),color ='gray30', size=3,hjust=0)
    


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
###2)  include transversion

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
    msites<-unique(mutSites$pos)
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
    pdf(paste0("Output/MF_PID/All.",monkey,".stock.pdf"), width = 7, height = length(ovDF)*2)
    do.call(grid.arrange, c(Plot, ncol=1))
    dev.off()
    
    tbweek<-tbs$tb[tbs$ids==monkey]
    ggplot(data=Mm, aes(x=Week, y=MF))+
        ylab("Mutation frequency")+xlab("Weeks")+
        facet_wrap(~ Position, nrow=rown, ncol=5)+
        geom_point(size=2, color=cols2[1])+theme_bw()+ggtitle(paste(monkey))+
        geom_vline(xintercept=tbweek, col="blue")+
        ggsave(paste0("Output/Timeseries_PID/", monkey, ".all.pdf"), heigh=rown*2, width =10)
}



## Plot with averages
for (i in 1:length(monkeys2)){
    print(ids[i])
    sample<-monkeys2[[i]]
    sample = sample[order(sample[,'Week'],-sample[,'Run']),]
    sample = sample[!duplicated(sample$Week),]
    ovDF<-Overview[as.vector(sample$File.name)]
    monkey<-names(monkeys2)[i]
    Plot<-list()
    Plot[[1]]<-Ps2
    for (j in 1:length(ovDF)){
        df<-ovDF[[j]]
        week<-sample$Week[j]
        ave<-format(round(mean(df$freq.mutations.ref, na.rm=T)*100,2),nsmall = 2)
        div<-format(round(mean(df$freq.mutations, na.rm=T)*100,2),nsmall = 2)
        
        if (week<=17) col=cols2[3]
        else {col=cols2[1]}
        Plot[[(j+1)]]<-ggplot(data=df, aes(x=pos, y=freq.Ts.ref))+
            ylab("Mutation frequency")+xlab("")+ylim(0,1)+
            geom_point(size=0.7, color=col)+theme_bw()+
            ggtitle(paste0(monkey," Week ",week))+theme(plot.title = element_text(size=12))+
            annotate(geom='text', x=600, y=0.96, label=paste0("MF: ", ave,"%"),color ='gray30', size=3,hjust=0)+
            annotate(geom='text', x=600, y=0.82, label=paste0("Diversity: ", div,"%"),color ='gray30', size=3,hjust=0)
    }
    pdf(paste0("Output/MF_PID/Averages.",monkey,".all.pdf"), width = 7, height = length(ovDF)*2)
    do.call(grid.arrange, c(Plot, ncol=1))
    dev.off()
}





#Find the overlapping positions
poss<-lapply(PosType, function(x) x$mutation<-paste0("pos.",x$pos,".",x$type))
sites<-unique(unlist(poss))
sites<-sites[order(sites)]

positions<-unique(unlist(unname(Pos)))
positions<-positions[order(positions)]
cpos<-data.frame(Pos=positions)
for (i in 1:length(monkeys2)){
    vec<-Pos[[i]]
    cpos[,names(monkeys2)[i]]<-apply(cpos['Pos'],2, function(x) ifelse(x %in% vec, "Y", "N"))
}

cpos

#write.csv(cpos, "Output/HighMutfreq_overtime.csv")
cpos$codon<-apply(cpos["Pos"], 1, function(x) if(x%%3==0) x=3 else x=x%%3)


df<-read.csv("Output/Overview.ref.csv", stringsAsFactors = F, row.names = 1)
df<-df[df$pos %in% positions,]
colnames(cpos)[1]<-"pos"
cpos<-merge(cpos, df, by="pos")

write.csv(cpos, "Output/HighMutfreq_sites_all.csv")


####### 
#Stock position high (how changed?)
sites<-stock$pos[stock$freq.Ts.ref>0.1]
sites<-unique(sites)
sites<-sites[!is.na(sites)]
sites<-sites[order(sites)]
df1<-data.frame(pos=sites, type="Ts")

sites1<-stock$pos[stock$freq.transv1.ref>0.1]
sites1<-unique(sites1)
sites1<-sites1[!is.na(sites1)]
sites1<- sites1[order(sites1)]
df2<-data.frame(pos=sites1, type="Tv1")

sites2<-stock$pos[stock$freq.transv2.ref>0.1]
sites2<-unique(sites2)
sites2<-sites2[!is.na(sites2)]
sites2<-sites2[order(sites2)]
df3<-data.frame(pos=sites2, type="Tv2")

mutSites<-rbind(df1,rbind(df2,df3))
msites<-unique(c(sites,sites1,sites2))
msites<-msites[order(msites)]

#select only the positions in msites
mfDF1<-stock[stock$pos %in% sites,]
mfDF2<-stock[stock$pos %in% sites1,]
mfDF3<-stock[stock$pos %in% sites2,]

#Transition
Ts<-mfDF1[,c("pos","ref","freq.Ts.ref")]
colnames(Ts)[3]<-"freq"
Ts$Type<-"Ts"
Ts$mutation<-paste0("pos.",Ts$pos,'.',"Ts")
#Tranv1
Tv1<-mfDF2[,c("pos","ref","freq.transv1.ref")]
colnames(Tv1)[3]<-"freq"
Tv1$Type<-"Tv1"
Tv1$mutation<-paste0("pos.",Tv1$pos,'.',"Tv1")

#Tranv2
Tv2<-mfDF3[,c("pos","ref","freq.transv2.ref")]
colnames(Tv2)[3]<-"freq"
Tv2$Type<-"Tv2"
Tv2$mutation<-paste0("pos.",Tv2$pos,'.',"Tv2")

stockHigh<-rbind(Ts, Tv1,Tv2)


Highs<-list()
for (i in 1:length(monkeys2)){
    sample<-monkeys2[[i]]
    sample = sample[order(sample[,'Week'],-sample[,'Run']),]
    sample = sample[!duplicated(sample$Week),]
    ovDF<-Overview[as.vector(sample$File.name)]
    monkey<-names(monkeys2)[i]
    
    #Transition
    hTs<-lapply(ovDF, function(x) x[x$pos %in% Ts$pos,c("pos","ref","freq.Ts.ref")])
    tTs<-unname(sapply(hTs, `[[`, "freq.Ts.ref"))
    tTs<-data.frame(t(tTs))
    nt1<-hTs[[1]][,"ref"]
    sites<-hTs[[1]]["pos"]
    colnames(tTs)<-paste0("Pos ",sites$pos,' Ts')
    tTs<-InsertRow(tTs, Ts$freq, RowNum=1)
    tTs$Week<-c(0,sample$Week)
    
    #Tsm<-melt(tTs, id.vars = "Week")
    #colnames(Tsm)[2:3]<-c("Position", "MF")
    #rown<-ceiling(length(sites)/5)
    
    #Tranv1
    hTv1<-lapply(ovDF, function(x) x[x$pos %in% Tv1$pos,c("pos","ref","freq.transv1.ref")])
    tTv1<-unname(sapply(hTv1, `[[`, "freq.transv1.ref"))
    tTv1<-data.frame(t(tTv1))
    nt1<-hTv1[[1]][,"ref"]
    sites<-hTv1[[1]]["pos"]
    colnames(tTv1)<-paste0("Pos ",sites$pos,' Tv1')
    tTv1<-InsertRow(tTv1, Tv1$freq, RowNum=1)
    
    #Tranv2
    hTv2<-lapply(ovDF, function(x) x[x$pos %in% Tv2$pos,c("pos","ref","freq.transv2.ref")])
    tTv2<-unname(sapply(hTv2, `[[`, "freq.transv2.ref"))
    tTv2<-data.frame(t(tTv2))
    nt1<-hTv2[[1]][,"ref"]
    sites<-hTv2[[1]]["pos"]
    colnames(tTv2)<-paste0("Pos ",sites$pos,' Tv2')
    tTv2<-InsertRow(tTv2, Tv2$freq, RowNum=1)
    
    
    all<-cbind(tTs, tTv1,tTv2)
    
    colorder<-order(colnames(all))
    all<-all[,colorder]
    mutations<-colnames(all)[1:20]
    allm<-melt(all, id.vars = "Week")
    colnames(allm)[2:3]<-c("Position", "MF")
    rown<-ceiling((ncol(all)-1)/5)
    allm$Monkey<-monkey
    
    Highs[[i]]<-allm
    names(Highs)[i]<-monkey
    
    tbweek<-tbs$tb[tbs$ids==monkey]
    #ggplot(data=allm, aes(x=Week, y=MF))+
    #    ylab("Mutation frequency")+xlab("Weeks")+
    #    facet_wrap(~ Position, nrow=rown, ncol=5)+
    #    geom_point(size=2, color=cols2[1])+theme_bw()+ggtitle(paste(monkey))+
    #    geom_vline(xintercept=tbweek, col="blue")+
    #    ggsave(paste0("Output/Timeseries_PID/StockHigh.", monkey, ".all.pdf"), heigh=rown*2, width =10)
}    


colnames(tbs)[1]<-"Monkey"

for (i in 1: length(mutations)){
    mu<-mutations[i]
    Mu<-lapply(Highs, function(x) x<-x[x$Position==mu,])
    Mu2<-bind_rows(Mu)
    ggplot(data=Mu2, aes(x=Week, y=MF))+
        ylab("Mutation frequency")+xlab("")+
        geom_point(size=2, color=cols2[1])+theme_bw()+ggtitle(paste(Mu2$Position[1]))+
        geom_vline(data = tbs, mapping = aes(xintercept = tb), col="blue")+
        facet_wrap(~ Monkey, nrow=7, ncol=1)+
        ggsave(paste0("Output/Timeseries_PID/", Mu2$Position[1], ".pdf"), heigh=8, width =3)
        

}
