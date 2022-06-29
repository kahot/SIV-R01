library(reshape2)
library(colorspace)
library(DataCombine)
source("Rscripts/baseRscript.R")
cols2<-qualitative_hcl(7, palette="Dark3")
hcl_palettes("qualitative", n=7,plot = T)

SIVFiles_overview<-list.files("Output/OverviewF_PIDcon/",pattern=".csv")

#Files with low read deapth -> use unfiltered overview files
SIVfiles<-list.files("Output/Overview_PIDcon/", pattern="Run5|Run6|Run7")
SIVfiles<-SIVfiles[SIVfiles!="Run5_12_overview.csv"]
SIVfiles<-SIVfiles[SIVfiles!="Run5_13_overview.csv"]
SIVfiles<-SIVfiles[SIVfiles!="Run5_15_overview.csv"]
SIVfiles<-c(SIVfiles,"Run0_14_overview.csv")
#Remove Run6_18 (has only 16 seq) for now
SIVfiles<-SIVfiles[SIVfiles!="Run6_18_overview.csv"]

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


# Stock is Run0_17
stock<-Overview[[which(names(Overview)=="Run0_17")]]
#TB infection week
tbs<-read.csv("Data/TBinfection.week.csv",stringsAsFactors = F)

samples<-read.csv("Data/SamplesNoduplicates.csv",stringsAsFactors = F)
samples<-samples[samples$File.name!="Run6_18",]
samples$Week<-as.integer(samples$Week)
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
monkeys2<-monkeyList[tbs$ids]


#
#Look for sites with mutation freq (diversity) >10% 

####################################
#Pos<-list()
#PosType<-list()
for (i in 1:length(monkeys2)){
    monkey<-names(monkeys2)[i]
    print(monkey)
    sample<-monkeys2[[i]]
    sample = sample[order(sample[,'Week']),]
    ovDF<-Overview[as.vector(sample$File.name)]
    
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
    ##High freq sites vector
    #msites<-unique(c(sites,sites1,sites2))
    #msites<-msites[order(msites)]
    
    #select only the positions in msites
    mfDF<-lapply(ovDF, function(x) x<-x[x$pos %in% sites,])
    mfDF1<-lapply(ovDF, function(x) x<-x[x$pos %in% sites1,])
    mfDF2<-lapply(ovDF, function(x) x<-x[x$pos %in% sites2,])
    
    #Transition
    Ts<-unname(sapply(mfDF, `[[`, "freq.Ts.ref"))
    Ts<-data.frame(t(Ts))
    nt1<-mfDF[[1]][,"ref"]
    colnames(Ts)<-paste0("Pos_",sites, ' ',nt1, ' Ts')
    Ts<-InsertRow(Ts, c(stock$freq.Ts.ref[stock$pos %in% sites]), RowNum=1)
    Ts<-InsertRow(Ts, c(stock$freq.Ts.ref[stock$pos %in% sites]), RowNum=1)
    
    #Tranv1
    Tv1<-unname(sapply(mfDF1, `[[`, "freq.transv1.ref"))
    if (nrow(mfDF1[[1]])==1) Tv1<-data.frame(Tv1)
    if (nrow(mfDF1[[1]])>1) Tv1<-data.frame(t(Tv1)) 
    
    nt2<-mfDF1[[1]][,"ref"]
    colnames(Tv1)<-paste0("Pos_",sites1, ' ',nt2, ' Tv1')
    Tv1<-InsertRow(Tv1, c(stock$freq.transv1.ref[stock$pos %in% sites1]), RowNum=1)
    Tv1<-InsertRow(Tv1, c(stock$freq.transv1.ref[stock$pos %in% sites1]), RowNum=1)
    
    #Tranv2
    Tv2<-unname(sapply(mfDF2, `[[`, "freq.transv2.ref"))
    if (nrow(mfDF2[[1]])==1) Tv2<-data.frame(Tv2)
    if (nrow(mfDF2[[1]])>1) Tv2<-data.frame(t(Tv2))
    
    nt3<-mfDF2[[1]][,"ref"]
    colnames(Tv2)<-paste0("Pos_",sites2, ' ',nt3, ' Tv2')
    Tv2$Week<-sample$Week
    Tv2<-InsertRow(Tv2, c(stock$freq.transv2.ref[stock$pos %in% sites2],0), RowNum=1)
    Tv2<-InsertRow(Tv2, c(stock$freq.transv2.ref[stock$pos %in% sites2],0), RowNum=1)
    Tv2$Tissue<-c("Plasma","Tissue", sample$Tissue)
    
    M<-cbind(Ts, Tv1,Tv2)
    
    colorder<-order(colnames(M))
    M<-M[,colorder]
    
    #PosType[[i]]<-mutSites
    #names(PosType)[i]<-monkey
    #Pos[[i]]<-msites
    #names(Pos)[i]<-monkey
    
    M$Tissue<-ifelse(M$Tissue=="Plasma", "Plasma","Tissue")
    
    Mm<-melt(M, id.vars = c("Week", "Tissue"))
    colnames(Mm)[3:4]<-c("Position", "MF")
    rown<-ceiling((ncol(M)-1)/5)
    
    tbweek<-tbs$tb[tbs$ids==monkey]
    ggplot(data=Mm, aes(x=Week, y=MF, color=Tissue))+
        ylab("Mutation frequency")+xlab("Weeks")+
        facet_wrap(~ Position, nrow=rown, ncol=5)+
        geom_point(size=2)+theme_bw()+ggtitle(paste(monkey))+
        geom_line(stat="identity")+
        geom_vline(xintercept=tbweek, col="blue")
    ggsave(paste0("Output/Timeseries_PID/", monkey, ".Both.pdf"), heigh=rown*1.5, width =11.5)
}


#Sites of interest
diffpos<-c(475,494,500,515,518,578,648)
summary<-read.csv("Output/MF_PID/HighFreq/Summary_highFreq_sites.csv", row.names = 1)
summary2<-summary[summary$pos %in% diffpos,]


monkey3<-monkeys[c(1:6,10:11)]
coltis<-c("gray",cols2[c(1,3,2,4,6,7,5)])    
Frq<-list()

#Look at the primary mutation     
for (i in 1:length(monkey3)){
    monkey<-monkey3[i]
    print(monkey)
    sample<-monkeys2[[monkey]]
    sample = sample[order(sample[,'Week']),]
    ovDF<-Overview[as.vector(sample$File.name)]
    
    
    ts<-summary2$pos[summary2$Type.1=="Ts"]  
    tv1<-summary2$pos[summary2$Type.1=="Tv1"]  
    tv2<-summary2$pos[summary2$Type.1=="Tv2"]  
    
    #select only the positions in msites
    mfDF<-lapply(ovDF, function(x) x<-x[x$pos %in% ts,])
    mfDF1<-lapply(ovDF, function(x) x<-x[x$pos %in% tv1,])
    mfDF2<-lapply(ovDF, function(x) x<-x[x$pos %in% tv2,])
    
    #Transition
    Ts<-unname(sapply(mfDF, `[[`, "freq.Ts.ref"))
    Ts<-data.frame(t(Ts))
    nt1<-mfDF[[1]][,"ref"]
    colnames(Ts)<-paste0("Pos_",ts, ' ',nt1, ' Ts')
    Ts<-InsertRow(Ts, c(stock$freq.Ts.ref[stock$pos %in% ts]), RowNum=1)
    
    #Tranv1
    Tv1<-unname(sapply(mfDF1, `[[`, "freq.transv1.ref"))
    if (nrow(mfDF1[[1]])==1) Tv1<-data.frame(Tv1)
    if (nrow(mfDF1[[1]])>1) Tv1<-data.frame(t(Tv1)) 
    
    nt2<-mfDF1[[1]][,"ref"]
    colnames(Tv1)<-paste0("Pos_",tv1, ' ',nt2, ' Tv1')
    Tv1<-InsertRow(Tv1, c(stock$freq.transv1.ref[stock$pos %in% tv1]), RowNum=1)
    
    #Tranv2
    Tv2<-unname(sapply(mfDF2, `[[`, "freq.transv2.ref"))
    if (nrow(mfDF2[[1]])==1) Tv2<-data.frame(Tv2)
    if (nrow(mfDF2[[1]])>1) Tv2<-data.frame(t(Tv2))
    
    nt3<-mfDF2[[1]][,"ref"]
    colnames(Tv2)<-paste0("Pos_",tv2, ' ',nt3, ' Tv2')
    Tv2$Week<-sample$Week
    Tv2<-InsertRow(Tv2, c(stock$freq.transv2.ref[stock$pos %in% tv2],0), RowNum=1)
    Tv2$Tissue<-c("Plasma", sample$Tissue)
    Tv2$Type<-c("Plasma", sample$Tissue.type)
    
    M<-cbind(Ts, Tv1,Tv2)
    colorder<-order(colnames(M))
    M<-M[,colorder]
    #save as one file
    m2<-M
    m2$Monkey<-monkey
    m2$Run<-c(0,sample$Run)
    m2$File.name<-c("Run0_16",sample$File.name)
    
    Mut<-rbind(Mut, m2)
    
    Mm<-melt(M, id.vars = c("Week", "Tissue","Type"))
    colnames(Mm)[4:5]<-c("Position", "MF")
    rown<-ceiling((ncol(M)-1)/5)
    
    tbweek<-tbs$tb[tbs$ids==monkey]
    types<-as.character(unique(Mm$Type))
    types<-types[!(types=="Plasma"|types=="dLN"|types=="pLN"|types=="LN")]
    Mm$Type<-factor(Mm$Type, levels = c("Plasma","dLN","pLN","LN",types))
    Mm<-Mm[!is.na(Mm$MF),]
    Pm<-Mm[Mm$Type=="Plasma",]
    
    l<-nrow(sample[sample$Tissue!="Plasma"& sample$Week==max(sample$Week),])
    ggplot()+
        ylab("Mutation frequency")+xlab("Week")+
        scale_color_manual(values=coltis)+
        facet_wrap(~ Position, nrow=rown, ncol=3)+
        geom_line(data=Pm, aes(x=Week, y=MF),color="gray")+
        geom_point(data=Mm, aes(x=Week, y=MF, color=Type),size=2.5, position=position_jitter(width = 0.6), alpha=0.7)+theme_bw()+
        ggtitle(paste0(monkey, " (# tissues at necropsy=", l,")"))+
        geom_vline(xintercept=tbweek, col="blue")
    #ggsave(paste0("Output/MF_PID/HighFreq/Tissue.", monkey, ".pdf"), height=rown*1.8, width =7)

    Mm$Monkey<-monkey
    Frq[[i]]<-Mm
}

write.csv(Mut, "Output/QC/SixNTfreq.csv")

poss<-unique(Frq[[1]]['Position'])
poss<-paste(poss$Position)
for (i in 1:length(poss)){
    position<-poss[i]
    P<-lapply(Frq, function(x) x[x$Position==position,])
    Mm<-do.call(rbind, P)
    ggplot()+
        ylab("Mutation frequency")+xlab("Weeks")+
        scale_color_manual(values=coltis)+
        facet_wrap(~ Monkey, ncol=4)+
        geom_line(data=Mm[Mm$Type=="Plasma",], aes(x=Week, y=MF),color="gray")+
        geom_point(data=Mm, aes(x=Week, y=MF, color=Type),size=2.5, position=position_jitter(width = 0.7), alpha=0.7)+theme_bw()+
        ggtitle(paste0(position))
    posi<-substr(position, 1,7)
    ggsave(paste0("Output/MF_PID/HighFreq/Tissue.",posi, ".pdf"), heigh=rown*1.8, width =7.5)
    
}


#Use AAfreq
Run6<-read.csv("Output/QC/")


