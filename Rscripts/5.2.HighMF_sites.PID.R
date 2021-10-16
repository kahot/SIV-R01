library(reshape2)
library(gridExtra)
library(colorspace)
library(DataCombine)
source("Rscripts/baseRscript.R")
cols2<-qualitative_hcl(6, palette="Dark3")

# read the files saved in Overview_output:
#SIVFiles_overview<-list.files("Output/OverviewF_PID/",pattern="filtered.overview.csv")
#Overview<-list()
#for (i in 1:length(SIVFiles_overview)){ 
#        overviews<-read.csv(paste0("Output/OverviewF_PID/",SIVFiles_overview[i]),stringsAsFactors=F, row.names = 1)
#        Overview[[i]]<-overviews
#        names(Overview)[i]<-substr(paste(SIVFiles_overview[i]),start=1,stop=7)
#}

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


# Stock is Run0_17
stock<-Overview[[which(names(Overview)=="Run0_17")]]
#TB infection week
tbs<-read.csv("Data/TBinfection.week.csv",stringsAsFactors = F)

samples<-read.csv("Data/SamplesNoduplicates.csv",stringsAsFactors = F)
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


######################
#Look for sites with mutation freq (diversity) >10% 
Pos<-list()
PosType<-list()
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
    msites<-unique(c(sites,sites1,sites2))
    msites<-msites[order(msites)]
    
    
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
    
    #Tranv1
    Tv1<-unname(sapply(mfDF1, `[[`, "freq.transv1.ref"))
    if (nrow(mfDF1[[1]])==1) Tv1<-data.frame(Tv1)
    else { 
    Tv1<-data.frame(t(Tv1)) 
    }
    nt2<-mfDF1[[1]][,"ref"]
    colnames(Tv1)<-paste0("Pos_",sites1, ' ',nt2, ' Tv1')
    Tv1<-InsertRow(Tv1, c(stock$freq.transv1.ref[stock$pos %in% sites1]), RowNum=1)
    
    #Tranv2
    Tv2<-unname(sapply(mfDF2, `[[`, "freq.transv2.ref"))
    if (nrow(mfDF2[[1]])==1) Tv2<-data.frame(Tv2)
    else { 
        Tv2<-data.frame(t(Tv2))
        }
    
    nt3<-mfDF2[[1]][,"ref"]
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
    
    tbweek<-tbs$tb[tbs$ids==monkey]
    #ggplot(data=Mm, aes(x=Week, y=MF))+
    #    ylab("Mutation frequency")+xlab("Weeks")+
    #    facet_wrap(~ Position, nrow=rown, ncol=5)+
    #    geom_point(size=2, color=cols2[1])+theme_bw()+ggtitle(paste(monkey))+
    #    geom_vline(xintercept=tbweek, col="blue")
    #    ggsave(paste0("Output/Timeseries_PID/", monkey, ".con.all.pdf"), heigh=rown*1.5, width =11.5)
}


#Find the over lapping positions
positions<-unique(unlist(unname(Pos)))
positions<-positions[order(positions)]
cpos<-data.frame(Pos=positions)
for (i in 1:length(monkeys2)){
    vec<-Pos[[i]]
    
    cpos[,names(monkeys2)[i]]<-sapply(cpos$Pos,function(x) ifelse(x %in% vec, "Y", "N"))
}

cpos$codon<-apply(cpos["Pos"], 1, function(x) if(x%%3==0) x=3 else x=x%%3)

df<-read.csv("Output/Overview.ref.csv", stringsAsFactors = F, row.names = 1)
df<-df[,-16]
df<-df[df$pos %in% positions,]
colnames(cpos)[1]<-"pos"
cpos<-merge(cpos, df, by="pos")
cpos$occurence<-apply(cpos[,2:11], 1, function(x) length(x[x=="Y"]))

#Remove the codon position 3 and all synonymous
cpos2<-cpos[!(cpos$Type.r=="syn"&cpos$Type.tv1.r=="syn"&cpos$Type.tv2.r=="syn"),]
#68 positions
write.csv(cpos2, "Output/HighMutfreq_sites_all.csv")



#Find the common sites that appear in many samples

high<-cpos2[cpos2$occurence>5,] #22 sites
highpos<-high$pos

df<-stock

for (j in 1:length(highpos)){
    position<-highpos[j]
    #nucleotide at Ref, Ts, tv1, and tv2
    re<-df$ref[df$pos==position]
    Ts<-transition(re)
    tv1<-transv1(re)
    tv2<-transv2(re)
    
    #mut freq of 'position' over time in each monkey
    Freq<-data.frame()
    for (i in 1:length(monkeys2)){
        monkey<-names(monkeys2[i])
        sample<-monkeys2[[i]]
        sample<-rbind(samples[samples$Tissue=="Stock",], sample)
        sample$Tissue[1]<-"Plasma"
        sample = sample[order(sample[,'Week']),]
        ovDF<-Overview[as.vector(sample$File.name)]
        
        poss<-lapply(ovDF, function(x) x<-x[x$pos ==position,])
        
        Pos<-data.frame(Week=sample$Week, Tissue=sample$Tissue)
        Pos[,paste(Ts)]<-unname(sapply(poss, `[[`, "freq.Ts.ref"))
        Pos[,paste(tv1)]<-unname(sapply(poss, `[[`, "freq.transv1.ref"))
        Pos[,paste(tv2)]<-unname(sapply(poss, `[[`, "freq.transv2.ref"))
        Pos[,paste(re)]<-1-(Pos[,3]+Pos[,4]+Pos[5])
        
        Posm<-melt(Pos, id.vars = c("Week","Tissue"))
        colnames(Posm)[3:4]<-c("nuc", "MF")
        Posm$Monkey<-monkey
        Posm$Tbweek<-tbs$tb[tbs$ids==monkey]
        Freq<-rbind(Freq, Posm)
    }
    
    Freq$nuc<-factor(Freq$nuc,levels=c("a","c",'g','t'))
    #type of mutation and codon position
    #which nuc is the second highest in freq?
    
    type<-high$Type.r[high$pos==position]
    type1<-high$Type.tv1.r[high$pos==position]
    type2<-high$Type.tv2.r[high$pos==position]
    
    codon<-high$codon[high$pos==position]
    
    ggplot(data=Freq, aes(x=Week, y=MF, color=nuc))+
        ylab("Mutation frequency")+xlab("Week")+ylim(0,1)+
        facet_wrap(~ Monkey, nrow=5, ncol=2)+
        geom_point(size=1.5)+theme_bw()+
        scale_color_manual(values=cols2[c(1,3,5,2)], label=c("A","C",'G',"T"))+
        geom_path(data=Freq[Freq$Tissue=="Plasma",], aes(x=Week, y=MF))+
        theme(legend.title = element_blank())+
        geom_vline(aes(xintercept=Tbweek), col="blue")+
        ggtitle(paste0("Pos.",position," Codon=",codon," Type=",type," ", type1," ", type2) )+theme(plot.title = element_text(size=10), axis.title.y = element_text(size=8))
    ggsave(paste0("Output/MF_PID/HighFreq/Pos", position, "_overtime.pdf"), width = 8, height=8)
}
   

# High freq position summary
table(cpos2$codon)
# 1  2  3 
#23 23 16  
table(high$codon)
# 1  2  3 
#4 13  5 
table(cpos2$occurence)
# 1  2  3  4  5  6  7  8  9 10 
#20  6  5  6  3  1  4  4  6  7  
table(high$occurence)
#6  7  8  9 10 
#1  4  4  6  7 

table(high$makesCpG.r)
table(high$makesCpG.tvs.r)
table(high$RefAA)
table(high$MUTAA.r)

summary<-data.frame(pos=highpos)
for (j in 1:length(highpos)){
    position<-highpos[j]
    summary$ref[j]<-stock$ref[stock$pos==position]
    summary$WTAA[j]<-stock$WTAA[stock$pos==position]
    
    #monkey<-names(monkeys2[1])
    sample<-monkeys2[[1]]
    sample<-rbind(samples[samples$Tissue=="Stock",], sample)
    sample$Tissue[1]<-"Plasma"
    sample = sample[order(sample[,'Week']),]
    ovDF<-Overview[as.vector(sample$File.name)]
    
    poss<-lapply(ovDF, function(x) x<-x[x$pos ==position,])
    
    Pos<-data.frame(Week=sample$Week, Tissue=sample$Tissue)
    Pos$Ts<-unname(sapply(poss, `[[`, "freq.Ts.ref"))
    Pos$Tv1<-unname(sapply(poss, `[[`, "freq.transv1.ref"))
    Pos$Tv2<-unname(sapply(poss, `[[`, "freq.transv2.ref"))
    
    ave<-data.frame(colMeans(Pos[c("Ts","Tv1","Tv2")], na.rm=T))
    colnames(ave)[1]<-"freq"
    ave$Type<-rownames(ave)
    ave<-ave[order(ave$freq,decreasing = TRUE),]
    
    summary$Type.1[j]<-ave$Type[1]
    summary$Type.2[j]<-ave$Type[2]    
    
    if(summary$Type.1[j]=="Ts") summary$aa1[j]<-stock$MUTAA[stock$pos==position]
    if(summary$Type.2[j]=="Ts") summary$aa2[j]<-stock$MUTAA[stock$pos==position]
    if(summary$Type.1[j]=="Tv1") summary$aa1[j]<-stock$TVS1_AA[stock$pos==position]
    if(summary$Type.2[j]=="Tv1") summary$aa2[j]<-stock$TVS1_AA[stock$pos==position]
    if(summary$Type.1[j]=="Tv2") summary$aa1[j]<-stock$TVS2_AA[stock$pos==position]
    if(summary$Type.2[j]=="Tv2") summary$aa2[j]<-stock$TVS2_AA[stock$pos==position]
    
    if(summary$Type.1[j]=="Ts") summary$cpg1[j]<-stock$makesCpG[stock$pos==position]
    if(summary$Type.2[j]=="Ts")  summary$cpg2[j]<-stock$makesCpG[stock$pos==position]
    if(summary$Type.1[j]=="Tv1") summary$cpg1[j]<-stock$makesCpG.tv1[stock$pos==position]
    if(summary$Type.2[j]=="Tv1") summary$cpg2[j]<-stock$makesCpG.tv1[stock$pos==position]
    if(summary$Type.1[j]=="Tv2") summary$cpg1[j]<-stock$makesCpG.tv2[stock$pos==position]
    if(summary$Type.2[j]=="Tv2") summary$cpg2[j]<-stock$makesCpG.tv2[stock$pos==position]
    
    if(summary$Type.1[j]=="Ts")  summary$bigAA1[j]<-stock$bigAAChange[stock$pos==position]
    if(summary$Type.2[j]=="Ts")  summary$bigAA2[j]<-stock$bigAAChange[stock$pos==position]
    if(summary$Type.1[j]=="Tv1") summary$bigAA1[j]<-stock$bigAAChange.tv1[stock$pos==position]
    if(summary$Type.2[j]=="Tv1") summary$bigAA2[j]<-stock$bigAAChange.tv1[stock$pos==position]
    if(summary$Type.1[j]=="Tv2") summary$bigAA1[j]<-stock$bigAAChange.tv2[stock$pos==position]
    if(summary$Type.2[j]=="Tv2") summary$bigAA2[j]<-stock$bigAAChange.tv2[stock$pos==position]
    
    summary$occurence[j]<-cpos2$occurence[cpos2$pos==position]
    
}

write.csv(summary,"Output/MF_PID/HighFreq/Summary_highFreq_sites.csv")
