library(reshape2)
library(gridExtra)
library(colorspace)
library(DataCombine)
source("Rscripts/baseRscript.R")
cols2<-qualitative_hcl(6, palette="Dark3")

SIVfiles<-list.files("Output/OverviewF_PID/",pattern=".csv")

Overview<-list()
for (i in 1:length(SIVfiles)){ 
    overviews<-read.csv(paste0("Output/OverviewF_PID/",SIVfiles[i]),stringsAsFactors=F, row.names = 1)
    Overview[[i]]<-overviews
    names(Overview)[i]<-substr(paste(SIVfiles[i]),start=1,stop=7)
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
Ref<-read.csv("Output/Overview.ref.csv", stringsAsFactors = F, row.names = 1)


Pos<-list()
Sites<-list()
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
    df1$Mutation<-Ref$Type.r[Ref$pos %in% sites]
    
    sites1<-lapply(ovDF, function(x) x$pos[x$freq.transv1.ref>0.1])
    sites1<-unique(unlist(sites1))
    sites1<-sites1[!is.na(sites1)]
    sites1<- sites1[order(sites1)]
    if (length(sites1)>0) {
        df2<-data.frame(pos=sites1, type="Tv1")
        df2$Mutation<-Ref$Type.tv1.r[Ref$pos %in% sites1]
        }
    if (length(sites1)==0) df2<-data.frame(pos=NA, type="Tv1", Mutation=NA)
    
    sites2<-lapply(ovDF, function(x) x$pos[x$freq.transv2.ref>0.1])
    sites2<-unique(unlist(sites2))
    sites2<-sites2[!is.na(sites2)]
    sites2<-sites2[order(sites2)]
    if (length(sites2)>0) {
        df3<-data.frame(pos=sites2, type="Tv2")
        df3$Mutation<-Ref$Type.tv2.r[Ref$pos %in% sites2]
        }
    if (length(sites2)==0) df3<-data.frame(pos=NA, type="Tv2", Mutation=NA) 
    
    mutSites<-rbind(df1,rbind(df2,df3))
    
    mutSites<-mutSites[!is.na(mutSites$pos),]
    msites<-unique(c(sites,sites1,sites2))
    msites<-msites[order(msites)]
    mutSites$Monkey<-monkey
    Sites[[i]]<-mutSites
    Pos[[i]]<-msites
    
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
    if (nrow(mfDF1[[1]])>1) Tv1<-data.frame(t(Tv1)) 
    
    nt2<-mfDF1[[1]][,"ref"]
    colnames(Tv1)<-paste0("Pos_",sites1, ' ',nt2, ' Tv1')
    Tv1<-InsertRow(Tv1, c(stock$freq.transv1.ref[stock$pos %in% sites1]), RowNum=1)
    
    #Tranv2
    Tv2<-unname(sapply(mfDF2, `[[`, "freq.transv2.ref"))
    if (nrow(mfDF2[[1]])==1) Tv2<-data.frame(Tv2)
    if (nrow(mfDF2[[1]])>1) Tv2<-data.frame(t(Tv2))
    
    nt3<-mfDF2[[1]][,"ref"]
    colnames(Tv2)<-paste0("Pos_",sites2, ' ',nt3, ' Tv2')
    Tv2$Week<-sample$Week
    Tv2<-InsertRow(Tv2, c(stock$freq.transv2.ref[stock$pos %in% sites2],0), RowNum=1)
    
    M<-cbind(Ts, Tv1,Tv2)
    
    colorder<-order(colnames(M))
    M<-M[,colorder]
    
    Pos[[i]]<-msites
    names(Pos)[i]<-monkey
               
    Mm<-melt(M, id.vars = "Week")
    colnames(Mm)[2:3]<-c("Position", "MF")
    rown<-ceiling((ncol(M)-1)/5)
    
    tbweek<-tbs$tb[tbs$ids==monkey]
    ggplot(data=Mm, aes(x=Week, y=MF))+
        ylab("Mutation frequency")+xlab("Weeks")+
        facet_wrap(~ Position, nrow=rown, ncol=5)+
        geom_point(size=2, color=cols2[1])+theme_bw()+ggtitle(paste(monkey))+
        geom_vline(xintercept=tbweek, col="blue")
    ggsave(paste0("Output/Timeseries_PID/", monkey, ".allsites.persite.pdf"), heigh=rown*1.5, width =11.5)
}


#Find the overlapping positions
hPos<-do.call(rbind, Sites)
#Remove the synonymous mutations
hPos<-hPos[hPos$Mutation!="syn",]

#create an ID for mutations
hPos$ID<-paste0(hPos$pos,"_",hPos$type)
positions<-unique(hPos$ID) #53

cpos<-data.frame(Mutation=positions)
cpos$Pos<-as.integer(gsub("_.*",'',cpos$Mutation))
cpos$Type<-gsub("\\d\\d\\d_",'',cpos$Mutation)

cpos<-cpos[order(cpos$Pos),]

counts<-data.frame(table(hPos$ID))
colnames(counts)<-c("Mutation","Freq")

cpos<-merge(cpos,counts, by="Mutation")

write.csv(cpos, "Output/MF_PID/HighFreq/HighMF_sites.csv") #53 total


#Find the common sites that appear in >5 samples
high<-cpos[cpos$Freq>5,] #14 mutations
highpos<-unique(high$Pos) #13 sites


df<-stock

for (j in 1:length(highpos)){
    position<-highpos[j]
    #nucleotide at Ref, Ts, tv1, and tv2
    re<-Ref$ref[Ref$pos==position]
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
    
    type<-high$Type.r[high$pos==position]
    type1<-high$Type.tv1.r[high$pos==position]
    type2<-high$Type.tv2.r[high$pos==position]
    
    ggplot(data=Freq, aes(x=Week, y=MF, color=nuc))+
        ylab("Mutation frequency")+xlab("Week")+ylim(0,1)+
        facet_wrap(~ Monkey, nrow=4, ncol=3)+
        geom_point(size=1.5)+theme_bw()+
        scale_color_manual(values=cols2[c(1,3,5,2)], label=c("A","C",'G',"T"))+
        geom_path(data=Freq[Freq$Tissue=="Plasma",], aes(x=Week, y=MF))+
        theme(legend.title = element_blank())+
        geom_vline(aes(xintercept=Tbweek), col="blue")+
        ggtitle(paste0("Pos.",position) )+theme(plot.title = element_text(size=10), axis.title.y = element_text(size=8))
    ggsave(paste0("Output/MF_PID/HighFreq/Pos", position, "_overtime.pdf"), width = 12, height=8)
}
   


table(cpos$Freq, cpos$Type)
# 1  2  3  5  6  8  10 11 
#23  8  4  3  4  1   4  6
#6 mutations occurred in all monkeys. 1 occurred in 10. 

#highpos
#428 is not interesting. Remove it

highpos<-highpos[highpos!=428]

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
    
    if(summary$Type.1[j]=="Ts")  summary$mutNT[j]<-transition(summary$ref[j])
    if(summary$Type.1[j]=="Tv1") summary$mutNT[j]<-transv1(summary$ref[j])
    if(summary$Type.1[j]=="Tv2") summary$mutNT[j]<-transv2(summary$ref[j])
    
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
    
    if (position==326|position==328|position==431|position==539) { 
                fq<-cpos[cpos$Pos==position,]
                summary$occurence[j]<-cpos$Freq[cpos$Pos==position & cpos$Type==summary$Type.1[j]]}
    else {summary$occurence[j]<-cpos$Freq[cpos$Pos==position]}
    
}


summary$AApos<-ceiling(summary$pos/3)
write.csv(summary,"Output/MF_PID/HighFreq/Summary_highFreq_sites.csv")

