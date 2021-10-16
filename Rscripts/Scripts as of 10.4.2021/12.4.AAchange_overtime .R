library(ggplot2)
library(reshape2)
library(gridExtra)
library(colorspace)
library(DataCombine)
source("Rscripts/baseRscript2.R")
cols2<-qualitative_hcl(6, palette="Dark3")

#List the aa translated files
AAs<-list.files("Output/AA/Genotype/",pattern="AA.14variants.csv$")
AAs<-AAs[AAs!="Run4_18.AA.14variants.csv"]
AA11<-list()
for (i in 1:length(AAs)){
    df<-read.csv(paste0("Output/AA/Genotype/",AAs[i]), row.names = 1,colClasses = "character" )
    remove<-which(colnames(df) %in% c("pos.110","pos.123","pos.219"))
    df<-df[,-remove]
    AA11[[i]]<-df
    names(AA11)[i]<-substr(AAs[i], 1,7)
}

aa<-read.csv("Output/Genotype11.info.csv", stringsAsFactors = F, row.names = 1)

mutations<-read.csv("Output/MF_PID/HighFreq/Summary_highFreq_sites.csv", stringsAsFactors = F, row.names = 1)
mutations<-mutations[mutations$pos %in% aa$sites,]
aalabels<-apply(mutations[,c("WTAA","aa1")], 1, function(x) paste0(x["WTAA"], "\U2192",x["aa1"]))

stock<-AA11[["Run0_17"]]

tbs<-read.csv("Data/TBinfection.week.csv",stringsAsFactors = F)
samples<-read.csv("Data/SamplesNoduplicates.csv",stringsAsFactors = F)
samples$Week<-as.integer(samples$Week)

samples<-samples[samples$File.name!="Run4_18",]
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


for (j in 1:length(AAsites)){
    position<-paste0("pos.",aa$AAsites[j])

    # AA freq at each position over time in each monkey
    Freq<-data.frame()
    for (i in 1:length(monkeys2)){
        monkey<-names(monkeys2[i])
        sample<-monkeys2[[i]]
        sample<-rbind(samples[samples$Tissue=="Stock",], sample)
        sample$Tissue[1]<-"Plasma"
        sample = sample[order(sample[,'Week']),]
        DF<-AA11[as.vector(sample$File.name)]
        
        aaFreq<-lapply(DF, function(x) data.frame(table(x[,colnames(x)==position])))
        aaFreq<-lapply(aaFreq, function(x) x<-x[x$Var1!="X",])
        aaFreq2<-lapply(aaFreq, function(x) apply(x["Freq"], 1, function(y) y/sum( x$Freq)))
        aaFreq<-mapply(cbind, aaFreq, "Prop"= aaFreq2,SIMPLIFY=F)
        
        #which amino acids are present?
        aaList<-unlist(unname(sapply(aaFreq, `[[`, "Var1")))
        aaList<-unique(as.character(aaList))
        
        Pos<-data.frame(Week=sample$Week, Tissue=sample$Tissue)
        
        for (l in 1: length(aaList)){
            L<-lapply(aaFreq, function(x) x$Prop[x$Var1==aaList[l]])
            L[sapply(L, length)==0]<-0
            Pos[,aaList[l]]<-unlist(L)
        }
        
        
        Posm<-melt(Pos, id.vars = c("Week","Tissue"))
        colnames(Posm)[3:4]<-c("AA", "Freq")
        Posm$Monkey<-monkey
        Posm$Tbweek<-tbs$tb[tbs$ids==monkey]
        Freq<-rbind(Freq, Posm)
    }
    
    ggplot(data=Freq, aes(x=Week, y=Freq, color=AA))+
        ylab("Mutation frequency")+xlab("Week")+ylim(0,1)+
        facet_wrap(~ Monkey, nrow=5, ncol=2)+
        geom_point(size=1.5)+theme_bw()+
        #scale_color_manual(values=cols2[c(1,3,5,2)], label=c("A","C",'G',"T"))+
        geom_path(data=Freq[Freq$Tissue=="Plasma",], aes(x=Week, y=Freq))+
        theme(legend.title = element_blank())+
        geom_vline(aes(xintercept=Tbweek), col="blue")+
        ggtitle(paste(position, aalabels[j]))+theme(plot.title = element_text(size=10), axis.title.y = element_text(size=8))
        ggsave(paste0("Output/AA/Overtime/", position,"_overtime.pdf"), width = 8, height=8)
}
   


SIVFiles_overview<-list.files("Output/OverviewF_PID/",pattern="filtered.overview.csv")

Overview<-list()
for (i in 1:length(SIVFiles_overview)){ 
    overviews<-read.csv(paste0("Output/OverviewF_PID/",SIVFiles_overview[i]),stringsAsFactors=F, row.names = 1)
    Overview[[i]]<-overviews
    names(Overview)[i]<-substr(paste(SIVFiles_overview[i]),start=1,stop=7)
}

df<-Overview[["Run3_14"]]

df2<-df[df$pos>=328&df$pos<331,]
aa<-read.csv("Output/AA/Run3_14AA.seq.csv", colClasses = "character")
