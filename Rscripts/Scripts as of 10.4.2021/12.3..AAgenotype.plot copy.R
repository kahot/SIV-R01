#Genotype with 10 sites
library(ggplot2)
library(gridExtra)
library(colorspace)
library(seqinr)
library(Rsamtools)
library(GenomicAlignments)
library(stringr)
library(reshape2)
library(ggpattern)
source("Rscripts/baseRscript2.R")
cols<-qualitative_hcl(6, palette="Dark3")

### 

sites<-c(305,328,340,368,373,431,475,515,536,539)
#AA sites
AAsites<-ceiling(sites/3)

#Ref genotype
ref<-read.csv("Output/Overview.ref.csv", row.names = 1,stringsAsFactors = F)
ref<-ref[ref$pos %in% sites,]
WTaa<-toupper(paste0(ref$RefAA, collapse = ''))


mutations<-read.csv("Output/MF_PID/HighFreq/Summary_highFreq_sites.csv", stringsAsFactors = F, row.names = 1)
mutations<-mutations[order(mutations$occurence, decreasing = T),]
mutations<-mutations[mutations$pos %in% sites,]
#The mutated AA
muAA<-mutations$aa1  
#[1] "K" "A" "T" "G" "D" "Q" "N" "G" "S" "T"

#mutations$codon<-sapply(mutations$pos, function(x) if(x%%3==0) x=3 else x=x%%3)

HS<-read.csv("Output/HighMutfreq_sites_all.csv",stringsAsFactors = F, row.names = 1)

mutations<-merge(mutations, HS[,c(1:12)])
mutations$AApos<-AAsites
write.csv(mutations,"Output/AA/10Mutations_summary.csv")

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


stock<-read.csv("Output/AA/Genotype.freq/Run0_17.AA.10gtypeFreq.csv", stringsAsFactors = F, row.names = 1)
stock$Freq<-as.integer(stock$Freq)
stock$Prop<-stock$Freq/sum(stock$Freq)
stock$Week<-0
stock$Tissue<-"Stock"


files<-list.files("Output/AA/Genotype.freq/", pattern="gtypeFreq.csv")

Gtype<-list()
for (i in 1: length(files)){
    df<-read.csv(paste0("Output/AA/Genotype.freq/",files[i]), stringsAsFactors = F, row.names = 1)
    fname<-substring(files[i],1,7)
    Gtype[[i]]<-df
    names(Gtype)[i]<-fname
}

labels<-c("WT","Singleton","Doubleton")
legends <-c(WTaa,"Singleton","Doubleton")

m=1
#for (m in 1:length(monkeys2)){
    #Select the monkey
    sample<-monkeys2[[m]]
    sample<-sample[order(sample$Week),]
    gfiles<-Gtype[as.vector(sample$File.name)]
    monkey<-names(monkeys2)[m]
    tbweek<-tbs$tb[tbs$ids==monkey]
    
    
    #Createa a vector of all unique genotypes for each monkey
    gtypes<-unique(stock$Var1)
    DF<-list()
    DF[[1]]<-stock
    for (j in 1:length(gfiles)){
        df2<-gfiles[[j]]
        df2$Freq<-as.integer(df2$Freq)
        df2$Prop<-df2$Freq/sum(df2$Freq)
        df2$Week<-sample$Week[j]
        df2$Tissue<-paste0(sample$Tissue2[j])
        
        gtypes<-c(gtypes,unique(df2$Var1))
        DF[[(j+1)]]<-df2
    }
    gtypes<-unique(gtypes) 
    
    GT2<-data.frame(Genotype=gtypes)
    Gdata<-data.frame()
    for (j in 1:length(DF)){
        df<-DF[[j]]
        colnames(df)[1]<-"Genotype"
        merged<-merge(GT2,df[,c(1,2,4)], by="Genotype", all=T)
        merged[is.na(merged)]<-0
        merged$Week<-df$Week[1]
        merged$Tissue<-df$Tissue[1]
        
        Gdata<-rbind(Gdata,merged)
    }
    
    #write.csv(Gdata,paste0("Output/Genotype/Freq/",monkey,".GenotypeFreqNoDoubleton.csv"))
    
    #select plasma only
    plasma<-Gdata[Gdata$Tissue=="Plasma"|Gdata$Tissue=="Stock",]
    
    #ggplot(data=plasma,aes(x=Week, y=Prop, fill=Genotype))+
    #    geom_area()+theme(legend.position = "none")+
    #    geom_vline(xintercept=tbweek, col="blue", size=0.3)+ylab("Proportion")+
    #    ggtitle(monkey)
    #
    #ggsave(paste0("Output/AA/",monkey,".Genotype_AreaStacked.pdf"),width = 9,height = 5)
    
    #which genotypes have the specific mutation?
    plasma$Genotype[plasma$Genotype=="WT"]<-WTaa
    
    for (k in 1:length(AAsites)){
        mutant<-muAA[k]
        plasma$mut<-sapply(plasma$Genotype, function(x) {
                  if (!(x=="Singleton"|x=="Doubleton")) {x=unlist(strsplit(x,"", fixed=T))
                                                        ifelse (x[k]==mutant, "Y","N")}
                    else "N"})
        
        ggplot(data=plasma,aes(x=Week, y=Prop, fill=Genotype, pattern=mut))+
            geom_area_pattern(color = "gray80",size=0.05, 
                              pattern_fill = "gray60",
                              pattern_angle = 45,
                              pattern_density = 0.1,
                              pattern_spacing = 0.025,
                              pattern_key_scale_factor = 0.6)+
            theme_bw()+
            theme(legend.title = element_blank())+ylab("Frequency")+
            scale_pattern_manual(values = c(Y = "stripe", N = "none") ) +
            geom_vline(xintercept=tbweek, col="blue", size=0.3)+ylab("Proportion")+
            ggtitle(paste(monkey,"AApos",AAsites[k],"mutant",mutant ) )+ 
            scale_fill_discrete(breaks=legends, labels=labels)+
            guides(pattern = guide_legend(override.aes = list(fill = "white")),
                   fill = guide_legend(override.aes = list(pattern = "none")))
        
        ggsave(paste0("Output/AA/",monkey,".",k,".AA.AreaStacked.pdf"),width = 11,height = 5)
    }
    
    mutFreq<-data.frame(Week=c(0,sample$Week[sample$Tissue2=="Plasma"]))
    #use aggregate
    for (k in 1:length(AAsites)){
        mutant<-muAA[k]
        
        plasma$mut<-sapply(plasma$Genotype, function(x) {
            if (!(x=="Singleton"|x=="Doubleton")) {x=unlist(strsplit(x,"", fixed=T))
            ifelse (x[k]==mutant, "Y","N")}
            else "N"})
        mutFreq[m,mutant]<-sum(plasma$Prop[plasma$mut=="Y"])
        #Freq of mutants at the end
        
        
    }
    
    #stacked bar chart of all files      
    Gdata$label<-paste0("Wk ",Gdata$Week," ",Gdata$Tissue)
    L<-unique(Gdata$label)
    Gdata$label<-factor(Gdata$label,levels=L)
    ggplot()+
        geom_bar(data=Gdata,aes(x=label, y=Prop, fill=Genotype), stat="identity")+
        ylab("Proportion")+
        theme(legend.position = "none",axis.title.x = element_blank(),axis.text.x = element_text(angle=45, hjust=1))
    ggsave(paste0("Output/Genotype/Freq/",monkey,".Genotype_freq.stackbar.pdf"),width = 12,height = 5)
    
    #interactive streamhrgraph
    plasma$Date<-as.Date(plasma$Week*7, origin="2020-01-01")
    pp<-streamgraph(plasma, key="Genotype", value="Prop", date="Date", height="300px", width="800px")%>%
        sg_axis_x(tick_interval = 1,tick_unit="month",tick_format ="%m%d")%>%
        sg_legend(show=TRUE, label="Genotype: ")
    
    saveWidget(pp,file=paste0("Output/Genotype/Streamplot.",monkey,".html"))
    
}

##adding a legend for certain types for 1 monkey 

#1. A21918
legends <- c("WT","Singleton","Doubleton")
ggplot(data=plasma,aes(x=Week, y=Prop, fill=Genotype))+
    geom_area()+
    geom_vline(xintercept=tbweek, col="blue")+ylab("Proportion")+
    ggtitle(monkey)+
    scale_fill_discrete(breaks=legends)
ggsave(paste0("Output/Genotype/Freq/",monkey,"areaStacked_withLegend.pdf"),width=10,height=5)



