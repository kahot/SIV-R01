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

#sites<-c(305,328,340,368,373,431,475,515,536,539)

sites<-c(297,305,326,340,368,373,431,485,515,536,539)

#AA sites
AAsites<-ceiling(sites/3)
#AA: 102 110 114 123 125 144 159 172 179 180
#NT: 305 328 340 368 373 431 475 515 536 539
#    "K" "A" "T" "G" "D" "Q" "G" "T" "S" "N"

#Ref genotype
ref<-read.csv("Output/Overview.ref.csv", row.names = 1,stringsAsFactors = F)
ref<-ref[ref$pos %in% sites,]
WTaa<-toupper(paste0(ref$RefAA, collapse = ''))


mutations<-read.csv("Output/MF_PID/HighFreq/Summary_highFreq_sites.csv", stringsAsFactors = F, row.names = 1)
mutations<-mutations[order(mutations$occurence, decreasing = T),]
mutations<-mutations[mutations$pos %in% sites,]
mutations<-mutations[order(mutations$pos),]

#The mutated AA
muAA<-mutations$aa1  
#[1] "K" "A" "T" "G" "D" "Q" "G" "T" "S" "N"


apply(mutations[,c("WTAA","aa1")], 1, function(x) paste0(x["WTAA"], "\U2192",x["aa1"]))
paste0(mutations$WTAA,"->",mutations$aa1)


#mutations$codon<-sapply(mutations$pos, function(x) if(x%%3==0) x=3 else x=x%%3)

HS<-read.csv("Output/HighMutfreq_sites_all.csv",stringsAsFactors = F, row.names = 1)

mutations<-merge(mutations, HS[,c(1:12)])
mutations$AApos<-AAsites
#write.csv(mutations,"Output/AA/10Mutations_summary.csv")



####
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


stock<-read.csv("Output/AA/Genotype.freq/Run0_17.AA.11gtypeFreq.csv", stringsAsFactors = F, row.names = 1)
stock$Freq<-as.integer(stock$Freq)
stock$Prop<-stock$Freq/sum(stock$Freq)
stock$Week<-0
stock$Tissue<-"Stock"


files<-list.files("Output/AA/Genotype.freq/", pattern="11gtypeFreq.csv")

Gtype<-list()
for (i in 1: length(files)){
    df<-read.csv(paste0("Output/AA/Genotype.freq/",files[i]), stringsAsFactors = F, row.names = 1)
    fname<-substring(files[i],1,7)
    Gtype[[i]]<-df
    names(Gtype)[i]<-fname
}

labels<-c("WT","Singleton","Doubleton")
legends <-c(WTaa,"Singleton","Doubleton")


color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
color<-color[color!="white"]
color<-color[!color %in% c("turquoise", "skyblue", "ivory")]
color<-c("turquoise", "skyblue", "ivory", color)

GT40<-read.csv("Output/AA/MostcommonGenotypes_40.csv", row.names = 1, stringsAsFactors = F)
gt40<-GT40$Genotype
gt40<-gt40[gt40!="WT"]



MUT<-data.frame(Monkey=monkeys)
for (m in 2:length(monkeys2)){
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
    
    
    #select plasma only
    plasma<-Gdata[Gdata$Tissue=="Plasma"|Gdata$Tissue=="Stock",]
    uniquegt<-unique(plasma$Genotype)
    uniquegt<-uniquegt[!(uniquegt %in% c("WT","Singleton","Doubleton", gt40))]
    plasma$Genotype<-factor(plasma$Genotype, levels=c("WT","Singleton","Doubleton", gt40, paste(uniquegt)))
    labels=c("WT","Singleton","Doubleton", gt40[1:5])
    ggplot(data=plasma,aes(x=Week, y=Prop, fill=Genotype))+
        geom_area()+
        geom_vline(xintercept=tbweek, col="blue", size=0.3)+ylab("Proportion")+
        ggtitle(monkey)+
        #scale_fill_manual(values = color)+
        scale_fill_discrete(type=color,breaks=labels)
    ggsave(paste0("Output/AA/11GenotypeChange_",monkey,".pdf"),width = 9,height = 5)
    
    #which genotypes have the specific mutation?
    plasma$Genotype<-as.character(plasma$Genotype)
    plasma$Genotype[plasma$Genotype=="WT"]<-WTaa
    
    plasma$Genotype2<-factor(plasma$Genotype, levels=c(WTaa,"Singleton","Doubleton", gt40, paste(uniquegt)))
    labels<-c("WT","Singleton","Doubleton", gt40[1:5])
    breaks<-c(WTaa,"Singleton","Doubleton", gt40[1:5])
    
    plots<-list()
    for (k in 1:length(AAsites)){
        mutant<-muAA[k]
        plasma$mut<-sapply(plasma$Genotype, function(x) {
                  if (!(x=="Singleton"|x=="Doubleton")) {x=unlist(strsplit(x,"", fixed=T))
                                                        ifelse (x[k]==mutant, "Y","N")}
                    else "N"})
        #proportion of "Y" at the end
        last<-plasma[plasma$Week==max(unique(plasma$Week)),]
        y<-sum(last$Freq[last$mut=="Y"], na.rm=T)/sum(last$Freq, na.rm=T)*100
        
       
        plots[[k]]<-
            ggplot(data=plasma,aes(x=Week, y=Prop, fill=Genotype2, pattern=mut))+
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
            ggtitle(paste0(monkey," AApos:",AAsites[k]," mutant:",mutant, " ", round(y, digits=1), "%") )+ 
            scale_fill_discrete(type=color,breaks=breaks, labels= labels)+
            guides(pattern = guide_legend(override.aes = list(fill = "white")),
                   fill = guide_legend(override.aes = list(pattern = "none")))
        
            #proportion of "Y" at the end
            #last<-plasma[plasma$Week==max(unique(plasma$Week)),]
            #y<-sum(last$Freq[last$mut=="Y"])/sum(last$Freq)*100
            MUT[m,(k+1)]<-y
    }
    
    pdf(paste0("Output/AA/Figures/",monkey,".11genotypeFreq.pdf"), width = 20, height = 20)
    do.call(grid.arrange, c(plots, ncol=2))
    dev.off()
    
}

colnames(MUT)[2:11]<-paste0("AA.",AAsites)
write.csv(MUT, "Output/AA/Final.MutAAFreq.10sites.csv")

#remove the moneky 34219
MUT<-MUT[!is.na(MUT$AA.102),]
Mutm<-melt(MUT, id.vars="Monkey")
colnames(Mutm)[2:3]<-c("AApos","Freq")
ggplot(data=Mutm, aes(x=AApos, y=Freq, fill=Monkey))+
    geom_bar(position=position_dodge(width=0.8), stat="identity", color="gray50")+
    theme_bw()+
    ylab("Frequency")+xlab('')+
    theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())+
    geom_vline(xintercept = c(1:10)+0.5, color="gray70", size=0.3)
ggsave("Output/AA/FinalAAfreq.summary.pdf", width = 8, height = 4)



##adding a legend for certain types for 1 monkey 

#1. A21918
legends <- c("WT","Singleton","Doubleton")
ggplot(data=plasma,aes(x=Week, y=Prop, fill=Genotype))+
    geom_area()+
    geom_vline(xintercept=tbweek, col="blue")+ylab("Proportion")+
    ggtitle(monkey)+
    scale_fill_discrete(breaks=legends)
ggsave(paste0("Output/Genotype/Freq/",monkey,"areaStacked_withLegend.pdf"),width=10,height=5)



