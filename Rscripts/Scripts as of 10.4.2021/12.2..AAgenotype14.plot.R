#Genotype with 13 sites
library(ggplot2)
library(gridExtra)
library(colorspace)
#library(seqinr)
#library(Rsamtools)
#library(GenomicAlignments)
#library(stringr)
library(reshape2)
library(ggpattern)
source("Rscripts/baseRscript2.R")
cols<-qualitative_hcl(6, palette="Dark3")

### 
aa14<-read.csv("Output/AA/14Genotypes.csv")
#aa13<-aa14[aa14$sites!=657]

sites<-aa14$sites
#AA sites
AAsites<-aa14$AAsites
#AA: 99 102 109 110 114 123 125 144 162 172 179 180 183 219
#NT: 297 305 326 328 340 368 373 431 485 515 536 539 548 657

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
# "E" "K" "A" "A" "T" "G" "D" "Q" "R" "T" "S" "N" "N" "Y"



apply(mutations[,c("WTAA","aa1")], 1, function(x) paste0(x["WTAA"], "\U2192",x["aa1"]))
#"D→E" "R→K" "V→A" "P→A" "A→T" "E→G" "N→D" "P→Q" "K→R" "R→T" "N→S" "T→N" "S→N" "Y→Y" 

#mutations$codon<-sapply(mutations$pos, function(x) if(x%%3==0) x=3 else x=x%%3)

HS<-read.csv("Output/HighMutfreq_sites_all.csv",stringsAsFactors = F, row.names = 1)

mutations<-merge(mutations, HS[,c(1:12)])
mutations$AApos<-AAsites
write.csv(mutations,"Output/AA/14Mutations_summary.csv")



####
tbs<-read.csv("Data/TBinfection.week.csv",stringsAsFactors = F)
samples<-read.csv("Data/SamplesNoduplicates.csv",stringsAsFactors = F)
samples$Week<-as.integer(samples$Week)

samples<-samples[samples$File.name!="Run4_18",]
#remove 34019 (no plasma sample)
samples<-samples[samples$Monkey!="A34019",]
samples<-samples[samples$Monkey!="A24018",]

list.animal<-split(samples, samples$Monkey)

#remove 34019 (no plasma sample)

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


stock<-read.csv("Output/AA/Genotype.freq/Run0_17.AA.14gtypeFreq.csv", stringsAsFactors = F, row.names = 1)
stock$Freq<-as.integer(stock$Freq)
stock$Prop<-stock$Freq/sum(stock$Freq)
stock$Week<-0
stock$Tissue<-"Stock"


files<-list.files("Output/AA/Genotype.freq/", pattern="14gtypeFreq.csv")

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

GT35<-read.csv("Output/AA/MostcommonGenotypes14_top35.csv", row.names = 1, stringsAsFactors = F)
gt35<-GT35$Genotype
gt35<-gt35[gt35!="WT"]


MUT<-data.frame(Monkey=monkeys)
Init<-data.frame(Monkey=monkeys)
for (m in 1:length(monkeyList)){
    #Select the monkey
    sample<-monkeyList[[m]]
    sample<-sample[order(sample$Week),]
    gfiles<-Gtype[as.vector(sample$File.name)]
    monkey<-names(monkeyList)[m]
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
    uniquegt<-uniquegt[!(uniquegt %in% c("WT","Singleton","Doubleton", gt35))]
    plasma$Genotype<-factor(plasma$Genotype, levels=c("WT","Singleton","Doubleton", gt35, paste(uniquegt)))
    labels=c("WT","Singleton","Doubleton", gt35[1:5])
    #ggplot(data=plasma,aes(x=Week, y=Prop, fill=Genotype))+
    #    geom_area()+
    #    geom_vline(xintercept=tbweek, col="blue", size=0.3)+ylab("Proportion")+
    #    ggtitle(monkey)+
    #    #scale_fill_manual(values = color)+
    #    scale_fill_discrete(type=color,breaks=labels)
    #ggsave(paste0("Output/AA/14GenotypeChange_",monkey,".pdf"),width = 9,height = 5)
    
    #which genotypes have the specific mutation?
    plasma$Genotype<-as.character(plasma$Genotype)
    plasma$Genotype[plasma$Genotype=="WT"]<-WTaa
    
    plasma$Genotype2<-factor(plasma$Genotype, levels=c(WTaa,"Singleton","Doubleton", gt35, paste(uniquegt)))
    labels<-c("WT","Singleton","Doubleton", gt35[1:5])
    breaks<-c(WTaa,"Singleton","Doubleton", gt35[1:5])
    
    plots<-list()
    for (k in 1:length(AAsites)){
    
        mutant<-muAA[k]
        plasma$mut<-sapply(plasma$Genotype, function(x) {
                  if (!(x=="Singleton"|x=="Doubleton")) {x=unlist(strsplit(x,"", fixed=T))
                                                        ifelse (x[k]==mutant, "Y","N")}
                    else "N"})
        #proportion of "Y" at the end and in stock
        last<-plasma[plasma$Week==max(unique(plasma$Week)),]
        y<-sum(last$Freq[last$mut=="Y"])/sum(last$Freq)*100
        stk<-plasma[plasma$Week==0,]
        b<-sum(stk$Freq[stk$mut=="Y"])/sum(stk$Freq)*100
        
        plots[[k]]<-ggplot(data=plasma,aes(x=Week, y=Prop, fill=Genotype2, pattern=mut))+
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
        
        #ggsave(paste0("Output/AA/",monkey,".pos110.T.genotypeFreq.pdf"),width = 11,height = 5)

        #ggsave(paste0("Output/AA/",monkey,".",k,".genotypeFreq.pdf"),width = 11,height = 5)
        MUT[m,(k+1)]<-y
        Init[m,(k+1)]<-b
    }
    
    #pdf(paste0("Output/AA/Figures/",monkey,".genotypeFreq.pdf"), width = 20, height = 28)
    #do.call(grid.arrange, c(plots, ncol=2))
    #dev.off()
    
}

colnames(MUT)[2:ncol(MUT)]<-paste0("AA.",AAsites)
colnames(Init)[2:ncol(Init)]<-paste0("AA.",AAsites)

MUT<-data.frame(rbind(Init[1,],MUT))
MUT$Monkey[1]<-"Stock"

write.csv(MUT, "Output/AA/Final.MutAAFreq.14sites.csv")

MUT<-read.csv("Output/AA/Final.MutAAFreq.14sites.csv", row.names = 1, stringsAsFactors = )
Mutm<-melt(MUT, id.vars="Monkey")
colnames(Mutm)[2:3]<-c("AApos","Freq")


gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}
cols = gg_color_hue(10)
cols<-c("gray60", cols[1:9])
Mutm$Monkey<-factor(Mutm$Monkey, levels=c("Stock",monkeys))

ggplot(data=Mutm, aes(x=AApos, y=Freq, fill=Monkey))+
    geom_bar(position=position_dodge(width=0.8), stat="identity", color="gray50")+
    theme_bw()+
    ylab("Frequency")+xlab('')+
    theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())+
    geom_vline(xintercept = c(1:14)+0.5, color="gray70", size=0.3)+
    scale_fill_manual(values=cols)
ggsave("Output/AA/FinalAAfreq14.summary.pdf", width = 10, height = 4)


## Based on the above, remove AAsites 110, 123, 219

Mutm<-Mutm[!(Mutm$AApos %in% c("AA.110","AA.123","AA.219")),]
mutations<-mutations[mutations$pos %in% sites,]
mutations<-mutations[!(mutations$AApos %in% c(110,123,219)),]

aalabels<-apply(mutations[,c("WTAA","aa1")], 1, function(x) paste0(x["WTAA"], "\U2192",x["aa1"]))

xlabels<-paste0("AA.",mutations$AApos,"\n",aalabels)


ggplot(data=Mutm, aes(x=AApos, y=Freq, fill=Monkey))+
    geom_bar(position=position_dodge(width=0.8), stat="identity", color="gray50")+
    theme_bw()+
    ylab("Frequency")+xlab('')+
    theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())+
    geom_vline(xintercept = c(1:14)+0.5, color="gray70", size=0.3)+
    scale_fill_manual(values=cols)+
scale_x_discrete(labels = xlabels)

ggsave("Output/AA/FinalAAfreq11.summary.png", width = 10, height = 4, units = "in", dpi=300)
ggsave("Output/AA/FinalAAfreq11.summary.pdf", width = 10, height = 4)
