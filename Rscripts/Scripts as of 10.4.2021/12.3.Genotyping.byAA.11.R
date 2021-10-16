# Based on 12.2, use 11 sites 
library(ggplot2)
library(gridExtra)
library(colorspace)
library(reshape2)
source("Rscripts/baseRscript2.R")


aa<-read.csv("Output/AA/14Genotypes.csv")
aa<-aa[!(aa$AAsites==110|aa$AAsites==123|aa$AAsites==219),]

sites<-aa$sites
AAsites<-aa$AAsites

## Calculate the genotype frequencies 
files<-list.files("Output/AA/Genotype/", pattern="\\.14variants.csv")

ref<-read.csv("Output/Overview_PIDcon/Run0_17_overview.csv",row.names = 1,stringsAsFactors = F)

#select the 14 sites of interest
ref<-ref[ref$pos %in% sites,]
WTaa<-toupper(paste0(ref$WTAA, collapse = '')) #Stock and ref genome are the same 
#DRVANPKRNTS

#exclude Run4_18 for now (too many missing sites)
files<-files[files!="Run4_18.AA.14variants.csv"]
genotypeCount<-data.frame(file=sub(".AA.14variants.csv",'',files))

for (i in 1: length(files)){
    df<-read.csv(paste0("Output/AA/Genotype/",files[i]), stringsAsFactors = F, row.names = 1, colClasses = "character")
    remove<-which(colnames(df) %in% c("pos.110","pos.123","pos.219"))
    df<-df[,-remove]
    fname<-substring(files[i],1,7)
    
    #calculate the freq of each genotype
    df$Genotype<-apply(df[1:ncol(df)], 1, function(x) paste0(x,collapse = ''))
    genotypeCount$Total[i]<-length(unique(df$Genotype))
    
    #remove the sequences with unkown aa  ('X')  
    df<-df[!grepl("X", df$Genotype, fixed=T),]
    genotypeCount$no_Xs[i]<-length(unique(df$Genotype))
    
    #Genotype freq.
    gt<-data.frame(table(df$Genotype))
    gt<-gt[order(gt$Freq, decreasing = T),]
    gt$Var1<-as.character(gt$Var1)
    
    w<-which(gt$Var1==WTaa)
    #Replace the wildtype (Reference) with WT
    gt$Var1[w]<-"WT"
    genotypeCount$TotalCleanReads[i]<-nrow(df)
    genotypeCount$WTcount[i]<-ifelse(length(gt$Freq[gt$Var1=="WT"])==0, 0,gt$Freq[gt$Var1=="WT"])
    
    #remove singletons
    gt2<-gt[gt$Freq!=1,] 
    s<-nrow(gt[gt$Freq==1,])
    gt2[nrow(gt2)+1,]<-c("Singleton",s)
    gt2$Freq<-as.integer(gt2$Freq)
    #gt2$ID<-fname
    genotypeCount$Singleton[i]<-s
    
    #remove doubletons
    gt3<-gt2[gt2$Freq!=2,] 
    d<-nrow(gt2[gt2$Freq==2,])
    gt3[nrow(gt3)+1,]<-c("Doubleton",d)
    gt3$Freq<-as.integer(gt3$Freq)
    gt3$ID<-fname
    genotypeCount$Doubleton[i]<-d
    
    write.csv(gt3,paste0("Output/AA/Genotype.freq/",fname,".AA.11gtypeFreq.csv"))
    print(fname)
}



genotypeCount$PercentWT<-genotypeCount$WTcount/genotypeCount$TotalCleanReads*100
genotypeCount$PercentSingleton<-genotypeCount$Singleton/genotypeCount$TotalCleanReads*100
genotypeCount$PercentDoubleton<-genotypeCount$Doubleton/genotypeCount$TotalCleanReads*100

write.csv(genotypeCount,"Output/AA/AA.Genotype11_count_summary.csv")

#look at the proportion of WT,singletons and doubletons
gt<-genotypeCount[,c(1,8:10)]
gtm<-melt(gt)
gtm$variable<-factor(gtm$variable, levels=c("PercentDoubleton","PercentSingleton","PercentWT"))
gtm<-gtm[order(gtm$variable),]
ggplot(gtm,aes(x=file,y=value, fill=variable))+
    geom_bar(position="stack",stat='identity')+
    scale_fill_manual(values=cols[c(5,3,1)], labels=c("% Doubleton","% Singleton","% WT"))+
    theme_bw()+ylab('')+xlab('')+
    theme(axis.text.x = element_text(angle=90, size=5), legend.title = element_blank())
ggsave("Output/AA/Proprtion.Single.Double.11genotypes.pdf",height = 3,width = 8)    
    



freqfiles<-list.files("Output/AA/Genotype.freq/", pattern='11gtypeFreq.csv')

#find common genotypes across files
gtList<-list()

GT<-data.frame(Genotype="WT")
for (i in 1:length(freqfiles)){

    fname<-substring(freqfiles[i],1,7)
    #for (i in 1:length(gfiles)){
    df<-read.csv(paste0("Output/AA/Genotype.freq/",freqfiles[i]), row.name=1,stringsAsFactors = F)
    colnames(df)[1]<-"Genotype"
    #select top 50 genotypes
    df<-df[!(df$Genotype=="Singleton"|df$Genotype=="Doubleton"),]
    df<-df[order(df$Freq, decreasing = T),]
    
    if (nrow(df)>50) df<-df[1:50,]
    
    GT<-merge(GT,df[,1:2], by="Genotype", all=T)
    colnames(GT)[i+1]<-fname
    
    gtList[[i]]<-GT
    names(gtList)[i]<-fname
    print(fname)
}

#Count the occurrence of each genotype
GT$Count<-rowSums(!is.na(GT[-which(names(GT)=="Genotype")]))

GT<-GT[order(GT$Count,decreasing = T),]

write.csv(GT, "Output/AA/Top50_CommonGenotypes11.csv")

hist(GT$Count)
GT2<-GT[GT$Count>=10,] #77

#in more than 20 files
GT3<-GT[GT$Count>=20,] #34

write.csv(GT3,"Output/AA/MostcommonGenotypes11_top34.csv")

ggplot(data=GT, aes(x=Count))+
    geom_histogram(fill="lightblue", bins=30, color="gray60")+
    xlab("Genotype appearance out of 92 samples")+ylab('')
ggsave("Output/AA/Genotype11.freq.pdf", width = 6, height = 3)    

ggplot(data=GT2, aes(x=Count))+
    geom_histogram(fill="lightblue", bins=30, color="gray60")+
    xlab("Genotype appearance out of 92 samples")+ylab('')
ggsave("Output/AA/Genotype11.freq.zoom.pdf", width = 5, height = 3)    



#### Plots
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

#Ref genotype
ref<-read.csv("Output/Overview.ref.csv", row.names = 1,stringsAsFactors = F)
ref<-ref[ref$pos %in% aa$sites,]
WTaa<-toupper(paste0(ref$RefAA, collapse = ''))


labels<-c("WT","Singleton","Doubleton")
legends <-c(WTaa,"Singleton","Doubleton")


color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
color<-color[color!="white"]
color<-color[!color %in% c("turquoise", "skyblue", "ivory")]
color<-c("turquoise", "skyblue", "ivory", color)

GT35<-read.csv("Output/AA/MostcommonGenotypes11_top34.csv", row.names = 1, stringsAsFactors = F)
gt35<-GT35$Genotype
gt35<-gt35[gt35!="WT"]

mutations<-read.csv("Output/MF_PID/HighFreq/Summary_highFreq_sites.csv", stringsAsFactors = F, row.names = 1)
mutations<-mutations[mutations$pos %in% aa$sites,]

aa<-merge(aa, mutations[,c("pos","aa1")], by.x = "sites",by.y = "pos")

#write.csv(aa,"Output/Genotype11.info.csv")

for (m in 1:length(monkeys2)){
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
    uniquegt<-uniquegt[!(uniquegt %in% c("WT","Singleton","Doubleton", gt35))]
    plasma$Genotype<-factor(plasma$Genotype, levels=c("WT","Singleton","Doubleton", gt35, paste(uniquegt)))
    labels=c("WT","Singleton","Doubleton", gt35[1:5])
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
    
    plasma$Genotype2<-factor(plasma$Genotype, levels=c(WTaa,"Singleton","Doubleton", gt35, paste(uniquegt)))
    labels<-c("WT","Singleton","Doubleton", gt35[1:5])
    breaks<-c(WTaa,"Singleton","Doubleton", gt35[1:5])
    
    plots<-list()
    for (k in 1:nrow(aa)){
        
        mutant<-aa$aa1[k]
        plasma$mut<-sapply(plasma$Genotype, function(x) {
            if (!(x=="Singleton"|x=="Doubleton")) {x=unlist(strsplit(x,"", fixed=T))
            ifelse (x[k]==mutant, "Y","N")}
            else "N"})
        #proportion of "Y" at the end and in stock
        last<-plasma[plasma$Week==max(unique(plasma$Week)),]
        y<-sum(last$Freq[last$mut=="Y"])/sum(last$Freq)*100
        
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

        #ggsave(paste0("Output/AA/",monkey,".",k,".genotypeFreq.pdf"),width = 11,height = 5)
       
    }
    
    pdf(paste0("Output/AA/Figures/",monkey,".genotype11Freq.pdf"), width = 20, height = 24)
    do.call(grid.arrange, c(plots, ncol=2))
    dev.off()
    
}


table(mutations$aa1)
table(mutations$WTAA)
table(mutations$bigAA1)
table(mutations$cpg1)
table(mutations$ref)

mutations<-read.csv("Output/MF_PID/HighFreq/Summary_highFreq_sites.csv", stringsAsFactors = F, row.names = 1)
table(mutations$aa1)
table(mutations$WTAA)
table(mutations$bigAA1)
table(mutations$cpg1)
table(mutations$ref)
