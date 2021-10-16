#Genotype with 10 sites
library(ggplot2)
library(gridExtra)
library(colorspace)
library(seqinr)
library(Rsamtools)
library(GenomicAlignments)
library(stringr)
library(reshape2)
source("Rscripts/baseRscript2.R")
cols<-qualitative_hcl(6, palette="Dark3")

#List the aa translated files
AAs<-list.files("Output/AA/",pattern="AA.seq.csv$")

#Interesting 14 sites 
sites<-c(297,305,326,328,340,368,373,431,485,515,536,539,548,657)

#AA sites
AAsites<-ceiling(sites/3)
site.cols<-paste0("pos.", AAsites)

aa14<-data.frame(sites=sites)
aa14$AAsites<-AAsites
#write.csv(aa14,"Output/AA/14Genotypes.csv", row.names = F)

for (i in 1:length(AAs)){
    aa<-read.csv(paste0("Output/AA/",AAs[i]), row.names = 1,colClasses = "character" )
    #select the 10 sites of interest
    col.num <- which(colnames(aa) %in% site.cols)
    aa<-aa[,col.num]
    fname<-substr(AAs[i], 1,7)
    
    write.csv(aa,paste0("Output/AA/Genotype/", fname, ".AA.14variants.csv"))
    print(fname)
}


## Calculate the genotype frequencies 
files<-list.files("Output/AA/Genotype/", pattern="\\.14variants.csv")

#Ref genotype
#ref<-read.csv("Output/Overview.ref.csv", row.names = 1,stringsAsFactors = F,colClasses = "character" )
#"DRVPAENPKRNTSY"
#Stock genotype
ref<-read.csv("Output/Overview_PIDcon/Run0_17_overview.csv",row.names = 1,stringsAsFactors = F)

#select the 14 sites of interest
ref<-ref[ref$pos %in% sites,]
WTaa<-toupper(paste0(ref$WTAA, collapse = '')) #Stock and ref genome are the same 

#exclude Run4_18 for now (too many missing sites)
files<-files[files!="Run4_18.AA.14variants.csv"]
genotypeCount<-data.frame(file=sub(".AA.14variants.csv",'',files))

#remove Run4_18: all genotypes have "N", remove from the analysis for now.
for (i in 1: length(files)){
    df<-read.csv(paste0("Output/AA/Genotype/",files[i]), stringsAsFactors = F, row.names = 1, colClasses = "character")
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
    
    write.csv(gt3,paste0("Output/AA/Genotype.freq/",fname,".AA.14gtypeFreq.csv"))
    print(fname)
}



genotypeCount$PercentWT<-genotypeCount$WTcount/genotypeCount$TotalCleanReads*100
genotypeCount$PercentSingleton<-genotypeCount$Singleton/genotypeCount$TotalCleanReads*100
genotypeCount$PercentDoubleton<-genotypeCount$Doubleton/genotypeCount$TotalCleanReads*100

write.csv(genotypeCount,"Output/AA/AA.Genotype14_count_summary.csv")

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
ggsave("Output/AA/Proprtion.Single.Double.14genotypes.pdf",height = 3,width = 8)    
    



freqfiles<-list.files("Output/AA/Genotype.freq/", pattern='14gtypeFreq.csv')

#find common genotypes across files
gtList<-list()
#commonGt<-list()

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

write.csv(GT, "Output/AA/Top50_CommonGenotypes.csv")

hist(GT$Count)
GT2<-GT[GT$Count>=10,] #72

#in more than half 
GT3<-GT[GT$Count>=20,] #35

write.csv(GT3,"Output/AA/MostcommonGenotypes14_top35.csv")

