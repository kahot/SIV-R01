#Genotype with 11 sites
library(ggplot2)
library(gridExtra)
library(colorspace)
library(seqinr)
library(stringr)
library(reshape2)
cols<-qualitative_hcl(6, palette="Dark3")

#List the aa translated files
AAs<-list.files("Output/AA/",pattern="AA.seq.csv$")
AAs<-list.files("Output/AA/",pattern="^Run8.+AA.seq.csv$")


#Interesting 11 sites 
hsites<-read.csv("Output/MF_PID/HighFreq/Summary_highFreq_sites.csv", stringsAsFactors = F)
hsites<-hsites[hsites$pos!=485,]

sites<-hsites$pos

#AA sites
AAsites<-ceiling(sites/3)
#NT:297 305 326 328 340 373 431 536 539 548 
#AA: 99 102 109 110 114 125 144 179 180 183 

site.cols<-paste0("pos.", AAsites)

for (i in 1:length(AAs)){
    aa<-read.csv(paste0("Output/AA/",AAs[i]), row.names = 1,colClasses = "character" )
    #select the 11 sites of interest
    col.num <- which(colnames(aa) %in% site.cols)
    aa<-aa[,col.num]
    fname<-substr(AAs[i], 1,7)
    
    write.csv(aa,paste0("Output/AA/Genotype10/", fname, ".AA.10variants.csv"))
    print(fname)
}


## Calculate the AA genotype frequencies 
files<-list.files("Output/AA/Genotype10/", pattern=".10variants.csv")

#Ref genotype
ref<-read.csv("Output/Overview.ref.csv", row.names = 1,stringsAsFactors = F,colClasses = "character" )
#select the 10 sites of interest)
ref<-ref[ref$pos %in% sites,]
WTaa<-toupper(paste0(ref$RefAA, collapse = ''))

genotypeCount<-data.frame(file=sub(".AA.10variants.csv",'',files))

for (i in 1:length(files)){
    df<-read.csv(paste0("Output/AA/Genotype10/",files[i]), stringsAsFactors = F, row.names = 1, colClasses = "character")
    fname<-substring(files[i],1,7)
    
    #calculate the freq of each genotype
    df$Genotype<-apply(df[1:ncol(df)], 1, function(x) paste0(x,collapse = ''))
    genotypeCount$Total[i]<-length(unique(df$Genotype))
    
    #remove the sequences with unknown aa  ('X')  
    df2<-df[!grepl("X", df$Genotype, fixed=T),]
    genotypeCount$no_Xs[i]<-length(unique(df2$Genotype))
    
    #Genotype freq. including X (otherwise many have no genotypes).
    gt<-data.frame(table(df$Genotype))
    gt<-gt[order(gt$Freq, decreasing = T),]
    gt$Var1<-as.character(gt$Var1)
    
    w<-which(gt$Var1==WTaa)
    #Replace the wildtype (Reference) with WT
    gt$Var1[w]<-"WT"
    genotypeCount$TotalReads[i]<-nrow(df)
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
    
    write.csv(gt3,paste0("Output/AA/Genotype.freq10/",fname,".AA.10gtypeFreq.csv"))
    print(fname)
}

#genotypeCount$PercentWT<-genotypeCount$WTcount/genotypeCount$TotalCleanReads*100
genotypeCount$PercentSingleton<-genotypeCount$Singleton/genotypeCount$TotalReads*100
genotypeCount$PercentDoubleton<-genotypeCount$Doubleton/genotypeCount$TotalReads*100

write.csv(genotypeCount,"Output/AA/AA.Genotype10_count_summary.csv")


#### Find common genotypes across monkeys/samples
freqfiles<-list.files("Output/AA/Genotype.freq10/", pattern='gtypeFreq.csv')
gtList<-list()
GT<-data.frame(Genotype="WT")
for (i in 1:length(freqfiles)){
    
    fname<-substring(freqfiles[i],1,7)
    
    df<-read.csv(paste0("Output/AA/Genotype.freq10/",freqfiles[i]), row.name=1,stringsAsFactors = F)
    colnames(df)[1]<-"Genotype"
    #select top 70 genotypes
    df<-df[!(df$Genotype=="Singleton"|df$Genotype=="Doubleton"),]
    df<-df[order(df$Freq, decreasing = T),]
    
    if (nrow(df)>70) df<-df[1:70,]
    
    GT<-merge(GT,df[,1:2], by="Genotype", all=T)
    colnames(GT)[i+1]<-fname
    
    gtList[[i]]<-GT
    names(gtList)[i]<-fname
    print(fname)
}

#Count the occurrence of each genotype
GT$Freq<-rowSums(!is.na(GT[-which(names(GT)=="Genotype")]))

#Remove the genotype with X
remove<-grep("X",GT$Genotype)
GT<-GT[-remove,]

GT<-GT[order(GT$Freq,decreasing = T),]
write.csv(GT[1:100,], "Output/AA/Common.Genotype10_top100.csv")

hist(GT$Freq)
GT2<-GT[GT$Freq>=10,] #86

GT2<-GT[GT$Freq>=20,] #46

#Genotypes appeared in more than 20 files
write.csv(GT2,"Output/AA/Mostcommon.Genotype10.csv")

ggplot(data=GT, aes(x=Freq))+
    geom_histogram(fill="lightblue", bins=30, color="gray60")+
    xlab("Genotype Freq out of 108 samples")+ylab('No. of genotypes')
ggsave("Output/AA/Genotype10.freq.pdf", width = 6, height = 3)    



##################################################################
### The first 7 AA sites for including the tissue
sites<-sites[1:7]

#AA sites
AAsites<-ceiling(sites/3)
#NT:297 305 326 328 340 373 431
#AA: 99 102 109 110 114 125 144

site.cols<-paste0("pos.", AAsites)

for (i in 1:length(AAs)){
    aa<-read.csv(paste0("Output/AA/",AAs[i]), row.names = 1,colClasses = "character" )
    #select the 10 sites of interest
    col.num <- which(colnames(aa) %in% site.cols)
    aa<-aa[,col.num]
    fname<-substr(AAs[i], 1,7)
    
    write.csv(aa,paste0("Output/AA/Genotype/", fname, ".AA.7variants.csv"))
    print(fname)
}


files<-list.files("Output/AA/Genotype/",pattern="AA.7variants.csv$")
files<-list.files("Output/AA/Genotype/",pattern="^Run8.+AA.7variants.csv$")

genotypeCount<-data.frame(file=sub(".AA.7variants.csv",'',files))

ref<-ref[ref$pos %in% sites,]
WTaa<-toupper(paste0(ref$RefAA, collapse = ''))


for (i in 1:length(files)){
    df<-read.csv(paste0("Output/AA/Genotype/",files[i]), stringsAsFactors = F, row.names = 1, colClasses = "character")
    fname<-substring(files[i],1,7)
    
    #calculate the freq of each genotype
    df$Genotype<-apply(df[1:ncol(df)], 1, function(x) paste0(x,collapse = ''))
    genotypeCount$Total[i]<-length(unique(df$Genotype))
    
    #remove the sequences with unknown aa  ('X')  
    df2<-df[!grepl("X", df$Genotype, fixed=T),]
    genotypeCount$no_Xs[i]<-length(unique(df2$Genotype))
    
    #Genotype freq. all
    gt<-data.frame(table(df$Genotype))
    gt<-gt[order(gt$Freq, decreasing = T),]
    gt$Var1<-as.character(gt$Var1)
    
    w<-which(gt$Var1==WTaa)
    #Replace the wildtype (Reference) with WT
    gt$Var1[w]<-"WT"
    genotypeCount$TotalReads[i]<-nrow(df)
    genotypeCount$WTcount[i]<-ifelse(length(gt$Freq[gt$Var1=="WT"])==0, 0,gt$Freq[gt$Var1=="WT"])
    
    
    s<-nrow(gt[gt$Freq==1,])
    gt[nrow(gt)+1,]<-c("Singleton",s)
    #remove singletons and #remove those with X
    gt2<-gt[gt$Freq!=1,] 
    gt2<-gt2[!grepl("X", gt2$Var1, fixed=T),]
    
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
    
    write.csv(gt3,paste0("Output/AA/Genotype.freq7/",fname,".AA.7gtypeFreq.csv"))
    print(fname)
}

genotypeCount$PercentSingleton<-genotypeCount$Singleton/genotypeCount$TotalReads*100
genotypeCount$PercentDoubleton<-genotypeCount$Doubleton/genotypeCount$TotalReads*100

write.csv(genotypeCount,"Output/AA/AA.Genotype7_count_summary.csv")



#### Find common genotypes across monkeys/samples

freqfiles<-list.files("Output/AA/Genotype.freq7/", pattern='7gtypeFreq.csv')


gtList<-list()

GT<-data.frame(Genotype="WT")
for (i in 1:length(freqfiles)){
    
    fname<-substring(freqfiles[i],1,7)
    #for (i in 1:length(gfiles)){
    df<-read.csv(paste0("Output/AA/Genotype.freq7/",freqfiles[i]), row.name=1,stringsAsFactors = F)
    colnames(df)[1]<-"Genotype"
    #select top 70 genotypes
    df<-df[!(df$Genotype=="Singleton"|df$Genotype=="Doubleton"),]
    df<-df[order(df$Freq, decreasing = T),]
    
    if (nrow(df)>70) df<-df[1:70,]
    
    GT<-merge(GT,df[,1:2], by="Genotype", all=T)
    colnames(GT)[i+1]<-fname
    
    gtList[[i]]<-GT
    names(gtList)[i]<-fname
    print(fname)
}

#Count the occurrence of each genotype
GT$Freq<-rowSums(!is.na(GT[-which(names(GT)=="Genotype")]))

GT<-GT[order(GT$Freq,decreasing = T),]
write.csv(GT[1:100,], "Output/AA/Common7Genotypes_top100.csv")

hist(GT$Freq)
GT2<-GT[GT$Freq>=10,] #71 genotypes

GT2<-GT[GT$Freq>=20,] #35 genotypes

#Genotypes appeared in more than 20 files
write.csv(GT2,"Output/AA/Mostcommon7Genotypes_Top35.csv")



ggplot(data=GT, aes(x=Freq))+
    geom_histogram(fill="lightblue", bins=30, color="gray60")+
    xlab("Genotype Freq out of 108 samples")+ylab('No. of genotypes')
ggsave("Output/AA/Genotype7.freq.pdf", width = 6, height = 3)    



