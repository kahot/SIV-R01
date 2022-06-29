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

#Interesting 11 sites 
hsites<-read.csv("Output/MF_PID/HighFreq/HighFreq.sites.selected.csv", stringsAsFactors = F)
sites<-hsites$pos

#AA sites
AAsites<-ceiling(sites/3)
#NT:297 305 326 340 373 431 485 515 536 539 548 
#AA: 99 102 109 114 125 144 162 172 179 180 183 

site.cols<-paste0("pos.", AAsites)

for (i in 1:length(AAs)){
    aa<-read.csv(paste0("Output/AA/",AAs[i]), row.names = 1,colClasses = "character" )
    #select the 10 sites of interest
    col.num <- which(colnames(aa) %in% site.cols)
    aa<-aa[,col.num]
    fname<-substr(AAs[i], 1,7)
    
    write.csv(aa,paste0("Output/AA/Genotype/", fname, ".AA.11variants.csv"))
    print(fname)
}


## Calculate the AA genotype frequencies 
files<-list.files("Output/AA/Genotype/", pattern="\\.11variants.csv")

#Ref genotype
ref<-read.csv("Output/Overview.ref.csv", row.names = 1,stringsAsFactors = F,colClasses = "character" )
#select the 10 sites of interest)
ref<-ref[ref$pos %in% sites,]
WTaa<-toupper(paste0(ref$RefAA, collapse = ''))

#exclude Run4_18 for now (too many missing sites)
files<-files[files!="Run4_18.AA.11variants.csv"]
genotypeCount<-data.frame(file=sub(".AA.11variants.csv",'',files))

for (i in 1:length(files)){
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
    
    write.csv(gt3,paste0("Output/AA/Genotype.freq/",fname,".AA.11gtypeFreq.csv"))
    print(fname)
}



genotypeCount$PercentWT<-genotypeCount$WTcount/genotypeCount$TotalCleanReads*100
genotypeCount$PercentSingleton<-genotypeCount$Singleton/genotypeCount$TotalCleanReads*100
genotypeCount$PercentDoubleton<-genotypeCount$Doubleton/genotypeCount$TotalCleanReads*100

write.csv(genotypeCount,"Output/AA/AA.Genotype11_count_summary.csv")

###
### Run6_18 only has 16 sequences -do not calculate singleton and doubleton for this file

df<-read.csv("Output/AA/Genotype/Run6_18.AA.11variants.csv", stringsAsFactors = F, row.names = 1, colClasses = "character")
fname<-"Run6_18"
df$Genotype<-apply(df[1:ncol(df)], 1, function(x) paste0(x,collapse = ''))
df<-df[!grepl("X", df$Genotype, fixed=T),]
#Genotype freq.
gt<-data.frame(table(df$Genotype))
gt<-gt[order(gt$Freq, decreasing = T),]
gt$Var1<-as.character(gt$Var1)

w<-which(gt$Var1==WTaa)
gt$Var1[w]<-"WT"
gt$ID<-fname
write.csv(gt, "Output/AA/Genotype.freq/Run6_18.AA.11gtypeFreq.csv")

######


gt<-genotypeCount[,c(1,8:10)]
gtm<-melt(gt)
gtm$variable<-factor(gtm$variable, levels=c("PercentDoubleton","PercentSingleton","PercentWT"))
gtm<-gtm[order(gtm$variable),]
ggplot(gtm,aes(x=file,y=value, fill=variable))+
    geom_bar(position="stack",stat='identity')+
    scale_fill_manual(values=paste0(cols[c(2,4,1)],"CC"), labels=c("% Doubleton","% Singleton","% WT"))+
    theme_bw()+ylab('')+xlab('')+
    theme(axis.text.x = element_text(angle=90, size=5), legend.title = element_blank())
ggsave("Output/AA/Proprtion.Single.Double.11genotypes.pdf",height = 3,width = 8)    



#### Find common genotypes across monkeys/samples

freqfiles<-list.files("Output/AA/Genotype.freq/", pattern='11gtypeFreq.csv')


gtList<-list()

GT<-data.frame(Genotype="WT")
for (i in 1:length(freqfiles)){

    fname<-substring(freqfiles[i],1,7)
    #for (i in 1:length(gfiles)){
    df<-read.csv(paste0("Output/AA/Genotype.freq/",freqfiles[i]), row.name=1,stringsAsFactors = F)
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
GT$Count<-rowSums(!is.na(GT[-which(names(GT)=="Genotype")]))

GT<-GT[order(GT$Count,decreasing = T),]

write.csv(GT, "Output/AA/CommonGenotypes_top70.csv")

hist(GT$Count)
GT2<-GT[GT$Count>=10,] #100 genotypes

GT2<-GT[GT$Count>=20,] #54 genotypes

#Top40
GT3<-GT[1:40,]
write.csv(GT3,"Output/AA/MostcommonGenotypes_40.csv")



ggplot(data=GT, aes(x=Count))+
    geom_histogram(fill="lightblue", bins=30, color="gray60")+
    xlab("Genotype appearance out of 109 samples")+ylab('')
ggsave("Output/AA/Genotype11.freq.pdf", width = 6, height = 3)    

ggplot(data=GT2, aes(x=Count))+
    geom_histogram(fill="lightblue", bins=30, color="gray60")+
    xlab("Genotype appearance out of 109 samples")+ylab('')
ggsave("Output/AA/Genotype11.freq.zoom.pdf", width = 5, height = 3)    


