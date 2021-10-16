#Genotype with 10 sites
library(ggplot2)
library(gridExtra)
library(colorspace)
library(seqinr)
library(Rsamtools)
library(GenomicAlignments)
library(stringr)
source("Rscripts/baseRscript2.R")
cols<-qualitative_hcl(6, palette="Dark3")

#List the bam files where merged PID processed reads are mapped 
bams<-list.files("Output/bam_PID/",pattern=".bam$")
#bams<-list.files("Output/bam_PID/",pattern="^Run7_.*bam$")

#parameter setting for stackStringsFromBam
gr <- GRanges(seqnames="AY032751env",  IRanges(190, 691))

#Interesting 10 sites 
sites<-c(305,328,340,368,373,431,475,515,536,539)
pos<-sites-189

for (i in 1:length(bams)){
    index<-paste0("Output/bam_PID/",bams[i],'.bai')
    bamfile <- BamFile(paste0("Output/bam_PID/",bams[i]), index=index)
    
    fname<-substr(bams[i], 1,7)
    align<-stackStringsFromBam(bamfile, param=gr)
    

    highSites<-data.frame(pos=sites)
    align2<-as.list(as.character(align))
    
    a1<-lapply(align2,function(x) unlist(strsplit(x,'')))
    a2<-lapply(a1, function(x) x[pos])
    
    Sites1<-do.call(cbind.data.frame, a2)
    colnames(Sites1)<-1:ncol(Sites1)
    highSites<-cbind(highSites,Sites1)
    
    #Sites2<-do.call(rbind.data.frame, a2)
    #colnames(Sites2)<-sites
    
    write.csv(highSites,paste0("Output/Genotype/", fname, ".10variants.csv"))
    #write.csv(Sites2,paste0("Output/Genotype/Var/", fname, "variantsMatrix.csv"))
    
    print(fname)
}


## Calculate the genotype frequencies 
files<-list.files("Output/Genotype/", pattern="\\.10variants.csv")

#Ref genotype
ref<-read.csv("Output/Overview.ref.csv", row.names = 1,stringsAsFactors = F)
ref<-ref[ref$pos %in% sites,]
WT<-toupper(paste0(ref$ref, collapse = ''))


genotypeCount<-data.frame(file=sub(".10variants.csv",'',files))

files<-files[files!="Run4_18.10variants.csv"]
#remove Run4_18: all genotypes have "N", remove from the analysis for now.
for (i in 1: length(files)){
    df<-read.csv(paste0("Output/Genotype/",files[i]), stringsAsFactors = F, row.names = 1)
    fname<-substring(files[i],1,7)
    
    #calculate the freq of each genotype
    df2<-data.frame(seq=colnames(df)[2:ncol(df)])
    df2$Genotype<-apply(df[2:ncol(df)], 2, function(x) paste0(x,collapse = ''))
    genotypeCount$Total[i]<-length(unique(df2$Genotype))
    
    #remove the sequences with gpas ('+')  
    df2<-df2[!grepl("+", df2$Genotype, fixed=T),]
    genotypeCount$no_gaps[i]<-length(unique(df2$Genotype))
    
    #remove the sequences with N for now
    df2<-df2[!grepl("N", df2$Genotype, fixed=T),]
    genotypeCount$no_Ns[i]<-length(unique(df2$Genotype))
    
    #remove the sequences with '-' for now
    df2<-df2[!grepl("-", df2$Genotype, fixed=T),]
    genotypeCount$no_missing[i]<-length(unique(df2$Genotype))
    
    #Genotype freq.
    gt<-data.frame(table(df2$Genotype))
    gt<-gt[order(gt$Freq, decreasing = T),]
    gt$Var1<-as.character(gt$Var1)
    
    w<-which(gt$Var1==WT)
    #Replace the wildtype (Reference) with WT
    gt$Var1[w]<-"WT"
    genotypeCount$TotalCleanReads[i]<-nrow(df2)
    genotypeCount$WTcount[i]<-ifelse(length(gt$Freq[gt$Var1=="WT"])==0, 0,gt$Freq[gt$Var1=="WT"])
    
    #remove singletons
    gt2<-gt[gt$Freq!=1,] 
    s<-nrow(gt[gt$Freq==1,])
    gt2[nrow(gt2)+1,]<-c("Singleton",s)
    gt2$Freq<-as.integer(gt2$Freq)
    gt2$ID<-fname
    genotypeCount$Singleton[i]<-s
    
    write.csv(gt2,paste0("Output/Genotype/Freq/",fname,".10gtypeFreq.csv"))
    print(fname)
}



genotypeCount$PercentWT<-genotypeCount$WTcount/genotypeCount$TotalCleanReads*100
genotypeCount$PercentSingleton<-genotypeCount$Singleton/genotypeCount$TotalCleanReads*100

write.csv(genotypeCount,"Output/Genotype/Genotype10_count_summary.csv")

gt<-genotypeCount[,c(1,8:9)]
gtm<-melt(gt)
gtm$variable<-factor(gtm$variable, levels=c("PercentSingleton","PercentWT"))
gtm<-gtm[order(gtm$variable),]
ggplot(gtm,aes(x=file,y=value, fill=variable))+
    geom_bar(position="stack",stat='identity')+
    scale_fill_manual(values=cols[c(5,1)], labels=c("% Singleton","% WT"))+
    theme_bw()+ylab('')+xlab('')+
    theme(axis.text.x = element_text(angle=90, size=5), legend.title = element_blank())
ggsave("Output/Genotype/Freq.of.WT.Singleton.10genotypes.pdf",height = 3,width = 8)    
    

