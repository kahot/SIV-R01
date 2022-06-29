#Genotype with 11 sites (table position 328 for now)
library(gridExtra)
library(colorspace)
library(seqinr)
library(Rsamtools)
library(GenomicAlignments)
source("Rscripts/baseRscript.R")
cols<-qualitative_hcl(6, palette="Dark3")


#Interesting 11 sites 
hsites<-read.csv("Output/MF_PID/HighFreq/Summary_highFreq_sites.csv", stringsAsFactors = F)
#Remove 328 
hsites<-hsites[hsites$pos!=328,]
sites<-hsites$pos
#297 305 326 340 373 431 485 515 536 539 548

#Alignments are between pos 215 and 691
pos<-sites-214

fastas<-list.files("Output/PID_Con_Alignment/")

for (i in 1:length(fastas)){
    align<-readDNAStringSet(paste0("Output/PID_Con_Alignment/",fastas[i]), format="fasta")
    highSites<-data.frame(pos=sites)
    align2<-as.list(as.character(align))
    
    a1<-lapply(align2,function(x) unlist(strsplit(x,'')))
    a2<-lapply(a1, function(x) x[pos])
    
    Sites1<-do.call(cbind.data.frame, a2)
    colnames(Sites1)<-1:ncol(Sites1)
    highSites<-cbind(highSites,Sites1)
    
    #Sites2<-do.call(rbind.data.frame, a2)
    #colnames(Sites2)<-sites
    fname<-gsub(".Alignment.fasta","",fastas[i])
    write.csv(highSites,paste0("Output/Genotype/", fname, ".11variants.csv"))
   
    print(fname)
}



## Calculate the genotype frequencies 
files<-list.files("Output/Genotype/", pattern="\\.11variants.csv")
#files<-list.files("Output/Genotype/", pattern="^Run8.*\\.11variants.csv")

#Ref genotype
ref<-read.csv("Output/Overview.ref.csv", row.names = 1,stringsAsFactors = F)
ref<-ref[ref$pos %in% sites,]
WT<-toupper(paste0(ref$ref, collapse = ''))


genotypeCount<-data.frame(file=sub(".11variants.csv",'',files))

#files<-files[files!="Run4_18.11variants.csv"]
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
    
    if (nrow(df2)==0){
        genotypeCount$TotalCleanReads[i]<-NA
        genotypeCount$WTcount[i]<-NA
        genotypeCount$Singleton[i]<-NA
    }
    #Genotype freq.
    else{
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
        write.csv(gt2,paste0("Output/Genotype/Freq/",fname,".11gtypeFreq.csv"))

    }
    print(fname)
}

genotypeCount$PercentWT<-genotypeCount$WTcount/genotypeCount$TotalCleanReads*100
genotypeCount$PercentSingleton<-genotypeCount$Singleton/genotypeCount$TotalCleanReads*100

write.csv(genotypeCount,"Output/Genotype/Genotype11_count_summary.csv")

gt<-genotypeCount[,c(1,9:10)]
gtm<-melt(gt)
gtm$variable<-factor(gtm$variable, levels=c("PercentSingleton","PercentWT"))
gtm<-gtm[order(gtm$variable),]
ggplot(gtm,aes(x=file,y=value, fill=variable))+
    geom_bar(position="stack",stat='identity')+
    scale_fill_manual(values=cols[c(5,1)], labels=c("% Singleton","% WT"))+
    theme_bw()+ylab('')+xlab('')+
    theme(axis.text.x = element_text(angle=90, size=5), legend.title = element_blank())
ggsave("Output/Genotype/Freq.of.WT.Singleton.11genotypes.pdf",height = 4,width = 10)    
    

