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

#Interesting sites
sites<-c(305,326,328,340,368,373,431,475,485,500,515,536,539,548)

## Genotype freq files
files<-list.files("Output/Genotype/Var/", pattern="variantsMatrix.csv")

Var<-list()
for (i in 1: length(files)){
    df<-read.csv(paste0("Output/Genotype/Var/",files[i]), stringsAsFactors = F, row.names = 1)
    fname<-substring(files[i],1,7)
    Var[[i]]<-df
    names(Var)[i]<-fname
}

#Run4_20 has too many N to interpret at this time    

#Read the files by monkey
tbs<-read.csv("Data/TBinfection.week.csv",stringsAsFactors = F)
samples<-read.csv("Data/SamplesNoduplicates.csv",stringsAsFactors = F)
samples$Week<-as.integer(samples$Week)
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

monkeyList<-monkeyList[tbs$ids]
monkeys<-names(monkeyList)

#Start with stock
#Stock genotypes
stock<-Var[["Run0_17"]]
colnames(stock)<-sites

#work on up to first five sites
basecomb<-paste0("sites.1.",2:5)

#paste0(basecomb, "=character()")
stock2<-data.frame(seq=1:nrow(stock))

for (i in 1: length(basecomb)){
    endcol<-i+1
    stock2[basecomb[i]]<-apply(stock[1:endcol], 1, function(x) paste0(x,collapse = ''))
}
#calculate frequency of each nucleotide combination and save in a list
StockFq<-list()
StockCounts<-list()
for (i in 1: length(basecomb)){
    df<-stock2[,c(1,(i+1))]
    #remove the sequences with gpas ('+')  
    df<-df[!grepl("+", df[,2], fixed=T),]
    #remove the sequences with N for now
    df<-df[!grepl("N", df[,2], fixed=T),]
    #remove the sequences with -N for now
    df<-df[!grepl("-", df[,2], fixed=T),]
    #Frequency
    fq<-data.frame(table(df[,2]))
    fq$Stock<-fq$Freq/sum(fq$Freq)

    StockFq[[i]]<-fq[,c(1,3)]
    names(StockFq)[i]<-basecomb[i]
    StockCounts[[i]]<-fq[,c(1,2)]
    names(StockCounts)[i]<-basecomb[i]
}

# Try calculating paris of nucleotide freq.
stock3<-data.frame(seq=1:nrow(stock))
basepair<-paste0("sites.",1:4,".",2:5)

table(stock$`305`)
table(stock$`326`)
pair<-StockCounts[[1]]
pair$Prop<-pair$Freq/sum(pair$Freq)

for (i in 1: length(basepair)){
    stock3[basepair[i]]<-apply(stock[i:(i+1)], 1, function(x) paste0(x,collapse = ''))
}
#calculate frequency of each nucleotide combination and save in a list
StockFq2<-list()
for (i in 1: length(basepair)){
    df<-stock3[,c(1,(i+1))]
    #remove the sequences with gpas ('+')  
    df<-df[!grepl("+", df[,2], fixed=T),]
    #remove the sequences with N for now
    df<-df[!grepl("N", df[,2], fixed=T),]
    #remove the sequences with -N for now
    df<-df[!grepl("-", df[,2], fixed=T),]
    #Frequency
    fq<-data.frame(table(df[,2]))
    fq$Stock<-fq$Freq/sum(fq$Freq)
    fq1<-fq[,c(1,3)]
    StockFq2[[i]]<-fq[,c(1,3)]
    names(StockFq2)[i]<-basepair[i]
    StockCounts[[i]]<-fq[,c(1,2)]
    names(StockCounts)[i]<-basepair[i]
}




#For monkey1
#for (m in 1:length(monkeys2)){
#Select the monkey
m=1
    sample<-monkeys2[[m]]
    #sample<-sample[sample$Tissue,]
    sample<-sample[order(sample$Week),]
    varfiles<-Var[as.vector(sample$File.name)]
    monkey<-names(monkeyList)[m]
    #tbweek<-tbs$tb[tbs$ids==monkey]
    
    #Stock's information
    for (k in 1:length(basecomb)){
        fqs<-StockFq[[basecomb[k]]]
        counts<-StockCounts[[basecomb[k]]]
        for (j in 1: length(varfiles)){
            cname<-paste0("Week.",sample$Week[j]," ",sample$Tissue2[j])
            variant<-varfiles[[j]]
    
            df2<-data.frame(seq=1:nrow(variant))
            df2[,cname]<-apply(variant[1:(k+1)], 1, function(x) paste0(x,collapse = ''))
            df2<-df2[!grepl("+", df2[,2], fixed=T),]
            df2<-df2[!grepl("N", df2[,2], fixed=T),]
            df2<-df2[!grepl("-", df2[,2], fixed=T),]
            #Frequency
            fq2<-data.frame(table(df2[,2]))
            fq2[,cname]<-fq2$Freq/sum(fq2$Freq)
    
            fqs<-merge(fqs, fq2[,c("Var1",cname)], by="Var1", all=T)
            fqs[is.na(fqs)]<-0
            fq2<-fq2[,1:2]
            colnames(fq2)[2]<-cname
            counts<-merge(counts,fq2[,c("Var1",cname)],by="Var1", all=T)
            
        }
        
        counts[is.na(counts)]<-0
        write.csv(fqs,paste0("Output/Genotype/Var/Freq/", monkey, ".",basecomb[k],".csv"))
        
        #follow the WT GTC
        ct<-counts[counts$Var1=="GTC"|counts$Var1=="GTG"|counts$Var1=="GTT"|counts$Var1=="GTA", ]
        
        ctp<-data.frame(apply(ct[,-1],2,function(x) x/sum(x)*100 ))
        colnames(ctp)[1]<-"Stock"
        ctp$Nuc<-ct$Var1
        ctpm<-melt(ctp)
        ggplot(ctpm,aes(x=variable,y=value, fill=Nuc))+
            geom_bar(stat="identity")+
            theme_bw()+ylab("Genotype proportion")+
            theme(legend.title = element_blank())+
            ggtitle(paste0("pos 305 to ", sites[k+1]))+
            theme(axis.text.x = element_text(angle=90, hjust=1), axis.title.x = element_blank())
        
        
        
        
        
        fqsm<-melt(fqs)
        ggplot(fqsm,aes(x=variable,y=value, fill=Var1))+
            geom_bar(stat="identity")+
            theme_bw()+ylab("Genotype proportion")+
            theme(legend.title = element_blank())+
            ggtitle(paste0("pos 305 to ", sites[k+1]))+
            theme(axis.text.x = element_text(angle=90, hjust=1), axis.title.x = element_blank())
        ggsave(paste0("Output/Genotype/Var/Figures/", basecomb[k],".",monkey,".pdf"), width=6, height = 4)
    }

    
    #Pairs
    for (k in 1:length(basepair)){
        fqs<-StockFq2[[basepair[k]]]
        
        for (j in 1: length(varfiles)){
            cname<-paste0("Week.",sample$Week[j]," ",sample$Tissue2[j])
            variant<-varfiles[[j]]
            
            df2<-data.frame(seq=1:nrow(variant))
            df2[,cname]<-apply(variant[k:(k+1)], 1, function(x) paste0(x,collapse = ''))
            df2<-df2[!grepl("+", df2[,2], fixed=T),]
            df2<-df2[!grepl("N", df2[,2], fixed=T),]
            df2<-df2[!grepl("-", df2[,2], fixed=T),]
            #Frequency
            fq2<-data.frame(table(df2[,2]))
            fq2[,cname]<-fq2$Freq/sum(fq2$Freq)
            
            fqs<-merge(fqs, fq2[,c("Var1",cname)], by="Var1", all=T)
            fqs[is.na(fqs)]<-0
        }
        
        write.csv(fqs,paste0("Output/Genotype/Var/Freq/", monkey, ".",basepair[k],".csv"))
        fqsm<-melt(fqs)
        ggplot(fqsm,aes(x=variable,y=value, fill=Var1))+
            geom_bar(stat="identity")+
            theme_bw()+ylab("Genotype proportion")+
            theme(legend.title = element_blank())+
            ggtitle(paste0("pos 305 to ", sites[k+1]))+
            theme(axis.text.x = element_text(angle=90, hjust=1), axis.title.x = element_blank())
        ggsave(paste0("Output/Genotype/Var/Figures/",monkey,".",k,".", basepair[k],".pdf"), width=6, height = 4)
    }
    
    
    
    
    #for (i in 1: length(basecomb)){
    #        endcol<-i+1
    #        Seq[basecomb[i]]<-apply(variant[1:endcol], 1, function(x) paste0(x,collapse = ''))
    #    }
    #    
    for (i in 1:ncol(Seq)){
        df2<-Seq[,c(1,(i+1))]
        
        #remove the sequences with gpas ('+')  
        df2<-df2[!grepl("+", df2[,2], fixed=T),]
        #remove the sequences with N for now
        df2<-df2[!grepl("N", df2[,2], fixed=T),]
        #remove the sequences with -N for now
        df2<-df2[!grepl("-", df2[,2], fixed=T),]
        #Frequency
        fq2<-data.frame(table(df[,2]))
        fq2[,week]<-fq2$Freq/sum(fq2$Freq)
        
        df<-stock2[,c(1,(i+1))]
        #remove the sequences with gpas ('+')  
        df<-df[!grepl("+", df[,2], fixed=T),]
        #remove the sequences with N for now
        df<-df[!grepl("N", df[,2], fixed=T),]
        #remove the sequences with -N for now
        df<-df[!grepl("-", df[,2], fixed=T),]
        #Frequency
        fq<-data.frame(table(df[,2]))
        fq$Stock<-fq$Freq/sum(fq$Freq)
        
        fqs<-merge(fq[,c("Var1","Stock")], fq2[,c("Var1",week)], by="Var1", all=T)
        fqs[is.na(fqs)]<-0
        fqsm<-melt(fqs)
        ggplot(fqsm,aes(x=variable,y=value, fill=Var1))+
            geom_bar(stat="identity")
        
 }   
    

#Calculate D (linkage disequilibrium )
pa=0.7518487
pb=0.760883
D=0.1789225
D^2/pa*(1-pa)*pb*(1-pb)
pa*(1-pb)
pb*(1-pa)
