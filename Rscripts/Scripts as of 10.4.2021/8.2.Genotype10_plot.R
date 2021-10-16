library(ggplot2)
library(gridExtra)
library(colorspace)
library(streamgraph)
library(dplyr)
library(htmlwidgets)
source("Rscripts/baseRscript2.R")
cols<-qualitative_hcl(6, palette="Dark3")


## Genotype freq files
files<-list.files("Output/Genotype/Freq/", pattern="10gtypeFreq.csv")

Gtype<-list()
for (i in 1: length(files)){
    df<-read.csv(paste0("Output/Genotype/Freq/",files[i]), stringsAsFactors = F, row.names = 1)
    fname<-substring(files[i],1,7)
    Gtype[[i]]<-df
    names(Gtype)[i]<-fname
}
    

#create stream plot    
#library(ggstream)

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
monkeys<-names(monkeyList)
monkeys2<-monkeyList[tbs$ids]


#Stock genotypes
stock<-Gtype[["Run0_17"]]
d<-nrow(stock[stock$Freq==2,])
stock<-stock[stock$Freq!=2,]
stock[nrow(stock)+1,]<-c("Doubleton",d,"Run0_17")
stock$Freq<-as.integer(stock$Freq)
stock$Prop<-stock$Freq/sum(stock$Freq)
stock$Week<-0
stock$Tissue<-"Stock"

sites<-c(305,328,340,368,373,431,475,515,536,539)
pos<-sites-189

for (m in 1:length(monkeys2)){
    #Select the monkey
    sample<-monkeys2[[m]]
    sample<-sample[order(sample$Week),]
    gfiles<-Gtype[as.vector(sample$File.name)]
    monkey<-names(monkeys2)[m]
    tbweek<-tbs$tb[tbs$ids==monkey]
    
    
    #Createa a vector of all unique genotypes for each monkey
    #Eliminate singletons and doubletons
    gtypes<-c(unique(stock$Var1))
    DF<-list()
    DF[[1]]<-stock
    for (j in 1:length(gfiles)){
        #remove the doubletons
        df2<-gfiles[[j]]
        #count # of the doubletons and group them 
        d<-nrow(df2[df2$Freq==2,])
        df2<-df2[df2$Freq!=2,]
        df2[nrow(df2)+1,]<-c("Doubleton",d, paste(df2$ID[1]))
        df2$Freq<-as.integer(df2$Freq)
        df2$Prop<-df2$Freq/sum(df2$Freq)
        df2$Week<-sample$Week[j]
        df2$Tissue<-paste0(sample$Tissue2[j])
        
        gtypes<-c(gtypes,unique(df2$Var1))
        DF[[(j+1)]]<-df2
    }
    gtypes<-unique(gtypes) #161 (+singletons and doubletons) unique genotypes for 1st monkey (21918)
    
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
    
    gtorder<-unique(plasma$Genotype)
    gtorder<-gtorder[!(gtorder=="WT"|gtorder=="Singleton"|gtorder=="Doubleton")]
    plasma$Genotype<-factor(plasma$Genotype, levels=c("Singleton","Doubleton",gtorder,"WT"))
    #print legend
    
    #ggplot(data=plasma,aes(x=Week, y=Prop, fill=Genotype))+
    #    geom_area()+
    #    geom_vline(xintercept=tbweek, col="blue", size=0.3)+ylab("Proportion")+
    #    ggtitle(paste0(monkey,' Plasma'))+
    #    theme_bw()
    #ggsave(paste0("Output/Genotype/Freq/AreaStacked.",monkey,".10Genotype.Legends.pdf"), width =15,height = 5)
    legends <- c("WT","Singleton","Doubleton")
    ggplot(data=plasma,aes(x=Week, y=Prop, fill=Genotype))+
        geom_area(color="gray80",size=0.05)+
        geom_vline(xintercept=tbweek, col="blue", size=0.3)+ylab("Proportion")+
        ggtitle(paste0(monkey,' Plasma'))+
        theme_bw()+
        theme(legend.title = element_blank())+
        scale_fill_discrete(breaks=legends)
    ggsave(paste0("Output/Genotype/Freq/AreaStacked.",monkey,".10Genotype.pdf"),width = 9,height = 5)
    
    
    #stacked bar chart of all files      
    Gdata$label<-paste0("Wk ",Gdata$Week," ",Gdata$Tissue)
    L<-unique(Gdata$label)
    Gdata$label<-factor(Gdata$label,levels=L)
    ggplot()+
        geom_bar(data=Gdata,aes(x=label, y=Prop, fill=Genotype), stat="identity")+
        ylab("Proportion")+
        theme(legend.position = "none",axis.title.x = element_blank(),axis.text.x = element_text(angle=45, hjust=1))
    ggsave(paste0("Output/Genotype/Freq/BarStacked.",monkey,".10Genotype.pdf"),width = 12,height = 5)
    
    #interactive streamhrgraph
    plasma$Date<-as.Date(plasma$Week*7, origin="2020-01-01")
    pp<-streamgraph(plasma, key="Genotype", value="Prop", date="Date", height="300px", width="700px")%>%
        sg_axis_x(tick_interval = 1,tick_unit="month",tick_format ="%m%d")%>%
        sg_legend(show=TRUE, label="Genotype: ")
    
    saveWidget(pp,file=paste0("Output/Genotype/Freq/Streamplot/.",monkey,".plasama.html"))
    
}


##adding a legend for certain types for 1 monkey 



## Look for the common gneotypes among monkeys?

gtList<-list()
commonGt<-list()
for (m in 1:length(monkeys2)){
    sample<-monkeys2[[m]]
    sample<-sample[order(sample$Week),]
    sample<-sample[sample$File.name!="Run4_18",]
    gfiles<-Gtype[as.vector(sample$File.name)]
    monkey<-names(monkeys2)[m]
    tbweek<-tbs$tb[tbs$ids==monkey]
    
    GT<-data.frame(Genotype=c("WT"))
    for (i in 1:length(gfiles)){
        df<-gfiles[[i]]
        colnames(df)[1]<-"Genotype"
        GT<-merge(GT,df[,1:2], all=T)
        colnames(GT)[i+1]<-paste0(names(gfiles[i]),".",sample$Week[i])
    }
    gtList[[m]]<-GT
    names(gtList)[m]<-monkey
    
    
    #Createa a vector of top 20 genotypes from all weeks for each monkey 
    gt<-c()
    for (i in 1:length(gfiles)){
        df<-gfiles[[i]]
        colnames(df)[1]<-"Genotype"
        df<-df[!(df$Genotype=="Singleton"|df$Genotype=="Doubleton"),]
        df<-df[order(df$Freq, decreasing = T),]
        gt<-c(gt, df$Genotype[1:50])
    }
    gt2<-unique(gt)
    comm<-data.frame(Genotype=gt2[1:50])
    comm[paste(monkey)]<-"Y"
    commonGt[[m]]<-comm
    names(commonGt)[m]<-monkey
    
}

ref<-read.csv("Output/Overview.ref.csv", row.names = 1,stringsAsFactors = F)
ref<-ref[ref$pos %in% sites,]
WT<-toupper(paste0(ref$ref, collapse = ''))


for (i in 1:length(gtList)){
    df<-gtList[[i]]
    #number of available files for the monkey
    n<-ncol(df)-1
    #replace WT with seq.
    
    #remove singleton 
    df<-df[df$Genotype!="Singleton",]
    
    #replace WT with actual bases
    df$Genotype[df$Genotype=="WT"]<-WT
    
    df2<-data.frame(do.call(rbind, str_split(df$Genotype, '')))
    colnames(df2)<-sites
    df2$Genotype<-df$Genotype
    
    df2
    
   
    for (j in 1:n){
        dt<-data.frame(cbind(df[,c(j+1,j+2)], df2))
        dt2<-dt[!is.na(dt[,1])|!is.na(dt[,2]),]
        

    }

}



Common<-commonGt[[1]]

for (m in 2:length(commonGt)){
    df<-commonGt[[m]]
    Common<-merge(Common,df, by="Genotype", all=T )
}
Common$Total<-apply(Common[2:11], 1, function(x) length(which(x=="Y")) )

table(Common$Total)
#1  2  3  4  5  6  7  8  9 10 
#62 24  8  4 12 74  3  8 12  2 

Com<-Common[Common$Total!=1,]
Com$Genotype[Com$Total%in% c(8,9,10)]
#[1] "ACAAACAGAC" "ACAAGAACAC" "ACAAGAAGAA" "ACAAGAAGAC" "ACAAGAAGGA"
#[6] "ACAAGCAGAC" "ACAAGGAGAC" "ACGAAAACAC" "ACGAAAAGAC" "ACGAAAAGGA"
#[11] "ACGAACAGAC" "ACGAAGAGAC" "GCAAGAAGGA" "GCGAAAACAC" "GCGAAAAGAC"
#[16] "GCGAAAAGGA" "GCGAACAGGA" "GCGAAGAGAC" "GCGAAGAGGA" "GCGAGAAGAC"
#[21] "GCGAGAAGGA" "WT"       
write.csv(Com,"Output/Genotype/MostCommonGenotype_10genotypes.csv")


G<-read.csv("Output/Genotype/MostCommonGenotype_10genotypes.csv", stringsAsFactors = F, row.names = 1)
do.call(rbind, str_split(before$type, '_and_'))
