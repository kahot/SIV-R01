#Genotype with 3 sites
library(ggplot2)
library(gridExtra)
library(colorspace)
library(seqinr)
library(stringr)
library(reshape2)
library(ggpattern)
library(plotrix)
cols<-qualitative_hcl(6, palette="Dark3")

### 

#### Read sample info files
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
monkeyLis<-monkeyList[tbs$ids]
monkeys<-names(monkeyLis)


########  Look at associations between P144Q and other two mutations

Ps<-c("pos.99","pos.102","pos.144")
AApatterns<-expand.grid(c("D","E"),c("R","K"),c("P","Q"))
AAcombs<-apply(AApatterns, 1, function(x) paste0(x,collapse = ''))
mutAA<-c("E","K","Q")

#Reference AAs
ref<-read.csv("Output/Overview.ref.csv", row.names = 1,stringsAsFactors = F,colClasses = "character" )
ref<-ref[ref$pos %in% c("297", "305","431"),]
WTaa<-toupper(paste0(ref$RefAA, collapse = ''))



## Genotyping by 3 sites
files<-list.files("Output/AA/Genotype.freq/", pattern="3gtypeFreq.csv")

Gtype<-list()
for (i in 1: length(files)){
  df<-read.csv(paste0("Output/AA/Genotype.freq/",files[i]), stringsAsFactors = F, row.names = 1)
  fname<-substring(files[i],1,7)
  Gtype[[i]]<-df
  names(Gtype)[i]<-fname
}


stock<-read.csv("Output/AA/Genotype.freq/Run0_17.AA.3gtypeFreq.csv", stringsAsFactors = F, row.names = 1)
stock$Freq<-as.integer(stock$Freq)
stock$Prop<-stock$Freq/sum(stock$Freq)
stock$Week<-0
stock$Tissue<-"Stock"

####
plots<-list()
plots2<-list()
for (m in 1:length(monkeyLis)){
  #Select the monkey
    sample<-monkeyLis[[m]]
    sample<-sample[order(sample$Week),]
    gfiles<-Gtype[as.vector(sample$File.name)]
    monkey<-names(monkeyLis)[m]
    tbweek<-tbs$tb[tbs$ids==monkey]
    
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
    
    #which genotypes have the specific mutation?
    Gdata$Genotype<-as.character(Gdata$Genotype)
    Gdata$Genotype[Gdata$Genotype=="WT"]<-WTaa
    
    #association between 144Q & other mutations
    Gdata$p144q<-sapply(Gdata$Genotype, function(x) {
        if (!(x=="Singleton"|x=="Doubleton")) {x=unlist(strsplit(x,"", fixed=T))
        ifelse (x[3]==mutAA[3], "Y","N")}
        else "N"})
    Gdata$d99e<-sapply(Gdata$Genotype, function(x) {
        if (!(x=="Singleton"|x=="Doubleton")) {x=unlist(strsplit(x,"", fixed=T))
        ifelse (x[1]==mutAA[1], "Y","N")}
        else "N"})
    Gdata$r102k<-sapply(Gdata$Genotype, function(x) {
        if (!(x=="Singleton"|x=="Doubleton")) {x=unlist(strsplit(x,"", fixed=T))
        ifelse (x[2]==mutAA[2], "Y","N")}
        else "N"})
    
    #Proportion of d99e without p144q
    Gdata$id<-paste0(Gdata$Week,".",Gdata$Tissue)
    
    ids<-unique(Gdata$id)
    assoc<-data.frame(id=ids)
    assoc$week<-c(0,sample$Week)
    assoc$tissue<-c("Stock", sample$Tissue2)
    assoc2<-assoc
    for (w in 1:length(ids)){
        df<-Gdata[Gdata$id==ids[w],]
        assoc$D99E[w]<-sum(df$Freq[df$d99e=="Y"&df$p144q=="N"])/sum(df$Freq[df$p144q=="N"])*100
        assoc$R102K[w]<-sum(df$Freq[df$r102k=="Y"&df$p144q=="N"])/sum(df$Freq[df$p144q=="N"])*100
        #assoc$D99E[w]<-sum(df$Freq[df$d99e=="Y"&df$p144q=="N"])/sum(df$Freq)*100
        #assoc$R102K[w]<-sum(df$Freq[df$r102k=="Y"&df$p144q=="N"])/sum(df$Freq)*100
    }
    assoc$D99E[is.nan(assoc$D99E)]<-0
    assoc$R102K[is.nan(assoc$R102K)]<-0
    
    for (w in 1:length(ids)){
        df<-Gdata[Gdata$id==ids[w],]
        assoc2$D99E[w]<-sum(df$Freq[df$d99e=="N"&df$p144q=="Y"])/sum(df$Freq[df$p144q=="Y"])*100
        assoc2$R102K[w]<-sum(df$Freq[df$r102k=="N"&df$p144q=="Y"])/sum(df$Freq[df$p144q=="Y"])*100
        
        #assoc2$D99E[w]<-sum(df$Freq[df$d99e=="N"&df$p144q=="Y"])/sum(df$Freq)*100
        #assoc2$R102K[w]<-sum(df$Freq[df$r102k=="N"&df$p144q=="Y"])/sum(df$Freq)*100
    }
    assoc2$D99E[is.nan(assoc2$D99E)]<-0
    assoc2$R102K[is.nan(assoc2$R102K)]<-0
    
    assocm<-melt(assoc, id.vars = c("id","week","tissue"))
    assocm$id<-factor(assocm$id, levels=paste(assoc$id))
    colnames(assocm)[4]<-"Mutation"
    if (max(assocm$value)<50) ymax=50
    else ymax=max(assocm$value)
    
    plots[[m]]<-ggplot(data=assocm, aes(x=id, y=value, fill=Mutation))+
        geom_bar(stat="identity", position = position_dodge(width = .8),alpha=0.8)+
        theme(axis.text.x = element_text(angle=90, hjust=1))+
        ggtitle(monkey)+xlab('')+
        scale_fill_manual(values = cols[c(1,7)])+
        ylab("% mut w/o P144Q/nonP144Q")+ylim(0,ymax)
    
   
    assocm2<-melt(assoc2, id.vars = c("id","week","tissue"))
    assocm2$id<-factor(assocm2$id, levels=paste(assoc2$id))
    colnames(assocm2)[4]<-"Mutation"
    if (max(assocm2$value)<50) ymax=50
    else ymax=max(assocm2$value)
    
    plots2[[m]]<-ggplot(data=assocm2, aes(x=id, y=value, fill=Mutation))+
        geom_bar(stat="identity", position = position_dodge(width = .8), alpha=0.8)+
        theme(axis.text.x = element_text(angle=90, hjust=1))+
        ggtitle(monkey)+xlab('')+
        ylab("% P144Q w/o other mut/allP144Q")+ylim(0,ymax)+
        scale_fill_manual(values = cols[c(1,7)], labels=c("no D99E","no R102K"))
}


pdf(paste0("Output/AA/Analysis/NoP144Q_withOtherMutation2.pdf"), width = 7, height =30)
do.call(grid.arrange, c(plots, ncol=1))
dev.off()

#Find the opposite -P144Q present but not others
pdf(paste0("Output/AA/Analysis/P144Qpresent_noOthermutation2.pdf"), width = 7, height =30)
do.call(grid.arrange, c(plots2, ncol=1))
dev.off()








#######

## Create stacked area charts with 3genotypes
stock<-read.csv("Output/AA/Genotype.freq/Run0_17.AA.3gtypeFreq.csv", stringsAsFactors = F, row.names = 1)
#Add ERQ to keep the color same for all monkeys
stock[12,]<-c("ERQ", 0)
stock$Freq<-as.integer(stock$Freq)
stock$Prop<-stock$Freq/sum(stock$Freq)
stock$Week<-0
stock$Tissue<-"Stock"


cols<-qualitative_hcl(8, palette="Dark3")
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
color<-color[color!="white"]

colors<-c("gray", cols, color)
aalabs<-c("D99E","R102K","P144Q")

#MUT<-data.frame(Monkey=monkeys)
#Initial<-data.frame(Monkey=monkeys)


for (m in 1:length(monkeyLis)){
    #Select the monkey
    sample<-monkeyLis[[m]]
    sample<-sample[order(sample$Week),]
    gfiles<-Gtype[as.vector(sample$File.name)]
    monkey<-names(monkeyLis)[m]
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
        merged<-merge(GT2,df[,c(1,2,3,4)], by="Genotype", all=T)
        merged[is.na(merged)]<-0
        merged$Week<-df$Week[1]
        merged$Tissue<-df$Tissue[1]
        
        Gdata<-rbind(Gdata,merged)
    }
    
    
    #select plasma only
    plasma<-Gdata[Gdata$Tissue=="Plasma"|Gdata$Tissue=="Stock",]
    gts<-as.character(unique(plasma$Genotype)) #31
    gts<-gts[!(gts %in% AAcombs)]
    gts<-gts[gts !="Singleton"]
    
    plasma$Genotype<-factor(plasma$Genotype, levels=c(AAcombs, paste(gts), "Singleton"))
    
    ggplot(data=plasma,aes(x=Week, y=Prop, fill=Genotype))+
        geom_area()+
        geom_vline(xintercept=tbweek, col="blue", size=0.3)+ylab("Proportion")+
        ggtitle(monkey)+
        scale_fill_discrete(type=colors,breaks=AAcombs,labels=c("WT",AAcombs[2:8]))
    ggsave(paste0("Output/AA/3GenotypeChange_",monkey,".pdf"),width = 9,height = 5)
    
    #which genotypes have the specific mutation?
    plasma$Genotype2<-plasma$Genotype
    plasma$Genotype<-as.character(plasma$Genotype)
    
    
    plots<-list()
    for (k in 1:3){
        plasma$mut<-sapply(plasma$Genotype, function(x) {
            if (x!="Singleton") {x=unlist(strsplit(x,"", fixed=T))
                                 ifelse (x[k]==mutAA[k], "Y","N")}
            else "N"})
      
       
        plots[[k]]<-ggplot(data=plasma,aes(x=Week, y=Prop, fill=Genotype2, pattern=mut))+
            geom_area_pattern(color = "gray80",size=0.05, 
                              pattern_fill = "gray60",
                              pattern_angle = 45,
                              pattern_density = 0.1,
                              pattern_spacing = 0.025,
                              pattern_key_scale_factor = 0.6)+
            theme_bw()+
            theme(legend.title = element_blank())+
            scale_pattern_manual(values = c(Y = "stripe", N = "none") ) +
            geom_vline(xintercept=tbweek, col="blue", size=0.3)+ylab("Proportion")+
            ggtitle(paste(monkey,aalabs[k]))+ 
            scale_fill_discrete(breaks=AAcombs,labels=c("WT",AAcombs[2:8]), type=colors)+
            guides(pattern = guide_legend(override.aes = list(fill = "white")),
                   fill = guide_legend(override.aes = list(pattern = "none")))
       
    }
    
    pdf(paste0("Output/AA/Figures/",monkey,".3genotypeFreq.pdf"), width = 8, height = 10)
    do.call(grid.arrange, c(plots, ncol=1))
    dev.off()
    
}

   
    
 