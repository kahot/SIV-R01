#Genotype analysis for the last 3-4 sites with P144Q
library(ggplot2)
library(gridExtra)
library(colorspace)
library(seqinr)
library(stringr)
library(reshape2)
library(ggpattern)
library(plotrix)
cols<-qualitative_hcl(8, palette="Dark3")

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
ID<-c("A21918","A22117","A22217","A22317","A22517","A22617","A23918","A34119","A34219")
monkeyLis<-monkeyList[ID]
monkeys<-names(monkeyLis)


########  Look at associations between P144Q and last 4 mutations
hsites<-read.csv("Output/MF_PID/HighFreq/Summary_highFreq_sites.csv", stringsAsFactors = F)
hsites<-hsites[hsites$pos!=485,]
sites<-hsites$pos
#AA sites (11)
AAsites<-ceiling(sites/3)
AAsites<-AAsites[c(7,9,10)]
#NT:297 305 326 328 340 373 431 515 536 539 548 
#AA: 99 102 109 110 114 125 144 172 179 180 183 

Ps<-paste0("pos.", AAsites)

hsites2<-hsites[c(7,9,10),]
aapairs<-apply(hsites2[,c("WTAA","aa1")], 1, function(x) c(x["WTAA"], x["aa1"]))

AApatterns<-expand.grid(aapairs[,1],aapairs[,2],aapairs[,3])
AAcombs<-apply(AApatterns, 1, function(x) paste0(x,collapse = ''))
AAcombs[AAcombs=="PNT"]<-"WT"

mutAA<-c("Q","S","N")
sites<-hsites2$pos
AAsites<-ceiling(sites/3)

#Wildtype AA ("DRVPANPRNTS")
ref<-read.csv("Output/Overview.ref.csv", row.names = 1,stringsAsFactors = F,colClasses = "character" )
ref<-ref[ref$pos %in% sites,]
WTaa<-toupper(paste0(ref$RefAA, collapse = ''))

#Create 3genotype frequency file
files<-list.files("Output/AA/Genotype11/",pattern="AA.11variants.csv$")
#select only plasma
plasma<-samples[samples$Tissue=="Plasma",]
files<-files[substring(files,1,7) %in% plasma$File.name]
genotypeCount<-data.frame(file=sub(".AA.11variants.csv",'',files))
for (i in 1:length(files)){
    df<-read.csv(paste0("Output/AA/Genotype11/",files[i]), stringsAsFactors = F, row.names = 1, colClasses = "character")
    fname<-substring(files[i],1,7)
    #first 3 sties only
    df<-df[,colnames(df) %in% Ps]
    #calculate the freq of each genotype
    df$Genotype<-apply(df[1:ncol(df)], 1, function(x) paste0(x,collapse = ''))
    genotypeCount$Total[i]<-length(unique(df$Genotype))
    
    #remove the sequences with unknown aa  ('X')  
    df2<-df[!grepl("X", df$Genotype, fixed=T),]
    genotypeCount$no_Xs[i]<-length(unique(df2$Genotype))
    
    #Genotype freq. all
    gt<-data.frame(table(df2$Genotype))
    gt<-gt[order(gt$Freq, decreasing = T),]
    gt$Var1<-as.character(gt$Var1)
    
    #Replace the wildtype (Reference) with WT
    gt$Var1[which(gt$Var1==WTaa)]<-"WT"
    genotypeCount$TotalReads[i]<-nrow(df)
    genotypeCount$WTcount[i]<-ifelse(length(gt$Freq[gt$Var1=="WT"])==0, 0,gt$Freq[gt$Var1=="WT"])
    
    write.csv(gt,paste0("Output/AA/Genotype.freq7.9.10//",fname,".AA.freq.csv"))
    print(fname)
}

#write.csv(genotypeCount,"Output/AA/AA.Genotype3_count_summary.csv")


files<-list.files("Output/AA/Genotype.freq7.9.10/", pattern="freq.csv")
Gtype<-list()
for (i in 1: length(files)){
  df<-read.csv(paste0("Output/AA/Genotype.freq7.9.10/",files[i]), stringsAsFactors = F, row.names = 1)
  fname<-substring(files[i],1,7)
  Gtype[[i]]<-df
  names(Gtype)[i]<-fname
}


stock<-read.csv("Output/AA/Genotype.freq7.9.10/Run0_17.AA.freq.csv", stringsAsFactors = F, row.names = 1)
stock$Freq<-as.integer(stock$Freq)
stock$Prop<-stock$Freq/sum(stock$Freq)
stock$Week<-0
stock$Tissue<-"Stock"
write.csv(stock, "Output/AA/Analysis/pos144.179.180.wt_proportion.csv")

#####

cols<-qualitative_hcl(8, palette="Dark3")
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
color<-color[color!="white"]

colors<-c("gray", cols, color)
aalabs<-c("P144Q","N179S","T180N")

for (m in 1:length(monkeyLis)){
    #Select the monkey
    sample<-monkeyLis[[m]]
    #plasma only
    sample<-sample[sample$Tissue=="Plasma",]
    sample<-sample[order(sample$Week),]
    gfiles<-Gtype[as.vector(sample$File.name)]
    monkey<-names(monkeyLis)[m]
    tbweek<-tbs$tb[tbs$ids==monkey]
    gtypes<-unique(stock$Var1)
    
    #Createa an aggregated data frame of gtype freq. for each monkey
    DF<-list()
    DF[[1]]<-stock
    for (j in 1:length(gfiles)){
        df2<-gfiles[[j]]
        df2$Freq<-as.integer(df2$Freq)
        df2$Prop<-df2$Freq/sum(df2$Freq)
        df2$Week<-sample$Week[j]
        df2$Tissue<-paste0(sample$Tissue[j])
        
        gtypes<-c(gtypes,unique(df2$Var1))
        DF[[j+1]]<-df2
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
        Gdata<-rbind(Gdata,merged)
    }
    
    
    gts<-as.character(unique(Gdata$Genotype)) #31
    gts<-gts[!(gts %in% AAcombs)]
    
    Gdata$Genotype<-factor(Gdata$Genotype, levels=c(AAcombs, paste(gts)))
    
    #ggplot(data=Gdata,aes(x=Week, y=Prop, fill=Genotype))+
    #    geom_area()+
    #    geom_vline(xintercept=tbweek, col="blue", size=0.3)+ylab("Proportion")+
    #    ggtitle(monkey)+
    #    scale_fill_discrete(type=colors,breaks=c(AAcombs),labels=c(AAcombs[1:8]))
    #ggsave(paste0("Output/AA/Analysis/3GenotypeChange_",monkey,".pdf"),width = 9,height = 5)
    
    #which genotypes have the specific mutation?
    Gdata$Genotype2<-Gdata$Genotype
    Gdata$Genotype<-as.character(Gdata$Genotype)
    
    #Replace WT with DRP
    Gdata$Genotype[Gdata$Genotype=="WT"]<-"PNT"
    plots<-list()
    for (k in 1:3){
        #Track each fate of mutation
        Gdata$mut<-sapply(Gdata$Genotype, function(x) {
            x=unlist(strsplit(x,"", fixed=T))
            ifelse (x[k]==mutAA[k], "Y","N")})
        
        plots[[k]]<-ggplot(data=Gdata,aes(x=Week, y=Prop, fill=Genotype2, pattern=mut))+
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
            scale_fill_discrete(breaks=AAcombs,labels=AAcombs, type=colors,drop = FALSE)+
            guides(pattern = guide_legend(override.aes = list(fill = "white")),
                   fill = guide_legend(override.aes = list(pattern = "none")))
        
    }
    
    pdf(paste0("Output/AA/Analysis/",monkey,".144.179.180genotypeFreq.pdf"), width = 8, height = 10)
    do.call(grid.arrange, c(plots, ncol=1))
    dev.off()
    
}






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
    
    ids<-c("0.Stock", paste0(sample$Week,".",sample$Tissue2))
    assoc<-data.frame(id=ids)
    assoc$week<-c(0,sample$Week)
    assoc$tissue<-c("Stock", sample$Tissue2)
    assoc2<-assoc
    for (w in 1:length(ids)){
        df<-Gdata[Gdata$id==ids[w],]
        assoc$D99E[w]<-sum(df$Freq[df$d99e=="Y"&df$p144q=="N"])/sum(df$Freq[df$p144q=="N"])*100
        assoc$R102K[w]<-sum(df$Freq[df$r102k=="Y"&df$p144q=="N"])/sum(df$Freq[df$p144q=="N"])*100
    }
    assoc$D99E[is.nan(assoc$D99E)]<-0
    assoc$R102K[is.nan(assoc$R102K)]<-0
    
    
    assocm<-melt(assoc, id.vars = c("id","week","tissue"))
    assocm$id<-factor(assocm$id, levels=paste(assoc$id))
    colnames(assocm)[4]<-"Mutation"
    if (max(assocm$value)<50) ymax=50
    else ymax=max(assocm$value)
    
    #Plot % of D99E or R102K that occurs without P144Q (out of all nonP144Q muatations)  
    plots[[m]]<-ggplot(data=assocm, aes(x=id, y=value, fill=Mutation))+
        geom_bar(stat="identity", position = position_dodge(width = .8),alpha=0.8)+
        theme(axis.text.x = element_text(angle=90, hjust=1))+
        ggtitle(monkey)+xlab('')+
        scale_fill_manual(values = cols[c(1,7)])+
        ylab("% mut w/o P144Q/nonP144Q")+ylim(0,ymax)
    
   
    for (w in 1:length(ids)){
        df<-Gdata[Gdata$id==ids[w],]
        assoc2$D99E[w]<-sum(df$Freq[df$d99e=="N"&df$p144q=="Y"])/sum(df$Freq[df$p144q=="Y"])*100
        assoc2$R102K[w]<-sum(df$Freq[df$r102k=="N"&df$p144q=="Y"])/sum(df$Freq[df$p144q=="Y"])*100
    }
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
        scale_fill_manual(values = cols[c(4,8)], labels=c("no D99E","no R102K"))
}


pdf(paste0("Output/AA/Analysis/NoP144Q_withOther2Mutations.pdf"), width =15, height =10)
do.call(grid.arrange, c(plots, ncol=3))
dev.off()

#Find the opposite -P144Q present but not others
pdf(paste0("Output/AA/Analysis/P144Qpresent_noOthermutation2.pdf"), width = 15, height =10)
do.call(grid.arrange, c(plots2, ncol=3))
dev.off()



#######

