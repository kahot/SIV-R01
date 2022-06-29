library(reshape2)
library(colorspace)
library(DataCombine)
library(gridExtra)
source("Rscripts/baseRscript.R")
cols2<-qualitative_hcl(7, palette="Dark3")
#hcl_palettes("qualitative", n=8,plot = T)

#TB infection week
tbs<-read.csv("Data/TBinfection.week.csv",stringsAsFactors = F)

samples<-read.csv("Data/SamplesNoduplicates.csv",stringsAsFactors = F)
samples<-samples[samples$File.name!="Run6_18",]
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

#remove plasma-only or tissue-only monkeys
monkey3<-monkeys[!(monkeys=="A23918"|monkeys=="A24018"|monkeys=="A34019")]

### Sites of interest
diffpos<-c(475,494,500,515,518,578)
aapos<-ceiling(diffpos/3)
summary<-read.csv("Output/MF_PID/HighFreq/Summary_highFreq_sites.csv", row.names = 1)
summary2<-summary[summary$pos %in% diffpos,]
aaNames<-paste0(summary2$WTAA, aapos,summary2$aa1)

### Stock profile
stockaa<-read.csv(paste0("Output/AA/Run0_17AA.seq.csv"), stringsAsFactors = F, row.names = 1, colClasses = "character")
#select the 6 sites of interest
col.num <- which(colnames(stockaa) %in% paste0("pos.",aapos))
stockaa<-stockaa[,col.num]
stockaaFq<-apply(stockaa, 2, function(x) data.frame(table(x)))
stockaaProp<-lapply(stockaaFq, function(x) proportions(x[x$x!="X",2])) 

stockaaPerc<-data.frame(Week=0, Tissue="Plasma")
for (j in 1:length(stockaaFq)){
    df1<-stockaaFq[[j]]
    df1<-df1[df1$x!="X",]
    df1$Prop<-stockaaProp[[j]]*100
    if (summary2$aa1[j] %in% df1$x) stockaaPerc[1,aaNames[j]]<-df1$Prop[df1$x==summary2$aa1[j]]
    if (!(summary2$aa1[j] %in% df1$x)) stockaaPerc[1,aaNames[j]]<-0
}




## Calculate aa freq and create plots
AAs<-list.files("Output/AA/",pattern="AA.seq.csv$")
aaFreq6<-list()

#prep for plotting
#color assignment
mycolors<-c("gray","red", cols2[c(3,2,1,4,7,5,6)])

AA2<-data.fram
aaProp.List<-list()
for (m in 1:length(monkey3)){
    monkey<-monkey3[m]
    #print(monkey)
    sample<-monkeys2[[monkey]]
    sample = sample[order(sample[,'Week']),]
    
    #create freq of amino acids for each monkey
    aaPerc<-data.frame(Week=sample$Week, Tissue=sample$Tissue.type)
    #colnames(aaPerc)[3:8]<-colnames(aa)
    
    for (i in 1:nrow(sample)){
        fname<-sample$File.name[i]
        aa<-read.csv(paste0("Output/AA/",fname,"AA.seq.csv"), stringsAsFactors = F, row.names = 1, colClasses = "character")
        
        #select the 6 sites of interest
        col.num <- which(colnames(aa) %in% paste0("pos.",aapos))
        aa<-aa[,col.num]
        
        #calculate the freq of each site
        aaFq<-apply(aa, 2, function(x) data.frame(table(x)))
        #Convert it to freq.proportions(x[x$x!="X",2])777
        aaProp<-lapply(aaFq, function(x) proportions(x[x$x!="X",2])) 
        
        for (j in 1:length(aaFq)){
            df1<-aaFq[[j]]
            df1<-df1[df1$x!="X",]
            df1$Prop<-aaProp[[j]]*100
            if (summary2$aa1[j] %in% df1$x) aaPerc[i,aaNames[j]]<-df1$Prop[df1$x==summary2$aa1[j]]
            if (!(summary2$aa1[j] %in% df1$x)) aaPerc[i,aaNames[j]]<-0
        }
    }
    
    aaPerc<-InsertRow(aaPerc,stockaaPerc[1,] ,1)
    aa2<-aaPerc
    aa2$Run<-c(0, sample$Run)
    aa2$File.name<-c("Run0_16", sample$File.name)
    aa2$Monkey<-c("Stock", sample$Monkey)
    AA2<-rbind(AA2, aa2)
    
    aaPm<-melt(aaPerc, id.vars=c("Week","Tissue"))
    colnames(aaPm)[3:4]<-c("Position", "MF")
    
    tbweek<-tbs$tb[tbs$ids==monkey]
    
    types<-as.character(unique(aaPm$Tissue))
    types<-types[!(types=="Plasma"|types=="dLN"|types=="pLN"|types=="LN")]
    aaPm$Tissue<-factor(aaPm$Tissue, levels = c("Plasma","dLN","pLN","LN",types))
    n<-length(levels(aaPm$Tissue))
    mycolors2<-mycolors[1:n]
    names(mycolors2)<-levels(aaPm$Tissue)
    
    Pm<-aaPm[aaPm$Tissue=="Plasma",]
    
    #l<-nrow(sample[sample$Tissue!="Plasma"& sample$Week==max(sample$Week),])
    ggplot()+
        ylab("AA frequency")+xlab("Week")+
        facet_wrap(~ Position, ncol=3)+
        geom_line(data=Pm, aes(x=Week, y=MF),color="gray")+
        geom_point(data=aaPm, aes(x=Week, y=MF, color=Tissue),size=2.5, position=position_jitter(width = 0.6), alpha=0.7)+theme_bw()+
        ggtitle(paste0(monkey ))+
        geom_vline(xintercept=tbweek, col="blue")+
        scale_color_manual(values=mycolors2, drop=TRUE, )+
        scale_fill_discrete(drop=TRUE,limits = levels(aaPm$Tissue))
    
    #ggsave(paste0("Output/MF_PID/HighFreq/Tissue.AA.", monkey, ".pdf"), heigh=rown*1.8, width =7)
    
    #save it to a list
    aaPm$Monkey<-monkey
    aaProp.List[[m]]<-aaPm
    names(aaProp.List)[m]<-monkey

    
    print(monkey)
}


AA2<-AA2[!duplicated(AA2),]
#write.csv(AA2, "Output/MF_PID/HighFreq/AA.HighFreq.Metadata.csv")


plots<-list()
for (i in 1:length(aaNames)){
    position<-aaNames[i]
    P<-lapply(aaProp.List, function(x) x[x$Position==position,])
    Mm<-do.call(rbind, P)
    plots[[i]]<-ggplot()+
        ylab("AA frequency")+xlab("Week")+
        scale_color_manual(values=coltis)+
        facet_wrap(~ Monkey, ncol=8)+
        geom_line(data=Mm[Mm$Tissue=="Plasma",], aes(x=Week, y=MF),color="gray")+
        geom_point(data=Mm, aes(x=Week, y=MF, color=Tissue),size=2.5, position=position_jitter(width = 0.7), alpha=0.7)+theme_bw()+
        ggtitle(paste0(position))
    #ggsave(paste0("Output/MF_PID/HighFreq/Tissue.AA2.",position, ".pdf"), heigh=2.5, width =14)
    
}


pdf("Output/MF_PID/HighFreq/Tissue.AA2.All.pdf", width = 14, height = 14)
do.call(grid.arrange, c(plots, ncol=1))
dev.off()


### Run and read depth comparison ###
depth<-read.csv("Output/ReadDeapth_all.csv", row.names = 1)
AA2<-merge(AA2,depth[,c("File.name","Average")], by="File.name")
write.csv(AA2, "Output/MF_PID/HighFreq/AA.HighFreq.Metadata.csv")


### Compare run by run 
### Sites of interest
diffpos<-c(475,494,500,515,518,578)
aapos<-ceiling(diffpos/3)
summary<-read.csv("Output/MF_PID/HighFreq/Summary_highFreq_sites.csv", row.names = 1)
summary2<-summary[summary$pos %in% diffpos,]
aaNames<-paste0(summary2$WTAA, aapos,summary2$aa1)



AA2<-read.csv("Output/MF_PID/HighFreq/AA.HighFreq.Metadata.csv", stringsAsFactors = F, row.names = 1)

table(AA2$Run[AA2$R159G>20])
run.summary<-list()
for (i in 1:6){
    run.summary[[i]]<-table(AA2$Run[AA2[,paste0(aaNames[i])]>5])
    print(run.summary[[i]])
}

#read depth summary
reads<-data.frame(Mut=aaNames)
for (i in 1:6){
    reads$Mean[i]<-mean(AA2$Average[AA2[,paste0(aaNames[i])]>5])
    reads$Median[i]<-median(AA2$Average[AA2[,paste0(aaNames[i])]>5])
    ran<-range(AA2$Average[AA2[,paste0(aaNames[i])]>5])
    reads$Range_low[i]<-ran[1]
    reads$Range_high[i]<-ran[2]
    df<-mean(AA2[AA2[,paste0(aaNames[i])]>5,aaNames[i]])
    reads$AA.mean[i]<-mean(AA2[AA2[,paste0(aaNames[i])]>5,aaNames[i]])
}

aggregate(AA2[4:9], by=list(AA2$Run), mean)
run6<-AA2[AA2$Run==6,]
run5<-AA2[AA2$Run==5,]
run7<-AA2[AA2$Run==7,]

write.csv(run6,"Output/QC/run6_")

AA2<-merge(samples[,c("File.name","Tissue2")], AA2, by="File.name")


sum<-data.frame(aggregate(AA2[,4:9], by=list(AA2$Run), mean, na.rm=T))
sum

