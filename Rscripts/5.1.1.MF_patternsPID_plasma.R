library(ggplot2)
library(gridExtra)
library(colorspace)
library(DataCombine)
cols<-qualitative_hcl(6, palette="Dark3")

# Read all Overview files:
filtered<-list.files("Output/OverviewF_PID_plasma/",pattern=".csv")
Overview1<-list()
for (i in 1:length(filtered)){ 
        id<-substr(paste(filtered[i]),start=1,stop=7)
        df<-read.csv(paste0("Output/OverviewF_PID_plasma/",filtered[i]),stringsAsFactors=F, row.names = 1)
        Overview1[[i]]<-df
        names(Overview1)[i]<-id
}


# Plot mutation freq across the genome for each week by monkeys

samples<-read.csv("Data/SamplesNoduplicates.csv",stringsAsFactors = F)
samples$Week<-as.integer(samples$Week)
samples$Tissue<-factor(samples$Tissue, levels=c("Stock","Plasma","LN","HLN","Lung","Unknown"))

#Animals with plasma samples
animals<-c("A21918","A22117","A22217","A22317","A22517","A22617","A23918","A34119","A34219")
samples<-samples[samples$Monkey %in% animals,] #95 files

list.animal<-split(samples, samples$Monkey)
#Remove tissue only A34019
list.animal<-list.animal[names(list.animal)!="A34019"]
monkeyList<-list()
k=1
for (i in 1:length(list.animal)){
    if (nrow(list.animal[[i]])>1){
        monkeyList[[k]]<-list.animal[[i]]
        names(monkeyList)[k]<-names(list.animal)[i]
        k=k+1
    }
}

#TB infection week
tbs<-read.csv("Data/TBinfection.week.csv",stringsAsFactors = F)
tbs<-tbs[tbs$ids %in% animals,]

#####   #Plot the stock virus ######
stock<-Overview1[[which(names(Overview1)=="Run0_17")]]
ave<-format(round(mean(stock$freq.mutations.ref, na.rm=T)*100,2),nsmall = 2)

Ps2<-ggplot(data=stock, aes(x=pos, y=freq.mutations.ref))+
    ylab("Mutation frequency")+xlab("")+ylim(0,1)+
    geom_point(size=0.7, color=cols[4])+theme_bw()+
    ggtitle("Stock")+theme(plot.title = element_text(size=12))+
    geom_rect(data=stock, inherit.aes=FALSE,
              aes(xmin=458,xmax=506,ymin=-Inf,ymax=Inf),
              fill="gray90",size=1,alpha=0.2)+
    annotate(geom='text', x=600, y=0.96, label=paste0("MF: ", ave,"%"),color ='gray30', size=3,hjust=0)



#Plot mutation freq (Divergence) #plasma only
Sum<-data.frame()
for (i in 1:length(monkeyList)){
    sample<-monkeyList[[i]]
    sample<-sample[sample$Tissue2=="Plasma",]
    sample<-sample[order(sample$Week),]
    
    ovDF<-Overview1[as.vector(sample$File.name)]
    monkey<-names(monkeyList)[i]
    tbweek<-tbs$tb[tbs$ids==monkey]
    
    Plot<-list()
    Plot[[1]]<-Ps2
    summary<-data.frame(File.name=sample$File.name)
    for (j in 1:length(ovDF)){
        df<-ovDF[[j]]
        week<-sample$Week[j]
        
        if (week<=tbweek|is.na(tbweek)) col=cols[3]
        else {col=cols[1]}
        ts<-ifelse(sample$Tissue2[j]=="Plasma", "",sample$Tissue2[j])
        ave<-round(mean(df$freq.mutations.ref, na.rm=T)*100,digit= 2)
        div<-round(mean(df$freq.mutations, na.rm=T)*100,2)
        
        summary$Ave.diversity[j]<-div
        summary$Ave.divergence[j]<-ave
        
        Plot[[j+1]]<-ggplot(data=df, aes(x=pos, y=freq.mutations.ref))+
            ylab("Mutation frequency")+xlab("")+ylim(0,1)+
            geom_point(size=0.7, color=col)+theme_bw()+
            ggtitle(paste0(monkey," Week ",week," ",ts))+
            theme(plot.title = element_text(size=11))+
            geom_rect(data=stock, inherit.aes=FALSE,
                      aes(xmin=458,xmax=506,ymin=-Inf,ymax=Inf),
                      fill="gray90",size=1,alpha=0.2)+
            annotate(geom='text', x=600, y=0.96, label=paste0("MF: ", ave,"%"),color ='gray30', size=3,hjust=0)+
            annotate(geom='text', x=600, y=0.82, label=paste0("Diversity: ", div,"%"),color ='gray30', size=3,hjust=0)
        
    }
    
    Sum<-rbind(Sum,summary)
    
    pdf(paste0("Output/MF_PID/MFsummary_trimmed/MF.",monkey,".pdf"), width = 7, height = length(ovDF)*2)
    do.call(grid.arrange, c(Plot, ncol=1))
    dev.off()
}


####
### Assessing the quality of each samples: Compare # PID reads vs. diversity
#remove the center portion for diversity comparison

reads<-read.csv("Output/ReadDeapth_all.csv", row.names = 1, stringsAsFactors = F)

Summary<-data.frame(File.name=names(Overview1))
for(i in 1:length(Overview1)){
    df<-Overview1[[i]]
    Summary$Divergence[i]<-round(mean(df$freq.mutations.ref, na.rm=T)*100,digit= 4)
    Summary$Diversity[i]<-round(mean(df$freq.mutations, na.rm=T)*100,4)
}

Summary<-merge(reads,Summary, by="File.name")
plot(Summary$Diversity~Summary$Average, pch=16)
cor.test(Summary$Diversity, Summary$Average, method="spearman")
#p-value = 0.6788
#rho -0.05540004  

Summary$File.name<-factor(Summary$File.name, levels=paste(Summary$File.name))
ggplot(Summary, aes(x=File.name, y=Diversity))+
    geom_point()+theme(axis.text.x = element_text(angle=90, size=5))
#ggsave("Output/MF_PID/AveDiversity_allsamples.pdf", width = 8, height = 5)
ggplot(Summary, aes(x=File.name, y=Average))+
    geom_point()+ylab('Read depth')+xlab('')+
    theme(axis.text.x = element_text(angle=90, size=5))
#ggsave("Output/MF_PID/ReadDepth_allsamples.pdf", width = 8, height = 5)
#No correlation between read depths and diversity


Summary<-merge(Summary, samples[,c("File.name","Monkey","Week","Coinfection")],by="File.name")
colnames(Summary)[4:5]<-c("Ave.reads","Max.read")
write.csv(Summary,"Output/Summary.Diversity_all_gapRemoved.csv")


############################
## Diversity Shift over time ##

morder<-c("A22517","A22617","A22117","A22317","A22217","A21918","A23918","A34119","A34219")
plots<-list()
for (i in 1:length(morder)){
    sample<-monkeyList[[morder[i]]]
    sample<-sample[sample$Tissue2=="Plasma",]
    sample<-sample[order(sample$Week),]
    ovDF<-Overview1[as.vector(sample$File.name)]
    monkey<-morder[i]
    tbweek<-tbs$tb[tbs$ids==monkey]
    diversity<-sample[,c(2,3,4)]
    monkey<-gsub("A",'',monkey)
    for (j in 1:length(ovDF)){
        df<-ovDF[[j]]
        diversity$Diversity[j]<-mean(df$freq.mutations, na.rm=T)*100
        diversity$MF[j]<-mean(df$freq.mutations.ref, na.rm=T)*100
    }
    diversity<-InsertRow(diversity,c("Stock","Stock",0,mean(stock$freq.mutations.ref, na.rm=T)*100,mean(stock$freq.mutations.ref, na.rm=T)*100 ),1)
    diversity$Week<-as.integer(diversity$Week)
    diversity$Diversity<-as.numeric(diversity$Diversity)
    diversity$MF<-as.numeric(diversity$MF)
    #diversity$Tissue<-factor(diversity$Tissue, levels=c("Stock","Plasma","LN","HLN","Unknown"))
    plots[[i]]<-
        ggplot()+
            ylab("Diversity")+xlab("")+
            geom_point(data=diversity, aes(x=Week, y=Diversity), color=4,size=2)+theme_bw()+
            ggtitle(paste0(monkey))+ylim(0.4,2.2)+
            theme(plot.title = element_text(size=11))+
            theme(panel.grid.major.x  = element_blank(),panel.grid.minor.x = element_blank())+
            geom_vline(xintercept=tbweek, col="deeppink")+
            scale_x_continuous(breaks=seq(0,32,2), limits = c(0,32))+
            geom_line(data=diversity,aes(x=Week, y=Diversity), color=4)
}
pdf("Output/MF_PID/Diversity_overtime_plasma_only.pdf", width = 8.5, height = 11)
do.call(grid.arrange, c(plots, ncol=2))
dev.off()

#median
plots<-list()
for (i in 1:length(morder)){
    sample<-monkeyList[[morder[i]]]
    sample<-sample[sample$Tissue2=="Plasma",]
    sample<-sample[order(sample$Week),]
    ovDF<-Overview1[as.vector(sample$File.name)]
    monkey<-morder[i]
    tbweek<-tbs$tb[tbs$ids==monkey]
    diversity<-sample[,c(2,3,4)]
    monkey<-gsub("A",'',monkey)
    for (j in 1:length(ovDF)){
        df<-ovDF[[j]]
        diversity$Diversity[j]<-median(df$freq.mutations, na.rm=T)*100
        diversity$MF[j]<-median(df$freq.mutations.ref, na.rm=T)*100
    }
    diversity<-InsertRow(diversity,c("Stock","Stock",0,median(stock$freq.mutations.ref, na.rm=T)*100,median(stock$freq.mutations.ref, na.rm=T)*100 ),1)
    diversity$Week<-as.integer(diversity$Week)
    diversity$Diversity<-as.numeric(diversity$Diversity)
    diversity$MF<-as.numeric(diversity$MF)
    #diversity$Tissue<-factor(diversity$Tissue, levels=c("Stock","Plasma","LN","HLN","Unknown"))
    plots[[i]]<-ggplot()+
        ylab("Diversity")+xlab("")+
        geom_point(data=diversity, aes(x=Week, y=Diversity), color=4,size=2)+theme_bw()+
        ggtitle(paste0(monkey))+ylim(0,0.5)+
        theme(plot.title = element_text(size=11))+
        theme(panel.grid.major.x  = element_blank(),panel.grid.minor.x = element_blank())+
        geom_vline(xintercept=tbweek, col="deeppink")+
        scale_x_continuous(breaks=seq(0,32,2), limits = c(0,32))+
        geom_line(data=diversity,aes(x=Week, y=Diversity), color=4)
}
    
    
pdf("Output/MF_PID/Diversity_overtime_plasma_median.pdf", width = 8.5, height = 11)
do.call(grid.arrange, c(plots, ncol=2))
dev.off()





### Overtime in one plot (Plasma only)
morder<-c("A22517","A22617","A22117","A22317","A22217","A21918","A23918","A34119","A34219")
monkeys3<-monkeyList[morder]

Div<-data.frame()
for (i in 1:length(morder)){
    sample<-monkeys3[[morder[i]]]
    sample<-sample[sample$Tissue2=="Plasma",]
    sample<-sample[order(sample$Week),]
    ovDF<-Overview1[as.vector(sample$File.name)]
    monkey<-names(monkeys3)[i]
    tbweek<-tbs$tb[tbs$ids==monkey]
    diversity<-sample[,c(4,5)]
    for (j in 1:length(ovDF)){
        df<-ovDF[[j]]
        diversity$Diversity[j]<-mean(df$freq.mutations, na.rm=T)*100
        diversity$MF[j]<-mean(df$freq.mutations.ref, na.rm=T)*100
    }
    diversity$Tissue<-as.character(diversity$Tissue)
    
    diversity<-InsertRow(diversity,c(0,"Stock",mean(stock$freq.mutations.ref, na.rm=T)*100,mean(stock$freq.mutations.ref, na.rm=T)*100 ),1)
    diversity$Monkey<-gsub("A",'',monkey)
    
    diversity$Week<-as.integer(diversity$Week)
    diversity$Diversity<-as.numeric(diversity$Diversity)
    diversity$MF<-as.numeric(diversity$MF)
    
    Div<-rbind(Div, diversity)
}

tbs$ids<-gsub("A",'',tbs$ids)


Div$coinfection<-"Y"
Div$coinfection[Div$Monkey %in% c("34219", "34119")]<-"N"
Div$coinfection<-factor(Div$coinfection, levels = c("Y","N"))


ggplot(data=Div, aes(x=Week, y=Diversity, color=Monkey))+
    ylab("% Diversity")+xlab("")+
    geom_point(size=2)+theme_bw()+
    geom_path(data=Div, aes(x=Week, y=Diversity,color=Monkey, linetype=coinfection))+
    ylim(0.3,2)+
    theme(panel.grid.major.x  = element_blank(),panel.grid.minor.x = element_blank())+
    geom_vline(data=tbs, aes(xintercept=tbs$tb,color=tbs$ids),linetype = "dashed", size=0.3)
    #scale_x_continuous(breaks=seq(0,32,2), limits = c(0,32))+
    
ggsave("Output/MF_PID/Diversity_overtime_all.pdf", width = 8, height = 4)



