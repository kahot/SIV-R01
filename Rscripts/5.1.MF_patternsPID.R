library(ggplot2)
library(gridExtra)
library(colorspace)
library(DataCombine)
#source("Rscripts/baseRscript2.R")
cols<-qualitative_hcl(6, palette="Dark3")

# read all Overview files:
#Use 'Filtered' overviews for those with over 100 average read 
depth<-read.csv("Output/ReadDeapth_all.csv", row.names = 1)
over100<-depth$File.name[depth$Average>100]

#(except for Run0_14 (500 reads) and some Run5 files with high read depths) 
filtered<-list.files("Output/OverviewF_PID/",pattern=".csv")

nonFiltered<-list.files("Output/Overview_PID/", pattern=".csv")

Overview<-list()
for (i in 1:length(filtered)){ 
        id<-substr(paste(filtered[i]),start=1,stop=7)
    
        if (id %in% over100) df<-read.csv(paste0("Output/OverviewF_PID/",filtered[i]),stringsAsFactors=F, row.names = 1)
        else {
            df<-read.csv(paste0("Output/Overview_PID/",nonFiltered[i]),stringsAsFactors=F, row.names = 1) 
            df[which(df$TotalReads<10), 7:16]<-NA #filter at 10
            }
        
        Overview[[i]]<-df
        names(Overview)[i]<-id
}


# Plot mutation freq across the genome for each week by monkeys

samples<-read.csv("Data/SamplesNoduplicates.csv",stringsAsFactors = F)
samples$Week<-as.integer(samples$Week)
samples$Tissue<-factor(samples$Tissue, levels=c("Stock","Plasma","LN","HLN","Lung","Unknown"))

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

#TB infection week
tbs<-read.csv("Data/TBinfection.week.csv",stringsAsFactors = F)

monkeys2<-monkeyList[tbs$ids]

#####   #Plot the stock virus ######
stock<-Overview[[which(names(Overview)=="Run0_17")]]
ave<-format(round(mean(stock$freq.mutations.ref, na.rm=T)*100,2),nsmall = 2)

Ps2<-ggplot(data=stock, aes(x=pos, y=freq.Ts.ref))+
    ylab("Mutation frequency")+xlab("")+ylim(0,1)+
    geom_point(size=0.7, color=cols[4])+theme_bw()+
    ggtitle("Stock")+theme(plot.title = element_text(size=12))+
    annotate(geom='text', x=600, y=0.96, label=paste0("MF: ", ave,"%"),color ='gray30', size=3,hjust=0)


#Plot mutation freq (Divergence)
Sum<-data.frame()
for (i in 1:length(monkeys2)){
    sample<-monkeys2[[i]]
    sample<-sample[order(sample$Week),]
    ovDF<-Overview[as.vector(sample$File.name)]
    monkey<-names(monkeys2)[i]
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
            annotate(geom='text', x=600, y=0.96, label=paste0("MF: ", ave,"%"),color ='gray30', size=3,hjust=0)+
            annotate(geom='text', x=600, y=0.82, label=paste0("Diversity: ", div,"%"),color ='gray30', size=3,hjust=0)
        
    }
    
    Sum<-rbind(Sum,summary)
    
    pdf(paste0("Output/MF_PID/MFsummary/MF.",monkey,".pdf"), width = 7, height = length(ovDF)*2)
    do.call(grid.arrange, c(Plot, ncol=1))
    dev.off()
}


#Plot minor variant freq (diversity)
for (i in 1:length(monkeys2)){
    sample<-monkeys2[[i]]
    sample<-sample[order(sample$Week),]
    ovDF<-Overview[as.vector(sample$File.name)]
    monkey<-names(monkeys2)[i]
    tbweek<-tbs$tb[tbs$ids==monkey]
    Plot<-list()
    Plot[[1]]<-Ps2
    for (j in 1:length(ovDF)){
        df<-ovDF[[j]]
        week<-sample$Week[j]
        if (week<=tbweek|is.na(tbweek)) col=cols[3]
        else {col=cols[1]}
        ts<-ifelse(sample$Tissue2[j]=="Plasma", "",sample$Tissue2[j])
        ave<-round(mean(df$freq.mutations.ref, na.rm=T)*100,digit= 2)
        div<-round(mean(df$freq.mutations, na.rm=T)*100,2)
        
        Plot[[j+1]]<-ggplot(data=df, aes(x=pos, y=freq.mutations*100))+
            ylab("Diversity")+xlab("")+ylim(0,50)+
            geom_point(size=0.7, color=col)+theme_bw()+
            ggtitle(paste0(monkey," Week ",week," ",ts))+
            theme(plot.title = element_text(size=11))+
            annotate(geom='text', x=600, y=47, label=paste0("MF: ", ave,"%"),color ='gray30', size=3,hjust=0)+
            annotate(geom='text', x=600, y=40, label=paste0("Diversity: ", div,"%"),color ='gray30', size=3,hjust=0)
        
    }
    pdf(paste0("Output/MF_PID/withRun7/MVF.",monkey,".pdf"), width = 7, height = length(ovDF)*2)
    do.call(grid.arrange, c(Plot, ncol=1))
    dev.off()
}


## Delete below after resequencing Run 2 samples ##

### LOOK at the read depth vs. trasnv freq (for low quality samples (Run2))#### 
## Run2 samples are messy in the middle (pos 485 - 550ish) ##
df<-Overview[["Run2_10"]]
mean(df$freq.Ts.ref, na.rm=T)
#0.01476899
mean(df$freq.transv.ref, na.rm=T)
#0.02486384

plot(df$freq.transv.ref)
plot(df$TotalReads)

pt<-list()
pt[[1]]<-ggplot(data=df, aes(x=pos, y=freq.mutations*100))+
    ylab("Diversity")+xlab("")+ylim(0,50)+
    geom_point(size=0.7, color=col)+theme_bw()
pt[[2]]<-ggplot(data=df, aes(x=pos, y=freq.Ts*100))+
    ylab("Diversity Ts")+xlab("")+ylim(0,50)+
    geom_point(size=0.7, color=col)+theme_bw()
pt[[3]]<-ggplot(data=df, aes(x=pos, y=freq.transv*100))+
    ylab("Diversity Tvs")+xlab("")+ylim(0,50)+
    geom_point(size=0.7, color=col)+theme_bw()

pt[[4]]<-ggplot(data=df, aes(x=pos, y=TotalReads))+
    ylab("Total Reads")+xlab("")+
    geom_point(size=0.7, color="blue")+theme_bw()

pdf(paste0("Output/MF_PID/withRun5.6/MVF.vs.Reads.A22117.week16.pdf"), width = 7, height = 8)
do.call(grid.arrange, c(pt, ncol=1))
dev.off()

pt<-list()
for (i in 1: length(Overview)){
    df<-Overview[[i]]
    filename<-names(Overview)[i]
    mk<-samples$Monkey[samples$File.name==filename]
    wk<-samples$Week[samples$File.name==filename]
    pt[[i]]<-ggplot(data=df, aes(x=pos, y=TotalReads))+
        ylab("Total Reads")+xlab("")+
        geom_point(size=0.7, color="blue")+theme_bw()+
        ggtitle(paste(filename, mk, "Week", wk))+
        theme(plot.title = element_text(size=11))
}
pdf(paste0("Output/MF_PID/withRun5.6/MVF.vs.Reads.allfiles.pdf"), width = 15, height = 97)
do.call(grid.arrange, c(pt, ncol=3))
dev.off()


####
### Assessing the quality of each samples: Compare # PID reads vs. diversity

reads<-read.csv("Output/ReadDeapth_all.csv", row.names = 1, stringsAsFactors = F)

Summary<-data.frame(File.name=names(Overview))
for(i in 1:length(Overview)){
    df<-Overview[[i]]
    Summary$Divergence[i]<-round(mean(df$freq.mutations.ref, na.rm=T)*100,digit= 2)
    Summary$Diversity[i]<-round(mean(df$freq.mutations, na.rm=T)*100,2)
}

Summary<-merge(reads,Summary, by="File.name")
plot(Summary$Diversity~Summary$Average, pch=16)
cor.test(Summary$Diversity, Summary$Average, method="kendall")
#p-value = 0.6702

Summary$File.name<-factor(Summary$File.name, levels=paste(Summary$File.name))
ggplot(Summary, aes(x=File.name, y=Diversity))+
    geom_point()
ggplot(Summary, aes(x=File.name, y=Average))+
    geom_point()+ylab('Read depth')+xlab('')

#No correlation between read depths and diversity

Summary<-merge(Summary, samples[,c("File.name","Monkey","Week","Coinfection")],by="File.name")
colnames(Summary)[4:5]<-c("Ave.reads","Max.read")
write.csv(Summary,"Output/Summary.Diversity_all.csv")





#### Compare the tissue and plasma diversity (%) at necropsy

monkeys
Nec<-data.frame()
Mid<-data.frame()
Div<-data.frame(Animal=monkeys)
for (i in 1:length(monkeys)){
    sum<-Summary[Summary$Monkey==monkeys[i],]
    
    #Midpoint tissue vs. plasma
    weeks<-sum$Week[sum$Tissue2=="Tissue"]
    if (length(weeks)==0) {
        Div$Mid_Plasma[i]<-NA
        Div$Mid_Tissue[i]<-NA
        Div$Nec_Plasma[i]<-NA
        Div$Nec_Tissue[i]<-NA
    }
    if (length(weeks)!=0) {
        if (length(unique(weeks))>1) {midpoint<-min(unique(weeks))
            mids<-(midpoint-2):(midpoint+2)
            sum1<-sum[sum$Week %in% mids,]
            Div$Mid_Plasma[i]<-sum1$Diversity[sum1$Tissue2=="Plasma"]
            Div$Mid_Tissue[i]<-sum1$Diversity[sum1$Tissue2=="Tissue"]
            Mid<-rbind(Mid, sum1)
        }
        else {
            Div$Mid_Plasma[i]<-NA
            Div$Mid_Tissue[i]<-NA
        }
        #Necropsy tissues vs. plasma
        sum2<-sum[sum$Week==max(sum$Week),]
        Nec<-rbind(Nec, sum2)
        div<-aggregate(sum2["Diversity"], by=list(sum2$Tissue2), mean)
        Div$Nec_Plasma[i]<-ifelse ("Plasma" %in% div$Group.1,  div$Diversity[div$Group.1=="Plasma"], NA)
        Div$Nec_Tissue[i]<-ifelse ("Tissue" %in% div$Group.1,  div$Diversity[div$Group.1=="Tissue"], NA)
    }
}

Div2<-Div[rowSums(is.na(Div))!=4,] 

Div2$Midpoint<-Div2$Mid_Plasma-Div2$Mid_Tissue
Div2$Necropsy<-Div2$Nec_Plasma-Div2$Nec_Tissue

co<-Summary[,c("Monkey","Coinfection")]
co<-co[!duplicated(co),]

Div2<-merge(Div2, co, by.x="Animal",by.y="Monkey")
Div2<-Div2[!is.na(Div2$Necropsy),]
write.csv(Div2, "Output/Diversity_tissue.vs.plasma.csv")

Div2m<-melt(Div2[,c("Animal","Midpoint","Necropsy","Coinfection")], id.vars=c("Animal","Coinfection"))
colnames(Div2m)[3:4]<-c("Timepoint","Diversity_difference") 

ggplot(Div2m, aes(x=Animal, y=Diversity_difference, fill=Timepoint, color=Coinfection))+
    geom_bar(stat="identity", position=position_dodge(width = 0.8))+
    ylab("Divrsity difference (Plasma-Tissues)")+xlab("")+
    theme_bw()+
    scale_fill_manual(values=cols[c(6,5)])+
    scale_color_manual(values=c("red","white"), 
            guide=guide_legend(override.aes=list(linetype=c(1,1),shape=c(NA,NA),color=c("red","white"), fill="gray90"))) +
    theme(panel.grid.major.x = element_blank())+
    geom_vline(xintercept = (1:(nrow(Div2)-1))+0.5, color="gray80", size=0.3)

ggsave("Output/MF_PID/Diversity_difference.Tissuesvs.Plasma.pdf", width = 7, height = 4)




############################
## Diversity Shift over time ##

monkeys2<-monkeyList[tbs$ids]
morder<-c("A22517","A22617","A22117","A22317","A22217","A21918","A23918","A34019","A34119","A34219")
plots<-list()
for (i in 1:length(morder)){
    sample<-monkeys2[[morder[i]]]
    sample<-sample[order(sample$Week),]
    ovDF<-Overview[as.vector(sample$File.name)]
    monkey<-morder[i]
    tbweek<-tbs$tb[tbs$ids==monkey]
    diversity<-sample[,c(2,3,5,7)]
    monkey<-gsub("A",'',monkey)
    for (j in 1:length(ovDF)){
        df<-ovDF[[j]]
        diversity$Diversity[j]<-mean(df$freq.mutations, na.rm=T)*100
        diversity$MF[j]<-mean(df$freq.mutations.ref, na.rm=T)*100
    }
    diversity<-InsertRow(diversity,c("Stock","Stock",0,"Stock",mean(stock$freq.mutations.ref, na.rm=T)*100,mean(stock$freq.mutations.ref, na.rm=T)*100 ),1)
    diversity$Week<-as.integer(diversity$Week)
    diversity$Diversity<-as.numeric(diversity$Diversity)
    diversity$MF<-as.numeric(diversity$MF)
    #diversity$Tissue<-factor(diversity$Tissue, levels=c("Stock","Plasma","LN","HLN","Unknown"))
    plots[[i]]<-ggplot()+
            ylab("Diversity")+xlab("")+
            geom_point(data=diversity, aes(x=Week, y=Diversity, color=Tissue), size=2)+theme_bw()+
            ggtitle(paste0(monkey))+ylim(0.4,2.2)+
            theme(plot.title = element_text(size=11))+
            scale_color_manual(values=cols[c(1,6,4,2,3)])+
            theme(panel.grid.major.x  = element_blank(),panel.grid.minor.x = element_blank())+
            geom_vline(xintercept=tbweek, col="deeppink")+
            scale_x_continuous(breaks=seq(0,32,2), limits = c(0,32))+
            geom_line(data=diversity[diversity$Tissue=="Plasma"|diversity$Tissue=="Stock",],aes(x=Week, y=Diversity), color=cols[6])
    
    
    
    
}
pdf("Output/MF_PID/withRun7/Diversity_overtime_all.pdf", width = 8.5, height = 11)
do.call(grid.arrange, c(plots, ncol=2))
dev.off()


### Overtime in one plot (Plasma only)

monkeys2<-monkeyList[tbs$ids]
#remove the monkey with Tissue only
monkeys<-names(monkeys2)
monkeys<-monkeys[monkeys!="A34019"]
monkeys3<-monkeys2[monkeys]

Div<-data.frame()
for (i in 1:length(monkeys)){
    sample<-monkeys2[[monkeys[i]]]
    sample<-sample[sample$Tissue=="Plasma",]
    sample<-sample[order(sample$Week),]
    ovDF<-Overview[as.vector(sample$File.name)]
    monkey<-monkeys[i]
    tbweek<-tbs$tb[tbs$ids==monkey]
    diversity<-sample[,c(3,5)]
    for (j in 1:length(ovDF)){
        df<-ovDF[[j]]
        diversity$Diversity[j]<-mean(df$freq.mutations, na.rm=T)*100
        diversity$MF[j]<-mean(df$freq.mutations.ref, na.rm=T)*100
    }
    
    diversity<-InsertRow(diversity,c(monkey,0,mean(stock$freq.mutations.ref, na.rm=T)*100,mean(stock$freq.mutations.ref, na.rm=T)*100 ),1)
    diversity$Monkey<-gsub("A",'',monkey)
    
    diversity$Week<-as.integer(diversity$Week)
    diversity$Diversity<-as.numeric(diversity$Diversity)
    diversity$MF<-as.numeric(diversity$MF)
    
    Div<-rbind(Div, diversity)
    
    
}

tbs2<-tbs[tbs$ids!="A34019",]

ggplot(data=Div, aes(x=Week, y=Diversity, color=Monkey))+
    ylab("% Diversity")+xlab("")+
    geom_point(size=2)+theme_bw()+
    geom_line()+
    #ylim(0.4,2.2)+
    #theme(plot.title = element_text(size=11))+
    #scale_color_manual(values=cols[c(1,6,4,2,3)])+
    theme(panel.grid.major.x  = element_blank(),panel.grid.minor.x = element_blank())+
geom_vline(data=tbs2, aes(xintercept=tbs2$tb,color=tbs2$ids),linetype = "dashed", size=0.5)
    #scale_x_continuous(breaks=seq(0,32,2), limits = c(0,32))+
    
ggsave("Output/MF_PID/Diversity_overtime_all.pdf", width = 10, height = 6)


plots[[i]]<-ggplot()+
    ylab("Diversity")+xlab("")+
    geom_point(data=diversity, aes(x=Week, y=Diversity, color=Tissue), size=2)+theme_bw()+
    ggtitle(paste0(monkey))+ylim(0.4,2.2)+
    theme(plot.title = element_text(size=11))+
    scale_color_manual(values=cols[c(1,6,4,2,3)])+
    theme(panel.grid.major.x  = element_blank(),panel.grid.minor.x = element_blank())+
    geom_vline(xintercept=tbweek, col="deeppink")+
    scale_x_continuous(breaks=seq(0,32,2), limits = c(0,32))+
    geom_line(data=diversity[diversity$Tissue=="Plasma"|diversity$Tissue=="Stock",],aes(x=Week, y=Diversity), color=cols[6])



