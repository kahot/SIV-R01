#Script to analyse the frequency data and associate with features
library(dplyr)
library(ggpubr)
library(reshape2)
library(colorspace)
source("Rscripts/baseRscript.R")
cols2<-qualitative_hcl(6, palette="Dark3")

#get the file name
SIVFiles_SeqDataP<-list.files("Output/SeqData_PID/",pattern="SeqData")


#where do the indels occur?
Ins<-data.frame(pos=190:681)
Del<-data.frame(pos=190:681)
Deletion<-list()
Indels<-list()
for (i in 1:length(SIVFiles_SeqDataP)){   
        id<-substr(paste(SIVFiles_SeqDataP[i]),start=9,stop=15)
        print(id)
        DF<-read.csv(paste0("Output/SeqData_PID/",SIVFiles_SeqDataP[i]),row.names = 1, stringsAsFactors=FALSE)
        #DF<-DF[!is.na(DF$a),]
        DF$del.percent<-DF$deletion/DF$TotalReads
        DF$ins.percent<-DF$insertion/DF$TotalReads
        
        df2<-DF[,c("pos","del.percent","ins.percent")]
        colnames(df2)[2:3]<-c("deletion","insertion")
        dfm<-melt(df2, id.vars="pos")
        dfm[dfm==0]<-NA
        ggplot(dfm, aes(x=pos, y=value, color=variable))+
            geom_point(size=.5)+
            scale_color_manual(values=cols2[c(2,5)])+
            facet_grid(rows = vars(variable))+
            theme_bw()+ylab('Indel observed')+
            xlab('Genome position')+
            guides(color = guide_legend(title = NULL))+
            ggtitle(id)
        ggsave(paste0("Output/Indels/", id, ".indel.observed.pdf"), height = 3, width = 6)
        
        Indels[[i]]<-df2
        names(Indels)[[i]]<-id
        
        Ins[,id]<-df2[,"insertion"]
        Del[,id]<-df2[,"deletion"]
}     

n<-length(SIVFiles_SeqDataP)
Ins$mean<-rowMeans(Ins[,2:(n+1)], na.rm=T)
Ins$count<-rowSums(Ins[,2:(n+1)]>0)
Del$mean<-rowMeans(Del[,2:(n+1)], na.rm=T)
Del$count<-rowSums(Del[,2:(n+1)]>0)

write.csv(Ins, "Output/Indels/Insertion.freq.csv")
write.csv(Del, "Output/Indels/Deletion.freq.csv")

###
Ins<-read.csv("Output/Indels/Insertion.freq.csv", stringsAsFactors = F, row.names = 1)
Del<-read.csv("Output/Indels/Deletion.freq.csv", stringsAsFactors = F, row.names = 1)

Ins$AApos<-ceiling(Ins$pos/3)

Insertions<-Ins[,c("pos","AApos","count","mean")]

insDF<-Ins[,c("pos","mean","count","AApos")]
colnames(insDF)[2:3]<-c("freq","No.of.samples")
#insDFm<-melt(insDF, id.vars="pos")

n<-length(SIVFiles_SeqDataP)

insDF$freq[insDF$freq==0]<-NA
insDF$No.of.samples[insDF$No.of.samples==0]<-NA
insDF$Prop.of.samples<-insDF$No.of.samples/n
max(insDF$Prop.of.samples, na.rm = T)
coeff<-100
ggplot(insDF, aes(x=pos))+
    geom_point(aes(y=Prop.of.samples*100),size=0.5, color=cols2[5])+
    geom_point(aes(y=freq*100), size=0.5, color=cols2[1])+
    scale_y_continuous(
        name="% sample with insertion", 
        sec.axis=sec_axis(~./coeff, name="Insertion freq"))+
    theme_bw()+
    xlab('Genome position')+
    annotate(geom='pointrange', x=600, y=100, ymin=50, ymax=50, color=cols2[5], size=.3)+
    annotate(geom='text', x=610, y=100, label="% sample", color=1, size=3, hjust=0)+
    annotate(geom='pointrange', x=600, y=94, ymin=50, ymax=50, color=cols2[1], size=.3)+
    annotate(geom='text', x=610, y=94, label="Average freq.", color=1, size=3, hjust=0)
ggsave("Output/Indels/Insertion_summary.pdf", width = 7, height = 3.5)  


ggplot(insDF, aes(x=AApos))+
    #geom_bar(aes(y=Prop.of.samples*100),stat="identity",color=paste0(cols2[5],"66"), fill=paste0(cols2[5],"66"))+
    geom_point(aes(y=freq*100), size=0.5, color=cols2[1])+
    theme_bw()+ylab('Insertion frequency')
    xlab('Amino acid position')
    annotate(geom='pointrange', x=600, y=100, ymin=50, ymax=50, color=cols2[5], size=.3)+
    annotate(geom='text', x=610, y=100, label="% sample", color=1, size=3, hjust=0)+
    annotate(geom='pointrange', x=600, y=92, ymin=46, ymax=46, color=cols2[1], size=.3)+
    annotate(geom='text', x=610, y=92, label="Average freq.", color=1, size=3, hjust=0)
ggsave("Output/Indels/Insertion_summary2.pdf", width = 7, height = 3.5)  

#AApos
ggplot(insDF, aes(x=pos))+
    geom_point(aes(y=Prop.of.samples*100),size=0.5, color=cols2[5])+
    geom_point(aes(y=freq*100), size=0.5, color=cols2[1])+
    scale_y_continuous(
        name="% sample with insertion", 
        sec.axis=sec_axis(~./coeff, name="Insertion freq"))+
    theme_bw()+
    xlab('Genome position')+
    annotate(geom='pointrange', x=600, y=100, ymin=50, ymax=50, color=cols2[5], size=.3)+
    annotate(geom='text', x=610, y=100, label="% sample", color=1, size=3, hjust=0)+
    annotate(geom='pointrange', x=600, y=94, ymin=50, ymax=50, color=cols2[1], size=.3)+
    annotate(geom='text', x=610, y=94, label="Average freq.", color=1, size=3, hjust=0)
ggsave("Output/Indels/Insertion_summary.pdf", width = 7, height = 3.5)  




####

delDF<-Del[,c("pos","mean","count")]
colnames(delDF)[2:3]<-c("freq","No.of.samples")

delDF$freq[delDF$freq==0]<-NA
delDF$No.of.samples[delDF$No.of.samples==0]<-NA
delDF$Prop.of.samples<-delDF$No.of.samples/n

coeff<-100


ggplot(delDF, aes(x=pos))+
    geom_bar(aes(y=Prop.of.samples*100),stat="identity",color=paste0(cols2[5],"66"), fill=paste0(cols2[5],"66"))+
    geom_point(aes(y=freq*100), size=0.5, color=cols2[1])+
    scale_y_continuous(
        name="% sample with deletion", 
        sec.axis=sec_axis(~./coeff, name="Deletion freq"))+
    theme_bw()+
    xlab('Genome position')+
    annotate(geom='pointrange', x=600, y=100, ymin=50, ymax=50, color=cols2[5], size=.3)+
    annotate(geom='text', x=610, y=100, label="% sample", color=1, size=3, hjust=0)+
    annotate(geom='pointrange', x=600, y=92, ymin=46, ymax=46, color=cols2[1], size=.3)+
    annotate(geom='text', x=610, y=92, label="Average freq.", color=1, size=3, hjust=0)
ggsave("Output/Indels/Deletion_summary2.pdf", width = 7, height = 3.5)  


####################
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




highIns<-Ins[Ins$mean>=0.005,]
highs<-highIns$pos

for (i in 1:length(monkeys2)){
    sample<-monkeys2[[i]]
    sample = sample[order(sample[,'Week'],-sample[,'Run']),]
    #sample = sample[!duplicated(sample$Week),]
    
    insdf<-Ins[,c('pos', sample$File.name)]

    insdf<-insdf[insdf$pos %in% highs,]
    
    insd<-data.frame(t(insdf))
    colnames(insd)<-insd[1,]
    insd<-insd[-1,]
    insd$File.name<-rownames(insd)
    insd<-merge(insd, sample[,c("File.name","Week")], by="File.name")
    insd<-insd[order(insd$Week),]
    insd<-insd[,-1]
    insdM<-melt(insd,id.vars = 'Week')
    colnames(insdM)[2:3]<-c("Position","Freq")
    monkey<-names(monkeys2)[i]
    
    tbweek<-tbs$tb[tbs$ids==monkey]
    
    ggplot(data=insdM, aes(x=Week, y=Freq))+
        ylab("Insertion frequency")+xlab("Week")+
        facet_wrap(~ Position, nrow=4, ncol=5)+
        geom_point(size=2, color=cols2[1])+theme_bw()+ggtitle(paste(monkey))+
        geom_vline(xintercept=tbweek, col="blue")
    ggsave(paste0("Output/Indels/Timeseries/Inser.freq_timeseries_",monkey, ".pdf"), width=7, height=4)
}


##deletion that appears in most files -> don't use
manyIns<-Ins[Ins$count>=105&!is.na(Ins$count),]
many<-manyIns$pos

for (i in 1:length(monkeys2)){
    sample<-monkeys2[[i]]
    sample = sample[order(sample[,'Week'],-sample[,'Run']),]
    sample = sample[!duplicated(sample$Week),]
    
    insdf<-Ins[,c('pos', sample$File.name)]
    
    insdf<-insdf[insdf$pos %in% many,]
    
    insd<-data.frame(t(insdf))
    colnames(insd)<-insd[1,]
    insd<-insd[-1,]
    insd$File.name<-rownames(insd)
    insd<-merge(insd, sample[,c("File.name","Week")], by="File.name")
    insd<-insd[order(insd$Week),]
    insd<-insd[,-1]
    insdM<-melt(insd,id.vars = 'Week')
    colnames(insdM)[2:3]<-c("Position","Freq")
    monkey<-names(monkeys2)[i]
    
    tbweek<-tbs$tb[tbs$ids==monkey]
    ggplot(data=insdM, aes(x=Week, y=Freq))+
        ylab("Insertion frequency")+xlab("Week")+
        facet_wrap(~ Position, nrow=1, ncol=6)+
        geom_point(size=2, color=cols2[1])+theme_bw()+ggtitle(paste(monkey))+
        geom_vline(xintercept=tbweek, col="blue")
    ggsave(paste0("Output/Indels/Timeseries/Common_Inser.freq_timeseries_",monkey, ".pdf"), width=7, height=2)
}



## deletion
highDel<-Del[Del$mean>=0.005,]
highd<-highDel$pos

for (i in 1:length(monkeys2)){
    sample<-monkeys2[[i]]
    sample = sample[order(sample[,'Week'],-sample[,'Run']),]
    sample = sample[!duplicated(sample$Week),]
    
    deldf<-Del[,c('pos', sample$File.name)]
    
    deldf<-deldf[deldf$pos %in% highd,]
    
    deld<-data.frame(t(deldf))
    colnames(deld)<-deld[1,]
    deld<-deld[-1,]
    deld$File.name<-rownames(deld)
    deld<-merge(deld, sample[,c("File.name","Week")], by="File.name")
    deld<-deld[order(deld$Week),]
    deld<-deld[,-1]
    deldM<-melt(deld,id.vars = 'Week')
    colnames(deldM)[2:3]<-c("Position","Freq")
    monkey<-names(monkeys2)[i]
    
    tbweek<-tbs$tb[tbs$ids==monkey]
    ggplot(data=deldM, aes(x=Week, y=Freq))+
        ylab("Insertion frequency")+xlab("Week")+
        facet_wrap(~ Position, nrow=2, ncol=5)+
        geom_point(size=2, color=cols2[1])+theme_bw()+ggtitle(paste(monkey))+
        geom_vline(xintercept=tbweek, col="blue")
    ggsave(paste0("Output/Indels/Timeseries/Del.freq_timeseries_",monkey, ".pdf"), width=7, height=4)
}
