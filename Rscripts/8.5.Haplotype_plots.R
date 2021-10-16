library(ggplot2)
library(reshape)
library(colorspace)
library(msa)
library(ape)
library(seqinr)
library(bios2mds)
library(DataCombine)
library(gridExtra)
source("Rscripts/baseRscript2.R")
cols2<-qualitative_hcl(6, palette="Dark3")


#Alignment fasta file list
fafiles<-list.files("Output/PID_Con_Alignment/",pattern=".fasta")

reference<-read.dna("Data/AY032751env.fasta", format = "fasta",as.character=TRUE)
ref<-reference[215:691]

samples<-read.csv("Data/SamplesNoduplicates.csv",stringsAsFactors = F)
samples$Week<-as.integer(samples$Week)
samples$Tissue<-factor(samples$Tissue, levels=c("Stock","Plasma","LN","HLN","Unknown"))

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

monkeyList<-monkeyList[tbs$ids]
Monkeys<-names(monkeyList)

for (i in 1:length(Monkeys)){
for (i in 1:5){
    monkeyID<-Monkeys[i]
    sample<-monkeyList[[monkeyID]]
    sample<-sample[order(sample$Week),]
    
    tbweek<-tbs$tb[tbs$ids==monkeyID]
    #monkey<-gsub("A",'',monkey)
    Plots<-list()
    for (j in 1:nrow(sample)){
        fname<-sample$File.name[j]
        fasta<-read.dna(paste0("Output/PID_Con_Alignment/",fname, ".Alignment.fasta"),format = "fasta",as.character=TRUE)
        if (is.list(fasta)) fasta<-data.frame(do.call(rbind,fasta))
        if (is.matrix(fasta)) fasta<-data.frame(fasta)
        
        #Take first 60 seqeunces
        if (nrow(fasta)>=60) fasta<-fasta[1:60,]

        Seq<-InsertRow(fasta,ref,1)
        #fa$ID<-paste0(fname,"_",1:nrow(fasta))
        rownames(Seq)<-c("ref", paste0(fname,"_",1:nrow(fasta)))
        
        for (k in 1:ncol(Seq)){
          Seq[2:nrow(Seq),k]<-sapply(Seq[2:nrow(Seq),k], function(x) if (x==Seq[1,k]) x<-"." else x=x)
        }
        Seq$id<-rownames(Seq)
        
        colnames(Seq)<-c(1:477,"id")
        seq2<-Seq[-1,]
        
        #seq2<-Seq[121:180,]
        
        
        Seqm<-melt(seq2, id.vars = "id")
        colnames(Seqm)[2:3]<-c("Pos","Nuc")
        nuc<-unique(Seqm$Nuc)
        if (!("-" %in% nuc)) colors<-c("gray90",cols)
        else colors<-c("white","gray90",cols)
        Seqm$Nuc<-factor(Seqm$Nuc, levels=c("-",".","a","c","g","t"))
        Seqm$Pos<-as.integer(Seqm$Pos)
        Plots[[j]]<-ggplot(Seqm, aes(x=Pos,y=id, fill=factor(Nuc)))+
          geom_tile()+ylab('')+xlab('')+
          theme(axis.text.x = element_text(angle = 90, size=6),
                axis.text.y = element_text(size=6),
                legend.title=element_blank(),
                plot.title = element_text(size=10))+
          #geom_vline(xintercept=c(.5,(seq(3, ncol(mutF), 3)+.5)), color="gray50")+
          scale_fill_manual(values=colors)+
          ggtitle(paste0(monkeyID," Week ",samples$Week[samples$File.name==fname]," (",fname,")" ))+
          scale_x_continuous(breaks=c(0,85,185,285,385,485), labels=c(215,300,400,500,600,700))
        
        #ggsave(paste0("Output/Haplotypes/",monkeyID,"_",fname,"..pdf"), height = 4, width = 10)
        
    }
    pdf(paste0("Output/Haplotypes/",monkeyID,".pdf"), width = 10, height = (4*j))
    do.call(grid.arrange, c(Plots, ncol=1))
    dev.off()
    
    
}
        
#Stock Figure       
fasta<-read.dna(paste0("Output/PID_Con_Alignment/Run0_17.Alignment.fasta"),format = "fasta",as.character=TRUE)
if (is.matrix(fasta)) fasta<-data.frame(fasta)

#remove the sequences with N
fa<-fasta[rowSums(fasta == "n")==0, , drop = FALSE]
#Take first 60 seqeunces
if (nrow(fa)>=60) fa<-fa[1:60,]
Seq<-InsertRow(fa,ref,1)
rownames(Seq)<-c("ref", paste0("Stock_",1:nrow(fa)))
  
for (k in 1:ncol(Seq)){
    Seq[2:nrow(Seq),k]<-sapply(Seq[2:nrow(Seq),k], function(x) if (x==Seq[1,k]) x<-"." else x=x)
}
Seq$id<-rownames(Seq)
  
colnames(Seq)<-c(1:477,"id")
seq2<-Seq[-1,]

Seqm<-melt(seq2, id.vars = "id")
colnames(Seqm)[2:3]<-c("Pos","Nuc")
nuc<-unique(Seqm$Nuc)
if (!("-" %in% nuc)) colors<-c("gray90",cols)
if (("-" %in% nuc)) colors<-c("white","gray90",cols)
Seqm$Nuc<-factor(Seqm$Nuc, levels=c("-",".","a","c","g","t"))
Seqm$Pos<-as.integer(Seqm$Pos)
ggplot(Seqm, aes(x=Pos,y=id, fill=factor(Nuc)))+
    geom_tile()+ylab('')+xlab('')+
    theme(axis.text.x = element_text(angle = 90, size=6),
          axis.text.y = element_text(size=6),
          legend.title=element_blank(),
          plot.title = element_text(size=10))+
    #geom_vline(xintercept=c(.5,(seq(3, ncol(mutF), 3)+.5)), color="gray50")+
    scale_fill_manual(values=colors)+
    ggtitle( "Stock")+
    scale_x_continuous(breaks=c(0,85,185,285,385,485), labels=c(215,300,400,500,600,700))
ggsave(paste0("Output/Haplotypes/Stock.pdf"), height = 4, width = 10)
  
