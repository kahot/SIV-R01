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

library(RColorBrewer)
#color for amino acids
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#pie(rep(1,n), col=sample(col_vector, n))
colors<-sample(col_vector, 21)

mycolors<-c("gray90", colors)



AA<-c(".","A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","*")
#Alignment fasta file list
fafiles<-list.files("Output/PID_Con_Alignment/",pattern=".fasta")
fafiles<-list.files("Output/PID_Con_Alignment/",pattern="^Run7")

for (i in 1: length(fafiles)){
    fname<-gsub(".Alignment.fasta","",fafiles[i])
    fasta<-read.dna(paste0("Output/PID_Con_Alignment/",fname, ".Alignment.fasta"),format = "fasta",as.character=TRUE)
    if (is.list(fasta)) fasta<-data.frame(do.call(rbind,fasta))
    if (is.matrix(fasta)) fasta<-data.frame(fasta)
    fasta<-fasta[,-c(1,2)]
    #NT position 217 to 690 (474 NT,158 AA)
    fasta<-fasta[,1:474]
    aa<-data.frame(t(apply(fasta, 1, function(x) seqinr::translate(x))))
    colnames(aa)<-paste0("pos.",(219/3):(690/3))
    write.csv(aa, paste0("Output/AA/", fname,"AA.seq.csv"))
}



reference<-read.dna("Data/AY032751env.fasta", format = "fasta",as.character=TRUE)
ref<-reference[217:690]
refAA<-seqinr::translate(ref)


samples<-read.csv("Data/SamplesNoduplicates.csv",stringsAsFactors = F)
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
tbs<-read.csv("Data/TBinfection.week.csv",stringsAsFactors = F)
monkeyList<-monkeyList[tbs$ids]
Monkeys<-names(monkeyList)

#color assingment
mycolors<-c("gray90", colors)
AAs<-factor(AA, levels=AA)

names(mycolors)<-levels(AAs)

#create haplotype plots with AA
for (i in 1:length(Monkeys)){
    monkeyID<-Monkeys[i]
    sample<-monkeyList[[monkeyID]]
    sample<-sample[order(sample$Week),]
  
    
    Plots<-list()
    for (j in 1:nrow(sample)){
      
        fname<-sample$File.name[j]
        Seq<-read.csv(paste0("Output/AA/", fname,"AA.seq.csv"), stringsAsFactors = F, row.names = 1, colClasses = c("character"))
        if (nrow(Seq)>=60)  Seq<-Seq[1:60,]
        for (k in 1:ncol(Seq)){
          Seq[,k]<-sapply(Seq[,k], function(x) if (x==refAA[k]) x<-"." else x=x)
        }
        
        colnames(Seq)<-1:ncol(Seq)
        Seq$id<-paste0(fname,"_",1:nrow(Seq))
        
        Seqm<-melt(Seq, id.vars = "id")
        colnames(Seqm)[2:3]<-c("Pos","AA")
        Seqm$AA<-factor(Seqm$AA, levels=AAs)
        Seqm$Pos<-as.integer(Seqm$Pos)
        
        Plots[[j]]<-
          ggplot(Seqm, aes(x=Pos,y=id, fill=factor(AA)))+
            geom_tile()+ylab('')+xlab('')+
            theme(axis.text.x = element_text(angle = 90, size=6),
                axis.text.y = element_text(size=6),
                legend.title=element_blank(),
                plot.title = element_text(size=10))+
            scale_colour_discrete(drop=TRUE,limits = levels(AAs))+
            scale_fill_manual(values=mycolors)+
            ggtitle(paste0(monkeyID," Week ",samples$Week[samples$File.name==fname]," (",fname,")" ))+
            scale_x_continuous(breaks=c(27,77,127), labels=c(100,150,200))
    
        }
    pdf(paste0("Output/Haplotypes/AA.",monkeyID,".pdf"), width = 10, height = (4*j))
    do.call(grid.arrange, c(Plots, ncol=1))
    dev.off()

}
    
    
    
    
    



