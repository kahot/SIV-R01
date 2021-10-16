library(ggplot2)
#find the unique seqeunces from a fasta file
funiq<-function(all.fa,uniq.fa)
{
    all<-unlist(all.fa)
    sequ.names <- names(uniq.fa)
    sequ <- NULL
    sequ[sequ.names] <- list(NULL)
    for(i in 1: length(uniq.fa))
    {
        sequ[[i]]<-names(all[which(all%in%(uniq.fa[[i]][1]))])
    }
    return(sequ)
}


library(seqinr)
all<-read.fasta("Run4_17_Protein.fasta",as.string = TRUE,seqtype="AA")
uniq<-read.fasta("uniqueRun4_17.fasta",as.string = TRUE,seqtype="AA")
nam<- funiq(all,uniq)
fcount<-lapply(nam,length)
fc<-data.frame(unlist(fcount))
fc$Seq<-rownames(fc)
fc$Seq<-substr(fc$Seq, 1, 11)
colnames(fc)[1]<-"Count"
fc<-fc[order(fc$Count, decreasing=T),]

nonsingle<-fc[!(fc$Count==1|fc$Count==2),]

single<-fc[fc$Count==1,]
twice<-fc[fc$Count==2,]

nonsingle[nrow(nonsingle)+1,]<-c(nrow(single),"Singleton")
nonsingle[nrow(nonsingle)+1,]<-c(nrow(twice),"Twice")
nonsingle$Count<-as.integer(nonsingle$Count)

nonsingle$Seq<-factor(nonsingle$Seq,levels=paste(nonsingle$Seq))

ggplot(data=nonsingle, aes(x=Seq, y=Count))+
    geom_bar(stat="identity")+
    theme(axis.text.x = element_text(angle=90, size=4))
ggsave("Output/Run4_17_uniqSeqcount.pdf", width = 10, height = 3)

all<-read.fasta("Run4_17_DNAalignedt.trimmmed.fasta",as.string = TRUE,seqtype="AA")
uniq<-read.fasta("uniqueRun4_17.fasta",as.string = TRUE,seqtype="AA")
nam<- funiq(all,uniq)
fcount<-lapply(nam,length)
fc<-data.frame(unlist(fcount))
fc$Seq<-rownames(fc)
fc$Seq<-substr(fc$Seq, 1, 11)
colnames(fc)[1]<-"Count"
fc<-fc[order(fc$Count, decreasing=T),]

