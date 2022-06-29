library(ggplot2)
library(reshape)
library(colorspace)
library(msa)
library(ape)
library(seqinr)
library(bios2mds)

source("Rscripts/baseRscript2.R")
source("Rscripts/alignment2Fasta.R")
#hcl_palettes(plot = TRUE)
cols2<-qualitative_hcl(6, palette="Dark3")
#cols2<-qualitative_hcl(5,palette="Set2")


# read the files saved in Overview_output:
samfiles<-list.files("Output/sam/",pattern=".sam")
coln<-c('QNAME','Flag','RefName','Pos','MapQ','cigar','MRNM','Mpos','isize','seq','Qual','tag1','tag2','tag3','tag4','tag5','tag6')

reference<-read.dna("Data/AY032751.fasta", format = "fasta",as.character=TRUE)
referF<-reference[1:480]
referF<-paste0(toupper(referF), collapse = '')
referR<-reference[400:720]
referR<-paste0(toupper(referR), collapse = '')

#variable positions
vars<-read.csv("Output/HighMutfreq_sites_all.csv",row.names = 1,stringsAsFactors = F)
#select the mutations that are present at least in 2 animals
vars$eval<-apply(vars[2:8],1, function(x) if(length(x[x=="Y"])>=2) 1 else 0)
vars<-vars[vars$eval==1,]

#create a vector of codons of all mutations
positions<-c()
for (i in 1:nrow(vars)){
    if (vars$codon[i]==1) codon<-c(vars$pos[i], vars$pos[i]+1, vars$pos[i]+2)
    if (vars$codon[i]==2) codon<-c(vars$pos[i]-1, vars$pos[i], vars$pos[i]+1)
    if (vars$codon[i]==3) codon<-c(vars$pos[i]-2, vars$pos[i]-1, vars$pos[i])
    positions<-c(positions, codon)
}
positions<-unique(positions)

### sample list
samples<-read.csv("Data/Samples.csv",stringsAsFactors = F)
ids<-c("A21918","A22317","A22517", "A22617","A22217","A23918","A22117")

for (id in 1:length(ids)){
    monkey<-ids[id]
    monkey1<-samples$File.name[samples$Monkey==monkey]
    files<-substr(samfiles, 1, 7)
    Is<-which(files %in% monkey1)
    #dir.create(paste0("Output/reads/",monkey,"/"))
    
    for (k in 1:length(Is)){
        i<-Is[k]
        sam<-read.table(paste0("Output/sam/",samfiles[i]),skip=3, col.names=coln, sep = "\t",fill=T, comment.char="",quote= "", stringsAsFactors = F)
        sam<-sam[,1:11]
        fname<-substr(samfiles[i],1,7)
        print(fname)
        sam<-subset(sam, MapQ>9&MapQ<61) 
    #    samF<-sam[sam$Pos<200,]
        samR<-sam[sam$Pos>400,]
    #    seqF<-samF[sample(nrow(samF), 20),]
        seqR<-samR[sample(nrow(samR), 20),]
    #    
    #    readsF<-seqF$seq
    #    readsF<-c(referF,readsF)
    #    dnaF<-DNAStringSet(readsF)
    #    names(dnaF)<-c("ref", paste0(fname,"_F_",1:20))
    #    alignF<-msa(dnaF)
    #    alignment2Fasta(alignF, paste0("Output/reads/", monkey,"/", fname, "_F_.fasta"))
    #    
        readsR<-seqR$seq
        readsR<-c(referR,readsR)
        dnaR<-DNAStringSet(readsR)
        names(dnaR)<-c("ref", paste0(fname,"_R_",1:20))
        alignR<-msa(dnaR)
        alignment2Fasta(alignR, paste0("Output/reads/", monkey,"/", fname, "_R_.fasta"))
    }
}
    


for (id in 2:length(ids)){
    monkey<-ids[id]

    #read the fasta files
    Forward<-list.files(paste0("Output/reads/", monkey,"/"), pattern="F_.fasta")
    Rev<-list.files(paste0("Output/reads/", monkey,"/"), pattern="R_.fasta")
    
    for (m in 1:length(Forward)){
        fname<-substr(Forward[m],1,7)
        readsF<-read.dna(paste0("Output/reads/",monkey,"/", Forward[m]), format = "fasta",as.character=TRUE)
        readsF<-data.frame(readsF, stringsAsFactors = F)
        ins<-which(readsF["ref",]=="-")
        k=1
        j=1
        for (i in 1:ncol(readsF)){
            if (!i %in% ins) {
                colnames(readsF)[i]<-k
                k<-k+1}
            if (i %in% ins){
                colnames(readsF)[i]<-paste0("g",j)
                j=j+1}
        }
        
        #readsR<-read.dna(paste0("Output/reads/",monkey,"/", Rev[m]), format = "fasta",as.character=TRUE)
        #readsR<-data.frame(readsR, stringsAsFactors = F)
        #ins<-which(readsR["ref",]=="-")
        #k=400
        #j=1
        #for (i in 1:ncol(readsR)){
        #    if (!i %in% ins) {
        #        colnames(readsR)[i]<-k
        #        k<-k+1}
        #    if (i %in% ins){
        #        colnames(readsR)[i]<-paste0("gr",j)
        #        j=j+1}
        #}
        
        posF<-positions[positions<=480]
        mutF<-readsF[,paste0(posF)]
        refF<-mutF["ref",]
        mutF<-mutF[-which(rownames(mutF)=="ref"),]
        mutF<-rbind(refF,mutF) 
        for (i in 1:ncol(mutF)){
            mutF[2:nrow(mutF),i]<-sapply(mutF[2:nrow(mutF),i], function(x) if (x==mutF[1,i]) x<-"." else x=x)
        }
        mutF$id<-rownames(mutF)
        mutFm<-melt(mutF, id.vars = "id")
        colnames(mutFm)[2:3]<-c("Pos","Nuc")
        mutFm$Nuc<-factor(mutFm$Nuc, levels=c("-",".","a","c","g","t"))
        ggplot(mutFm, aes(x=Pos,y=id, fill=factor(Nuc)))+
            geom_tile(color="gray70")+ylab('')+xlab('')+
            theme(axis.text.x = element_text(angle = 90, size=6),
                  axis.text.y = element_text(size=6),
                  legend.title=element_blank(),
                  plot.title = element_text(size=10))+
            geom_vline(xintercept=c(.5,(seq(3, ncol(mutF), 3)+.5)), color="gray50")+
            scale_fill_manual(values=c("white","gray90",cols2))+
            ggtitle(paste0(monkey," Week ",samples$Week[samples$File.name==fname]," (",fname,")" ))
        ggsave(paste0("Output/reads/",monkey,"/",fname,"F.pdf"), height = 3, width = 5)
        
        
        #posR<-positions[positions>=480]
        #mutR<-readsR[,paste(posR)]
        #refR<-mutR["ref",]
        #mutR<-mutR[-which(rownames(mutR)=="ref"),]
        #mutR<-rbind(refR,mutR) 
        #for (i in 1:ncol(mutR)){
        #    mutR[2:nrow(mutR),i]<-sapply(mutR[2:nrow(mutR),i], function(x) if (x==mutR[1,i]) x<-"." else x=x)
        #}
        #
        #mutR$id<-rownames(mutR)
        #mutRm<-melt(mutR, id.vars = "id")
        #colnames(mutRm)[2:3]<-c("Pos","Nuc")
        #mutRm$Nuc<-factor(mutRm$Nuc, levels=c("-",".","a","c","g","t"))
        #ggplot(mutRm, aes(x=Pos,y=id, fill=factor(Nuc)))+
        #    geom_tile(color="gray70")+ylab('')+xlab('')+
        #    theme(axis.text.x = element_text(angle = 90, size=6),
        #          axis.text.y = element_text(size=6),
        #          legend.title=element_blank(),
        #          plot.title = element_text(size=10))+
        #    geom_vline(xintercept=c(.5,(seq(3, ncol(mutR), 3)+.5)), color="gray50")+
        #    scale_fill_manual(values=c("white","gray90",cols2))+
        #    ggtitle(paste0(monkey," Week ",samples$Week[samples$File.name==fname]," (",fname,")" ))
        #ggsave(paste0("Output/reads/",monkey,"/",fname,"R.pdf"), height = 3, width = 4.5)
        
    }
}   
    




##########    
## A21918 reverse didn't align well:  adjust the start & try
#for (k in 1:length(Is)){ k=1,3,6

monkey=ids[6]
#monkey=ids[5]
monkey1<-samples$File.name[samples$Monkey==monkey]
files<-substr(samfiles, 1, 7)
Is<-which(files %in% monkey1)

referR2<-reference[400:720]
referR2<-paste0(toupper(referR2), collapse = '')
k=6
for (k in 1:length(Is)){
    i<-Is[k]
    sam<-read.table(paste0("Output/sam/",samfiles[i]),skip=3, col.names=coln, sep = "\t",fill=T, comment.char="",quote= "", stringsAsFactors = F)
    sam<-sam[,1:11]
    fname<-substr(samfiles[i],1,7)
    print(fname)
    sam<-subset(sam, MapQ>9&MapQ<61) 
    samR<-sam[sam$Pos>400,]
    #seqF<-samF[sample(nrow(samF), 20),]
    seqR<-samR[sample(nrow(samR), 20),]
    
    readsR<-seqR$seq
    readsR<-c(referR2,readsR)
    dnaR<-DNAStringSet(readsR)
    names(dnaR)<-c("ref", paste0(fname,"_R_",1:20))
    alignR<-msa(dnaR)
    alignment2Fasta(alignR, paste0("Output/reads/", monkey,"/", fname, "_R_.fasta"))
}
Rev<-list.files(paste0("Output/reads/", monkey,"/"), pattern="22_R_.fasta")
    
for (m in 1:length(Rev)){
    readsR<-read.dna(paste0("Output/reads/", monkey,"/",Rev[m]), format = "fasta",as.character=TRUE)
    readsR<-data.frame(readsR, stringsAsFactors = F)
    fname<-substr(Rev[m],1,7)
    ins<-which(readsR["ref",]=="-")
    k=400
    j=1
    for (i in 1:ncol(readsR)){
        if (!i %in% ins) {
            colnames(readsR)[i]<-k
            k<-k+1}
        if (i %in% ins){
            colnames(readsR)[i]<-paste0("gr",j)
            j=j+1}
    }
    
    posR<-positions[positions>=480]
    mutR<-readsR[,paste(posR)]
    refR<-mutR["ref",]
    mutR<-mutR[-which(rownames(mutR)=="ref"),]
    mutR<-rbind(refR,mutR) 
    for (i in 1:ncol(mutR)){
        mutR[2:nrow(mutR),i]<-sapply(mutR[2:nrow(mutR),i], function(x) if (x==mutR[1,i]) x<-"." else x=x)
    }
    
    mutR$id<-rownames(mutR)
    mutRm<-melt(mutR, id.vars = "id")
    colnames(mutRm)[2:3]<-c("Pos","Nuc")
    mutRm$Nuc<-factor(mutRm$Nuc, levels=c("-",".","a","c","g","t"))
    ggplot(mutRm, aes(x=Pos,y=id, fill=factor(Nuc)))+
        geom_tile(color="gray70")+ylab('')+xlab('')+
        theme(axis.text.x = element_text(angle = 90, size=6),
              axis.text.y = element_text(size=6),
              legend.title=element_blank(),
              plot.title = element_text(size=10))+
        geom_vline(xintercept=c(.5,(seq(3, ncol(mutR), 3)+.5)), color="gray50")+
        scale_fill_manual(values=c("white","gray90",cols2))+
        ggtitle(paste0(monkey," Week ",samples$Week[samples$File.name==fname]," (",fname,")" ))
    ggsave(paste0("Output/reads/",monkey,"/",fname,"R.pdf"), height = 3, width = 4.5)
    
    }
    
    
    
