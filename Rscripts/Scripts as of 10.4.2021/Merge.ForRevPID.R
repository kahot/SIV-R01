library(seqinr)
library(msa)
library(ape)
library(plyr)
library(stringr)

### Read the mapped sam file
fList<-list.files("~/programs/PID-master/Sam/",recursive = T, pattern="clipped.sam") 
coln<-c('QNAME','Flag','RefName','Pos','MapQ','cigar','MRNM','Mpos','isize','seq','Qual','tag1','tag2','tag3','tag4','tag5','tag6')

for ( i in 1:length(fList)){
    sam<-read.table(paste0("~/programs/PID-master/Sam/",fList[i]),skip=3, col.names=coln, sep = "\t",fill=T, comment.char="",quote= "", stringsAsFactors = F)
    sam<-sam[,1:11]
    while (grepl("@", sam[1,1], fixed = TRUE)){
        sam<-sam[-1,]
    }
    
    fname<-substr(fList[i],1,7)    
    
    sqnames<-unique(sam$QNAME)
    
    for (j in 1:length(sqnames)){
        seqs<-sam[sam$QNAME==sqnames[j],]
        
        #Reverse and Forward sequences
        f1<-unlist(strsplit(seqs$seq[1], split = ""))
        r1<-unlist(strsplit(seqs$seq[2], split = ""))
        
        #Check indels for forward 
        cigF<-seqs$cigar[1]
        gap=0
        if (grepl("I|D",cigF)) {
            if (grepl("I",cigF)){
                ins<-unlist(strsplit(cigF,"I"))
                ins<-head(ins,-1)
                for (k in 1:length(ins)){
                    m<-str_extract(ins[k], "(\\d+$)")
                    gap<-gap+as.numeric(m)
                }
            }
            if (grepl("D",cigF)){
                del<-unlist(strsplit(cigF,"D"))
                del<-head(del,-1)
                for (k in 1:length(del)){
                    d<-str_extract(del[k], "(\\d+$)")
                    gap<-gap-as.numeric(d)
                }
            } 
        }       
        
        
        startF<-as.integer(seqs$Pos[1]) #start pos F
        startR<-as.integer(seqs$Pos[2]) #start pos R
        end=startF+length(f1)-1  #ending pos of F
        end=end-gap
        #cat(j," ", "Overlap? ",end>=startR)
        
        if (end>=startR){
            #overlapping section
            ovF<-f1[(startR-startF+1+gap):length(f1)]
            ovR<-r1[1:length(ovF)]
            if (identical(ovF,ovR))  Cons<-c(f1,r1[(length(ovF)+1):length(r1)])
            else {
                #cat("Sequence=", j ," ")
                #print(ovF==ovR)
                #print(ovF)
                #print(ovR)
                Cons<-c(f1,r1[(length(ovF)+1):length(r1)])
            }
        }
        if (end<startR){
            n<-startR-end-1
            Cons<-c(f1,rep("N", times=n),r1)
            
        }
        
        #remove the sequences that have majority Ns
        n<-length(Cons[Cons=="N"])
        if (n<150){
            write.fasta(Cons, sqnames[j], file.out=paste0("Output/PID_Consensus/",fname, ".consensus.fasta"), open="a", nbchar = 700)
        }
        
    }
    print(fname)
}





