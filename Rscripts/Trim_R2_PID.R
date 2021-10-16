

fastas<-list.files("Output/PID_nonmerged/", pattern=".fasta")
i=1

ambnuc=c("y", "r", "w", "s", "k", "m", "b", "d", "h", "v", "n", "x")


#cutoff length
k<-130

for (i in 1:length(fastas)){
    r2<-read.fasta(paste0("Output/PID_nonmerged/",fastas[i]), as.string= TRUE)
    fname<-gsub('.fasta','',fastas[i])
    seqn<- data.frame(cbind(seqnam= paste(">", names(seqn), sep= ""), seqn= seqn), stringsAsFactors= FALSE) #
    
    for (j in 1:length(r2)){
        temp<- paste(r2[j])
        temp<-unlist(strsplit(temp,""))      
        temp<-temp[(k+1):length(temp)]
        temp[which(temp %in% ambnuc)]<-"n"
        temp<-paste(temp, collapse='')
        temp<-toupper(temp)
        
        write.fasta(temp, names(r2)[j], file.out=paste0("Output/PID_nonmerged/",fname,".trimmed.fasta"), open="a", nbchar = 400)
    }
    
}
seqn<-read.fasta(paste0("Output/PID_nonmerged/",fastas[i]), as.string= TRUE)
r2s<-read.fasta(paste0("Output/PID_nonmerged/",fastas[i]), as.string= TRUE)

j=1
ambcode=c("y", "r", "w", "s", "k", "m", "b", "d", "h", "v", "n", "x")
for (i in )
seq<-paste0(r2[j])
seqn<- data.frame(cbind(seqnam= paste(">", names(seqn), sep= ""), seqn= seqn), stringsAsFactors= FALSE) #
