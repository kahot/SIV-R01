source("Rscripts/baseRscript2.R")


# Trim Run4 R1 files 
samples<-read.csv("Data/SamplesNoduplicates.csv",stringsAsFactors = F)
run4<-samples[samples$Run==4,]
run4files<-run4$File.name

ambnuc=c("y", "r", "w", "s", "k", "m", "b", "d", "h", "v", "n", "x")


#Set cutoff length
k<-100


for (i in 1:length(run4files)){
    fastas<-list.files(paste0("~/programs/PID-master/",run4files[i]), pattern=glob2rx("r1_*.fasta"), recursive = T)
    fname<-run4files[i]
    for (f in 2:length(fastas)){
        r1<-read.fasta(paste0("~/programs/PID-master/",fname,"/",fastas[f]), as.string= TRUE)
        sname<-gsub('.fasta','',fastas[f])
        sname<-gsub("Out/",'',sname)
        #seqn<- data.frame(cbind(seqnam= paste(">", names(seqn), sep= ""), seqn= seqn), stringsAsFactors= FALSE) #
        
        for (j in 1:length(r1)){
            temp<- paste(r1[j])
            temp<-unlist(strsplit(temp,""))      
            temp<-temp[1:(length(temp)-k)]
            temp[which(temp %in% ambnuc)]<-"n"
            temp<-paste(temp, collapse='')
            temp<-toupper(temp)
            
            write.fasta(temp, names(r1)[j], file.out=paste0("~/programs/PID-master/",fname,"/Out/",sname,".trimmed.fasta"), open="a", nbchar = 400)
        }
        
        
    }
    
    
    
