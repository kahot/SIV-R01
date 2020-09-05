library(stringr)

#### 9/3/20 ####
# Read files directly from Basespace directory
fq<-list.files("~/programs/BaseSpace/B670_Run_3-142079938/", pattern="fastq.gz$",recursive = T) 
fq<-list.files("~/programs/BaseSpace/B670_Run_2-123312015/", pattern="fastq.gz$",recursive = T) 
fq<-list.files("~/programs/BaseSpace/Ambrose_B670-2-77614630/", pattern="fastq.gz$",recursive = T) 

#create a list of file names:
n<-seq(1, by = 2, len = (length(fq)/2))
fq2<-fq[n]
flist<-data.frame(matrix(nrow=23, ncol=1))
for (i in 1:length(fq2)){
    #choose the paired reads fastq files
    fa1<-fq2[i]
    fname1<-sub(".*/", "", fa1)
    fname<-paste0("Run1_", substr(fname1, 1,2))
    flist[i,1]<-fname
}
write.csv(flist,"Data/Run1_filenames.csv")
write.csv(flist,"Data/Run3_filenames.csv")
write.csv(flist,"Data/Run2_filenames.csv")

##### 
#Select Monkey  A22317 files

samples<-read.csv("Data/Samples.csv")

# number of unique monkeys
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
monkeys<-names(monkeyList)

samples2<-samples[samples$Monkey %in% monkeys,]

fileNames<-as.character(samples2$File.name)
fileNames<-gsub('$','.sh', fileNames)
fileConn<-file("Data/monkeyNames.txt")
writeLines(fileNames, fileConn)
close(fileConn)


write.table(fileNames, "Data/mokeyNames.txt")

#create the bash files to run bbmap and bwa
# read the template command text file:
cmmd<-readLines("Data/template/Bashxx.sh")

#choose the fastq files to be prrocessed
fq<-list.files("~/programs/BaseSpace/", pattern="fastq.gz$",recursive = T) 


#dir.create("Data/Bashscripts/")

#create vector of odd numbers:
n<-seq(1, by = 2, len = (length(fq)/2))
fq2<-fq[n]
for (i in 1:length(fq2)){
  #choose the paired reads fastq files
  fa1<-fq2[i]
  fa2<-gsub(pattern="R1",replace="R2",x=fa1)
  fname1<-sub(".*/", "", fa1)
  fname2<-substr(fa1,start=1,stop=10)
  if (fname2=="Ambrose_B6") fname<- paste0("Run1_", substr(fname1, 1,2))
  if (fname2=="B670_Run_2") fname<- paste0("Run2_", substr(fname1, 1,2))
  if (fname2=="B670_Run_3") fname<- paste0("Run3_", substr(fname1, 1,2))
  
  new<-gsub(pattern="10_S10_L001_R1_001.fastq", replace=paste0(fa1),x=cmmd)
  new<-gsub(pattern="10_S10_L001_R2_001.fastq", replace=paste0(fa2),x=new)
  new<-gsub(pattern="M10",replace=paste0(fname),x=new)
  writeLines(new, con=paste0("Data/Bashscripts/",fname,".sh"))

}





##### Map raw R2 to SIV2 to look for Primer IDs
tmp<-readLines("Data/template/map_rawFastq.sh")
fq<-list.files("Data/",recursive = T, pattern="R2_001.fastq.gz$") 
map<-c()
for (i in 1:length(fq)){
    fname1<-sub(".*/", "", fq[i])
    fname<-paste0(substr(fq[i],start=6,stop=10), "_", substr(fname1, 1,2))
    new<-gsub(pattern="10_S10_L001_R2_001.fastq", replace=fq[i],x=tmp)
    new<-gsub(pattern="Run_3_9",replace=paste0(fname),x=new)
    write(new,file="Data/Bashscripts/map.sh",append=TRUE)
   # writeLines(new, con=paste0("Data/Bashscripts/map.",fname,".sh"))
    
}


### run samtools -sort and index bam
cmmd<-readLines("Data/template/samtool_run.sh")
f<-fq<-list.files("Data/Bash2/") 
map<-c()
for (i in 1:length(f)){
    fname<-substr(f[i],start=1,stop=7)
    new<-gsub(pattern="Run1_2_",replace=paste0(fname),x=cmmd)
    write(new,file="Data/Bash2/map.sh",append=TRUE)
    
}
