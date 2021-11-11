## Non PID ##
# Create bash files to trim and map 'UNPROCESSED' raw fastq files 
library(stringr)

#Specify your Run name
run="Run7_"

#select the directory where fastq files are saved
dir<-"~/programs/BaseSpace/B670_Run_7_Repeat_2-292249961/"
dir2<-"FASTQ_Generation_2021-09-07_14_35_26Z-457696525/"

### !!!! 
#If the directory names are not simplified yet by running the PID pipeline (1.CreateBash_PID):
folders<-list.files(paste0(dir,dir2))
sink(paste0(dir,dir2,"rename_",run,".sh"))
for (i in 1: length(folders)){
  nam<-folders[i]
  newn<-substring(nam, 1,2)
  phrase<-paste('mv -n',nam,newn)
  cat(phrase)
  cat("\n")
}
sink(NULL)

#cd to the right directory & run 'bash rename.RunXX.sh'




## Step 1

dir<-"/Volumes/Kaho_Data/SIV_DATA/B670_Run_6-268458190/"
dir2<-"FASTQ_Generation_2021-06-04_22_44_27Z-424462038/"



# List all files to be processed
fq<-list.files(paste0(dir,dir2), pattern="fastq.gz$",recursive = T) 

#create a list of file names:
n<-seq(1, by = 2, len = (length(fq)/2))
fq2<-fq[n]
flist<-data.frame(matrix(nrow=length(n), ncol=1))
colnames(flist)<-"File.name"
for (i in 1:length(fq2)){
  #choose the paired reads fastq files
  fa1<-fq2[i]
  fname1<-sub(".*/", "", fa1)
  fname<-paste0(run, substr(fname1, 1,2))
  flist[i,1]<-fname
}
#write.csv(flist,paste0("Data/",run,"filenames.csv"))

runfiles<-flist$File.name

# Create bash files to run bbmap and bwa

# 1. Read the template text file:
temp<-readLines("Data/template/Bashxx.sh")

#dir.create("Data/Bashscripts/")

for (i in 1:length(runfiles)){
  #choose the paired reads fastq files
  fname<-runfiles[i]
  numb<-gsub(run,'',fname)
  
  fa1<-fq2[i]
  fa1<-sub(".*/", "", fa1)
  fa2<-gsub(pattern="R1",replace="R2",x=fa1)
  
  new<-gsub(pattern="dir1/dir2/", replace=paste0(dir,dir2),x=temp)
  new<-gsub(pattern="XX", numb, x=new)
  new<-gsub(pattern="10_S10_L001_R1_001.fastq", replace=paste0(fa1),x=new)
  new<-gsub(pattern="10_S10_L001_R2_001.fastq", replace=paste0(fa2),x=new)
  new<-gsub(pattern="M10",replace=paste0(fname),x=new)
  
  writeLines(new, con=paste0("Data/Bashscripts/",fname,".sh"))

}

############
### Run6 -filter and map with relaxed setting 

#Specify your Run name
run="Run6_"

#select the directory where fastq files are saved
dir<-"/Volumes/Kaho_Data/SIV_DATA/B670_Run_6-268458190/"
dir2<-"FASTQ_Generation_2021-06-04_22_44_27Z-424462038/"


# List all files to be processed
fq<-list.files(paste0(dir,dir2), pattern="fastq.gz$",recursive = T) 

#create a list of file names:
n<-seq(1, by = 2, len = (length(fq)/2))
fq2<-fq[n]
flist<-data.frame(matrix(nrow=length(n), ncol=1))
colnames(flist)<-"File.name"
for (i in 1:length(fq2)){
  #choose the paired reads fastq files
  fa1<-fq2[i]
  fname1<-sub(".*/", "", fa1)
  fname<-paste0(run, substr(fname1, 1,2))
  flist[i,1]<-fname
}

runfiles<-flist$File.name

# Create bash files
# Read the template text file:
temp<-readLines("Data/template/Bash_QC.sh")

#dir.create("Data/Bashscripts/")

for (i in 1:length(runfiles)){
  #choose the paired reads fastq files
  fname<-runfiles[i]
  numb<-gsub(run,'',fname)
  
  fa1<-fq2[i]
  fa1<-sub(".*/", "", fa1)
  fa2<-gsub(pattern="R1",replace="R2",x=fa1)
  
  new<-gsub(pattern="dir1/dir2/", replace=paste0(dir,dir2),x=temp)
  new<-gsub(pattern="XX", numb, x=new)
  new<-gsub(pattern="10_S10_L001_R1_001.fastq", replace=paste0(fa1),x=new)
  new<-gsub(pattern="10_S10_L001_R2_001.fastq", replace=paste0(fa2),x=new)
  new<-gsub(pattern="M10",replace=paste0(fname),x=new)
  
  writeLines(new, con=paste0("Data/Bashscripts/",fname,".qc.sh"))
  
}

#### New pipeline creation: filter the reads before running with PID Pipeline 

#Specify your Run name
run="Run6_"

#select the directory where fastq files are saved
dir<-"/Volumes/Kaho_Data/SIV_DATA/B670_Run_6-268458190/"
dir2<-"FASTQ_Generation_2021-06-04_22_44_27Z-424462038/"


# List all files to be processed
fq<-list.files(paste0(dir,dir2), pattern="fastq.gz$",recursive = T) 

#create a list of file names:
n<-seq(1, by = 2, len = (length(fq)/2))
fq2<-fq[n]
flist<-data.frame(matrix(nrow=length(n), ncol=1))
colnames(flist)<-"File.name"
for (i in 1:length(fq2)){
  #choose the paired reads fastq files
  fa1<-fq2[i]
  fname1<-sub(".*/", "", fa1)
  fname<-paste0(run, substr(fname1, 1,2))
  flist[i,1]<-fname
}

runfiles<-flist$File.name

# Create bash files
# Read the template text file:
temp<-readLines("Data/template/Bash_Trim.sh")

for (i in 1:length(runfiles)){
  #choose the paired reads fastq files
  fname<-runfiles[i]
  numb<-gsub(run,'',fname)
  
  fa1<-fq2[i]
  fa1<-sub(".*/", "", fa1)
  fa2<-gsub(pattern="R1",replace="R2",x=fa1)
  
  new<-gsub(pattern="dir1/dir2/", replace=paste0(dir,dir2),x=temp)
  new<-gsub(pattern="XX", numb, x=new)
  new<-gsub(pattern="10_S10_L001_R1_001.fastq", replace=paste0(fa1),x=new)
  new<-gsub(pattern="10_S10_L001_R2_001.fastq", replace=paste0(fa2),x=new)
  new<-gsub(pattern="Run6_1_",replace=paste0(fname),x=new)
  
  writeLines(new, con=paste0("~/programs/PID-master/Bashscripts/",fname,".2.sh"))
  
}


