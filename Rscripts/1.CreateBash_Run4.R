library(stringr)

samples<-read.csv("Data/SamplesNoduplicates.csv",stringsAsFactors = F)

run4<-samples[samples$Run==4,]
run4files<-run4$File.name

#read run4 template
temp<-readLines("Data/template/Run4_temp.sh")
    
for (i in 1:length(run4files)){
  fname<-run4files[i]
  numb<-gsub('Run4_','',fname)
  new<-gsub(pattern="Run4_18", replace=paste0(fname),x=temp)
  new<-gsub(pattern="XX", numb, x=new)
  writeLines(new, con=paste0("~/programs/PID-master/Bashscripts/",fname,".sh"))
}

########
## For Run3
run3<-samples[samples$Run==3,]
run3files<-run3$File.name

#read run3 template
temp<-readLines("Data/template/Run3_temp.sh")

for (i in 1:length(run3files)){
    fname<-run3files[i]
    numb<-gsub('Run3_','',fname)
    new<-gsub(pattern="Run4_18", replace=paste0(fname),x=temp)
    new<-gsub(pattern="XX", numb, x=new)
    writeLines(new, con=paste0("~/programs/PID-master/Bashscripts/",fname,".sh"))
}

##Run3
folders<-list.files("~/programs/BaseSpace/B670_Run_3-142079938/FASTQ_Generation_2019-10-03_13_19_45Z-193695834/")

sink("rename_Run3.txt")
for (i in 1: length(folders)){
    nam<-folders[i]
    newn<-substring(nam, 1,2)
    phrase<-paste('mv -n',nam,newn)
    print(phrase)
}
sink(NULL)





## For Run2
folders<-list.files("~/programs/BaseSpace/B670_Run_2-123312015/FASTQ_Generation_2019-03-29_02_49_40Z-171066940/")

sink("rename_Run2.txt")
for (i in 1: length(folders)){
    nam<-folders[i]
    newn<-substring(nam, 1,2)
    phrase<-paste('mv -n',nam,newn)
    print(phrase)
}
sink(NULL)


##
run2<-samples[samples$Run==2,]
run2files<-run2$File.name

#read run2 template
temp<-readLines("Data/template/Run2_temp.sh")

for (i in 1:length(run2files)){
    fname<-run2files[i]
    numb<-gsub('Run2_','',fname)
    new<-gsub(pattern="Run4_18", replace=paste0(fname),x=temp)
    new<-gsub(pattern="XX", numb, x=new)
    writeLines(new, con=paste0("~/programs/PID-master/Bashscripts/",fname,".sh"))
}


### Run0
folders<-list.files("~/programs/BaseSpace/Ambrose1-41042018/FASTQ_Generation_2018-05-21_03_30_27Z-96699097/")
sink("rename_Run0.txt")
for (i in 1: length(folders)){
    nam<-folders[i]
    newn<-substring(nam, 1,2)
    phrase<-paste('mv -n',nam,newn)
    print(phrase)
}
sink(NULL)


##
run0<-samples[samples$Run==0,]
run0files<-run0$File.name

#read run2 template
temp<-readLines("Data/template/Run0_temp.sh")
#dir.create("~/programs/PID-master/Bashscripts2")
for (i in 1:length(run0files)){
    fname<-run0files[i]
    numb<-gsub('Run0_','',fname)
    new<-gsub(pattern="Run4_18", replace=paste0(fname),x=temp)
    new<-gsub(pattern="XX", numb, x=new)
    writeLines(new, con=paste0("~/programs/PID-master/Bashscripts2/",fname,".sh"))
}



