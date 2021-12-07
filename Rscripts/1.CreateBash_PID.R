library(stringr)

# Process Raw FASTAQ files

#Specify your Run name
run="Run7_"

#select the directory where fastq files are saved
dir<-"~/programs/BaseSpace/B670_Run_7_Repeat_2-292249961/FASTQ_Generation_2021-09-07_14_35_26Z-457696525/"


#Other runs
run="Run6_"
dir<-"/Volumes/Kaho_Data/SIV_DATA/B670_Run_6-268458190/FASTQ_Generation_2021-06-04_22_44_27Z-424462038/"

run="Run5_"
dir<-"/Volumes/Kaho_Data/SIV_DATA/B670_run_5_and_Mouse_LAI_RT-215809595/FASTQ_Generation_2020-12-12_02_55_33Z-352361009/"

run="Run4_"
dir<-"/Volumes/Kaho_Data/SIV_DATA/Youya_Comp_Assay_and_B670_repeats-197798601/FASTQ_Generation_2020-09-21_12_33_39Z-317066848/"

run="Run3_"
dir<-"/Volumes/Kaho_Data/SIV_DATA/B670_Run_3-142079938/FASTQ_Generation_2019-10-03_13_19_45Z-193695834/"

run="Run2_"
dir<-"/Volumes/Kaho_Data/SIV_DATA/B670_Run_2-123312015/FASTQ_Generation_2019-03-29_02_49_40Z-171066940/"

run="Run0_"
dir<-"/Volumes/Kaho_Data/SIV_DATA/Ambrose1-41042018/FASTQ_Generation_2018-05-21_03_30_27Z-96699097/"



##################################################
## Step 1: Rename the folders for fastq files (simplify)

folders<-list.files(paste0(dir))

sink(paste0(dir,"rename_",run,".sh"))
for (i in 1: length(folders)){
    nam<-folders[i]
    newn<-substring(nam, 1,2)
    phrase<-paste('mv -n',nam,newn)
    cat(phrase)
    cat("\n")
}
sink(NULL)

#cd to the right directory & run 'bash rename.RunXX.sh'

##################################################
## Step 2:Process FASTQ files using PID pipeline + BBTools + BWA 

#Read the file names
fq<-list.files(dir, pattern="fastq.gz$",recursive = T) 
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

#read the appropriate template
temp<-readLines("Data/template/Run_TrimPID.sh")

#for Run0 and Run2, use Run_TrimPIDrun2.sh (different PRIMER ID/TCS file)
temp<-readLines("Data/template/Run_TrimPIDrun2.sh")

for (i in 1:length(runfiles)){
    #choose the paired reads fastq files
    fname<-runfiles[i]
    numb<-gsub(run,'',fname)
    
    fa1<-fq2[i]
    fa1<-sub(".*/", "", fa1)
    fa2<-gsub(pattern="R1",replace="R2",x=fa1)
    
    new<-gsub(pattern="dir/", replace=paste0(dir),x=temp)
    new<-gsub(pattern="XX", numb, x=new)
    new<-gsub(pattern="10_S10_L001_R1_001.fastq", replace=paste0(fa1),x=new)
    new<-gsub(pattern="10_S10_L001_R2_001.fastq", replace=paste0(fa2),x=new)
    new<-gsub(pattern="Run6_1_",replace=paste0(fname),x=new)
    
    writeLines(new, con=paste0("~/programs/PID-master/Bashscripts/",fname,".sh"))
    
}

## Step 2.2: Run the bash scripts in PID-master directory
#-------
#for f in Bashscripts/*.sh; do
#bash "$f" -H   || break 
#done
#-------


##################################################
## Step 3: Run Merge.ForRevPID.R first. This step is for post merging 
# Map PID-consensus fasta to reference using bwa (w/ relaxed setting)

files<-list.files("Output/PID_Consensus/", pattern = "Run3_")
#files<-list.files("Output/PID_Consensus/", pattern = ".fasta")

temp<-readLines("Data/template/Bash_mapConPID.sh")
for (i in 1:length(files)){
    fname<-files[i]
    fname<-gsub(".consensus.fasta",'',fname)
    new<-gsub(pattern="Run4_17", replace=paste0(fname),x=temp)
    writeLines(new, con=paste0("Data/Bash_PID/",fname,".sh"))
}

# Run the bash files in SIV-R01/ directory
#for f in Data/Bash_PID/*.sh; do
#bash "$f" -H   || break 
#done

#Outputs in bam_PID (bam and bam.bai)





#####  Scripts to rename fasta files ####
#files<-list.files("Output/PID_Consensus/", pattern = ".fa$")
#
#sink("Output/PID_Consensus/renameFasta.sh")
#for (i in 1: length(files)){
#    nam<-files[i]
#    newn<-gsub(".fa$","",nam)
#    phrase<-paste0('mv -n ',nam," New/",newn)
#    cat(phrase)
#    cat("\n")
#}
#sink(NULL)
#



#####  Scripts to run CoVaMa ####
files<-list.files("Output/PID_Consensus/", pattern = ".fasta")

temp<-readLines("Data/template/RunCoVaMa.sh")
for (i in 1:length(files)){
    fname<-files[i]
    fname<-gsub(".consensus.fasta",'',fname)
    new<-gsub(pattern="Run0_17", replace=paste0(fname),x=temp)
    writeLines(new, con=paste0("Data/Bash_CoVaMa/",fname,".covama.sh"))
}

