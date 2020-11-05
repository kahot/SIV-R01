#Script to analyse the frequency data and associate with features
library(dplyr)

source("Rscripts/BaseRscript2.R")

#get the file name
SIVFiles_SeqData<-list.files("Output/SeqData_PID/",pattern="SeqData")
#SIVFiles_SeqData<-SIVFiles_SeqData[s2]

#vec=c(29:31,33:36)
Overview<-list()
for (i in 1:length(SIVFiles_SeqData)){   
#for (i in vec){   
        id<-substr(paste(SIVFiles_SeqData[i]),start=9,stop=15)
        print(id)
        OverviewDF<-read.csv(paste0("Output/SeqData_PID/",SIVFiles_SeqData[i]),row.names = 1, stringsAsFactors=FALSE)
                
        TypeOfSite<-c() 
        TypeOfSite.tv1<-c()
        TypeOfSite.tv2<-c()
        TypeOfSite2<-c()
        TypeOfSite2.tv1<-c()
        TypeOfSite2.tv2<-c()
        
        for (codon in 1:(nrow(OverviewDF)/3)) { #modify based on reading frame
                positions <- c(codon*3-2,codon*3-1, codon*3)  #starting with codon 1
                WTcodon <- OverviewDF$MajNt[positions]  
                Refcodon<-OverviewDF$ref[positions] 
                if (is.na(WTcodon[1])|is.na(WTcodon[2])|is.na(WTcodon[3])){ 
                        WTcodon<-c('n','n','n')
                        mutant1codon<-c('n','n','n')
                        mutant2codon<-c('n','n','n')
                        mutant3codon<-c('n','n','n')
                        
                        mutant1codon.tv1 <- c('n','n','n')
                        mutant2codon.tv1 <- c('n','n','n')
                        mutant3codon.tv1 <-c('n','n','n')
                        
                        mutant1codon.tv2 <- c('n','n','n')
                        mutant2codon.tv2 <- c('n','n','n')
                        mutant3codon.tv2 <- c('n','n','n')
                    }
                    
                else {                        
                        mutant1codon <- c(transition(WTcodon[1]), WTcodon[2:3])  
                        mutant2codon <- c(WTcodon[1],transition(WTcodon[2]), WTcodon[3])
                        mutant3codon <- c(WTcodon[1:2], transition(WTcodon[3]))
                        
                        #transversion mutation to 'a' or 'c'
                        mutant1codon.tv1 <- c(transv1(WTcodon[1]), WTcodon[2:3]) 
                        mutant2codon.tv1 <- c(WTcodon[1],transv1(WTcodon[2]), WTcodon[3])
                        mutant3codon.tv1 <- c(WTcodon[1:2], transv1(WTcodon[3]))
                        #transversion mutation to 'g' or 't'
                        mutant1codon.tv2 <- c(transv2(WTcodon[1]), WTcodon[2:3])  
                        mutant2codon.tv2 <- c(WTcodon[1],transv2(WTcodon[2]), WTcodon[3])
                        mutant3codon.tv2 <- c(WTcodon[1:2], transv2(WTcodon[3]))
                }
                #compare to the ref seq
                mutant1codon2 <- c(transition(Refcodon[1]), Refcodon[2:3])  #If the first position has transistion mutation, it's labeld as mutatnt1codon.
                mutant2codon2 <- c(Refcodon[1],transition(Refcodon[2]), Refcodon[3])
                mutant3codon2 <- c(Refcodon[1:2], transition(Refcodon[3]))
                
                #transversion mutation to 'a' or 'c'
                mutant1codon2.tv1 <- c(transv1(Refcodon[1]), Refcodon[2:3]) 
                mutant2codon2.tv1 <- c(Refcodon[1],transv1(Refcodon[2]), Refcodon[3])
                mutant3codon2.tv1 <- c(Refcodon[1:2], transv1(Refcodon[3]))
                #transversion mutation to 'g' or 't'
                mutant1codon2.tv2 <- c(transv2(Refcodon[1]), Refcodon[2:3])  
                mutant2codon2.tv2 <- c(Refcodon[1],transv2(Refcodon[2]), Refcodon[3])
                mutant3codon2.tv2 <- c(Refcodon[1:2], transv2(Refcodon[3]))
                

                
                TypeOfSite<-c(TypeOfSite,typeofsitefunction(WTcodon,mutant1codon))
                TypeOfSite<-c(TypeOfSite,typeofsitefunction(WTcodon,mutant2codon))
                TypeOfSite<-c(TypeOfSite,typeofsitefunction(WTcodon,mutant3codon))
                
                TypeOfSite.tv1<-c(TypeOfSite.tv1,typeofsitefunction(WTcodon,mutant1codon.tv1))
                TypeOfSite.tv1<-c(TypeOfSite.tv1,typeofsitefunction(WTcodon,mutant2codon.tv1))
                TypeOfSite.tv1<-c(TypeOfSite.tv1,typeofsitefunction(WTcodon,mutant3codon.tv1))

                TypeOfSite.tv2<-c(TypeOfSite.tv2,typeofsitefunction(WTcodon,mutant1codon.tv2))
                TypeOfSite.tv2<-c(TypeOfSite.tv2,typeofsitefunction(WTcodon,mutant2codon.tv2))
                TypeOfSite.tv2<-c(TypeOfSite.tv2,typeofsitefunction(WTcodon,mutant3codon.tv2))
                
                
                TypeOfSite2<-c(TypeOfSite2,typeofsitefunction(Refcodon,mutant1codon2))
                TypeOfSite2<-c(TypeOfSite2,typeofsitefunction(Refcodon,mutant2codon2))
                TypeOfSite2<-c(TypeOfSite2,typeofsitefunction(Refcodon,mutant3codon2))
                
                TypeOfSite2.tv1<-c(TypeOfSite2.tv1,typeofsitefunction(Refcodon,mutant1codon2.tv1))
                TypeOfSite2.tv1<-c(TypeOfSite2.tv1,typeofsitefunction(Refcodon,mutant2codon2.tv1))
                TypeOfSite2.tv1<-c(TypeOfSite2.tv1,typeofsitefunction(Refcodon,mutant3codon2.tv1))
                
                TypeOfSite2.tv2<-c(TypeOfSite2.tv2,typeofsitefunction(Refcodon,mutant1codon2.tv2))
                TypeOfSite2.tv2<-c(TypeOfSite2.tv2,typeofsitefunction(Refcodon,mutant2codon2.tv2))
                TypeOfSite2.tv2<-c(TypeOfSite2.tv2,typeofsitefunction(Refcodon,mutant3codon2.tv2))     
        } 

        OverviewDF$Type<-TypeOfSite[1:length(OverviewDF$pos)]
        OverviewDF$Type.tv1<-TypeOfSite.tv1[1:length(OverviewDF$pos)]
        OverviewDF$Type.tv2<-TypeOfSite.tv2[1:length(OverviewDF$pos)]
        
        OverviewDF$Type.r<-TypeOfSite2[1:length(OverviewDF$pos)]
        OverviewDF$Type.tv1.r<-TypeOfSite2.tv1[1:length(OverviewDF$pos)]
        OverviewDF$Type.tv2.r<-TypeOfSite2.tv2[1:length(OverviewDF$pos)]
        
        Overview[[i]]<-OverviewDF[,-c(1:7)]
        
        names(Overview)[i]<-id   
}



###############################
#mutrates1<-read.csv("Data/HIVMutRates.csv")

Overview_sum<-list()
#for (i in vec){
for (i in 1:length(Overview)){

        OverviewDF<-Overview[[i]]
        id<-substr(paste(SIVFiles_SeqData[i]),start=9,stop=15)
        
        OverviewDF$WTAA<-""
        OverviewDF$MUTAA<-""
        OverviewDF$TVS1_AA<-""
        OverviewDF$TVS2_AA<-""
        OverviewDF$TVS2_AA<-""
    
        
        for (k in 1:nrow(OverviewDF)){
                   
                if (k%%3==1){
                        if (is.na(OverviewDF$MajNt[k])|is.na(OverviewDF$MajNt[k+1])|is.na(OverviewDF$MajNt[k+2])) { OverviewDF$WTAA[k]<-"NA"
                                OverviewDF$MUTAA[k]<-"NA"
                                OverviewDF$TVS1_AA[k]<-"NA"
                                OverviewDF$TVS2_AA[k]<-"NA"}
                        else {  OverviewDF$WTAA[k] = seqinr::translate(OverviewDF$MajNt[c(k,k+1,k+2)])
                                OverviewDF$MUTAA[k] = seqinr::translate(c(transition(OverviewDF$MajNt[k]),OverviewDF$MajNt[c(k+1,k+2)]))
                                OverviewDF$TVS1_AA[k] = seqinr::translate(c(transv1(OverviewDF$MajNt[k]),OverviewDF$MajNt[c(k+1,k+2)]))
                                OverviewDF$TVS2_AA[k] = seqinr::translate(c(transv2(OverviewDF$MajNt[k]),OverviewDF$MajNt[c(k+1,k+2)]))}
                        } 
                if (k%%3==2){
                        if (is.na(OverviewDF$MajNt[k-1])|is.na(OverviewDF$MajNt[k])|is.na(OverviewDF$MajNt[k+1]))  {OverviewDF$WTAA[k]<-"NA"
                        OverviewDF$MUTAA[k]<-"NA"
                        OverviewDF$TVS1_AA[k]<-"NA"
                        OverviewDF$TVS2_AA[k]<-"NA"}
                        else {  OverviewDF$WTAA[k] = seqinr::translate(OverviewDF$MajNt[c(k-1,k,k+1)])
                                OverviewDF$MUTAA[k] = seqinr::translate(c(OverviewDF$MajNt[c(k-1)],transition(OverviewDF$MajNt[k]),OverviewDF$MajNt[c(k+1)]))
                                OverviewDF$TVS1_AA[k] = seqinr::translate(c(OverviewDF$MajNt[c(k-1)],transv1(OverviewDF$MajNt[k]),OverviewDF$MajNt[c(k+1)]))
                                OverviewDF$TVS2_AA[k] = seqinr::translate(c(OverviewDF$MajNt[c(k-1)],transv2(OverviewDF$MajNt[k]),OverviewDF$MajNt[c(k+1)]))}
                        }
                if (k%%3==0){
                        if (is.na(OverviewDF$MajNt[k-2])|is.na(OverviewDF$MajNt[k-1])|is.na(OverviewDF$MajNt[k]))  {  OverviewDF$WTAA[k]<-"NA"
                        OverviewDF$MUTAA[k]<-"NA"
                        OverviewDF$TVS1_AA[k]<-"NA"
                        OverviewDF$TVS2_AA[k]<-"NA"}
                        else {  OverviewDF$WTAA[k] = seqinr::translate(OverviewDF$MajNt[c(k-2,k-1,k)])
                                OverviewDF$MUTAA[k] = seqinr::translate(c(OverviewDF$MajNt[c(k-2,k-1)],transition(OverviewDF$MajNt[k])))
                                OverviewDF$TVS1_AA[k] = seqinr::translate(c(OverviewDF$MajNt[c(k-2,k-1)],transv1(OverviewDF$MajNt[k])))
                                OverviewDF$TVS2_AA[k] = seqinr::translate(c(OverviewDF$MajNt[c(k-2,k-1)],transv2(OverviewDF$MajNt[k])))}
                }

        }
       #Add whether AA change is drastic & makes CpG
        OverviewDF$bigAAChange<-0
        OverviewDF$bigAAChange.tv1<-0
        OverviewDF$bigAAChange.tv2<-0
        OverviewDF$makesCpG <- 0
        OverviewDF$makesCpG.tvs <- 0
        OverviewDF$makesCpG.tv1 <- 0
        OverviewDF$makesCpG.tv2 <- 0
        
        for(j in 1:nrow(OverviewDF)){
                WT <- amCat(OverviewDF[j,'WTAA'])
                MUT <- amCat(OverviewDF[j,'MUTAA'])
                MUT1<-amCat(OverviewDF[j,'TVS1_AA'])
                MUT2<-amCat(OverviewDF[j,'TVS2_AA'])
                
                if (WT != MUT) OverviewDF$bigAAChange[j] <- 1
                if (WT != MUT1) OverviewDF$bigAAChange.tv1[j] <- 1
                if (WT != MUT2) OverviewDF$bigAAChange.tv2[j] <- 1
        
                trip <- OverviewDF$MajNt[c(j, j+1,j+2)]
                if (is.na(trip[1])|is.na(trip[2])|is.na(trip[3])) 
                        next
                        else{
                                if (trip[1] == "c" & trip[2] == "a" ) OverviewDF$makesCpG[j] <- 1 
                                if (trip[2] == "t" & trip[3] == "g")  OverviewDF$makesCpG[j] <- 1
                                if (trip[1] == "c" & (trip[2]=="c"|trip[2]=='t')) OverviewDF$makesCpG.tvs[j] <- 1
                                if (trip[3] == "g" & (trip[2]=="a"|trip[2]=="g")) OverviewDF$makesCpG.tvs[j] <- 1

                                if (trip[1] == "c" & (trip[2]=="c"|trip[2]=='t')) OverviewDF$makesCpG.tv2[j] <- 1                                
                                if (trip[3] == "g" & (trip[2]=="a"|trip[2]=="g")) OverviewDF$makesCpG.tv1[j] <- 1
                                
                        }
        } 
        
      
        write.csv(OverviewDF,paste0("Output/Overview_PID/",id,"_overview.csv"))
        
        #filter the sites with reads<10000
        remove<-which(OverviewDF$TotalReads<1000)
        OverviewDF[remove, 7:16]<-NA
        write.csv(OverviewDF,paste0("Output/OverviewF_PID/",id,"_filtered.overview.csv"))
}


#Filter at 5000 for Run4
SIVFiles_overview<-list.files("Output/Overview/",pattern="Run4.*_overview.csv")

for (i in 1:length(SIVFiles_overview)){ 
    overview<-read.csv(paste0("Output/Overview/",SIVFiles_overview[i]),stringsAsFactors=F, row.names = 1)
    remove<-which(overview$TotalReads<5000)
    overview[remove, 7:16]<-NA
    id<-substr(paste(SIVFiles_overview[i]),start=1,stop=7)
    write.csv(OverviewDF,paste0("Output/OverviewF/",id,"_filtered.overview.csv"))
}





#####################################
#Add the AA for ref
OverviewDF<-read.csv("Output/Overview/Run3_9__overview.csv",row.names = 1, stringsAsFactors=FALSE)
    
    OverviewDF$RefAA<-""
    OverviewDF$MUTAA.r<-""
    OverviewDF$TVS1_AA.r<-""
    OverviewDF$TVS2_AA.r<-""

    for (k in 1:nrow(OverviewDF)){
        if (k%%3==1){
            OverviewDF$RefAA[k] = seqinr::translate(OverviewDF$ref[c(k,k+1,k+2)])
            OverviewDF$MUTAA.r[k] = seqinr::translate(c(transition(OverviewDF$ref[k]),OverviewDF$ref[c(k+1,k+2)]))
            OverviewDF$TVS1_AA.r[k] =  seqinr::translate(c(transv1(OverviewDF$ref[k]),OverviewDF$ref[c(k+1,k+2)]))
            OverviewDF$TVS2_AA.r[k] =  seqinr::translate(c(transv2(OverviewDF$ref[k]),OverviewDF$ref[c(k+1,k+2)]))
        } 
        if (k%%3==2){
            OverviewDF$RefAA[k] =       seqinr::translate(OverviewDF$ref[c(k-1,k,k+1)])
            OverviewDF$MUTAA.r[k] =   seqinr::translate(c(OverviewDF$ref[c(k-1)],transition(OverviewDF$ref[k]),OverviewDF$ref[c(k+1)]))
            OverviewDF$TVS1_AA.r[k] = seqinr::translate(c(OverviewDF$ref[c(k-1)],transv1(OverviewDF$ref[k]),OverviewDF$ref[c(k+1)]))
            OverviewDF$TVS2_AA.r[k] = seqinr::translate(c(OverviewDF$ref[c(k-1)],transv2(OverviewDF$ref[k]),OverviewDF$ref[c(k+1)]))
        }
        if (k%%3==0){
            OverviewDF$RefAA[k] =     seqinr::translate(OverviewDF$ref[c(k-2,k-1,k)])
            OverviewDF$MUTAA.r[k] = seqinr::translate(c(OverviewDF$ref[c(k-2,k-1)],transition(OverviewDF$ref[k])))
            OverviewDF$TVS1_AA.r[k] = seqinr::translate(c(OverviewDF$ref[c(k-2,k-1)], transv1(OverviewDF$ref[k])))
            OverviewDF$TVS2_AA.r[k] = seqinr::translate(c(OverviewDF$ref[c(k-2,k-1)], transv2(OverviewDF$ref[k])))
        }
        
    }
    
    OverviewDF$bigAAChange.r<-0
    OverviewDF$bigAAChange.tv1.r<-0
    OverviewDF$bigAAChange.tv2.r<-0
    OverviewDF$makesCpG.r <- 0
    OverviewDF$makesCpG.tv1.r<- 0
    OverviewDF$makesCpG.tv2.r<- 0
    
    for(j in 1:nrow(OverviewDF)){
        WT <- amCat(OverviewDF[j,'RefAA'])
        MUT <- amCat(OverviewDF[j,'MUTAA.r'])
        MUT1<-amCat(OverviewDF[j,'TVS1_AA.r'])
        MUT2<-amCat(OverviewDF[j,'TVS2_AA.r'])
        
        if (WT != MUT) OverviewDF$bigAAChange.r[j] <- 1
        if (WT != MUT1) OverviewDF$bigAAChange.tv1.r[j] <- 1
        if (WT != MUT2) OverviewDF$bigAAChange.tv2.r[j] <- 1
        
        trip <- OverviewDF$ref[c(j, j+1,j+2)]
        if (is.na(trip[1])|is.na(trip[2])|is.na(trip[3])) 
            next
        else{
            if (trip[1] == "c" & trip[2] == "a" ) OverviewDF$makesCpG.r[j] <- 1 
            if (trip[2] == "t" & trip[3] == "g")  OverviewDF$makesCpG.r[j] <- 1
            if (trip[1] == "c" & (trip[2]=="c"|trip[2]=='t')) OverviewDF$makesCpG.tvs.r[j] <- 1
            if (trip[3] == "g" & (trip[2]=="a"|trip[2]=="g")) OverviewDF$makesCpG.tvs.r[j] <- 1
            
            if (trip[1] == "c" & (trip[2]=="c"|trip[2]=='t')) OverviewDF$makesCpG.tv2.r[j] <- 1                                
            if (trip[3] == "g" & (trip[2]=="a"|trip[2]=="g")) OverviewDF$makesCpG.tv1.r[j] <- 1
            
        }
    } 
    
OverviewDF<-OverviewDF[,c(1,4,20,21,22,34:44)]
write.csv(OverviewDF,"Output/Overview.ref.csv")
  

###############################       
# mean read depth
overviews<-list.files("Output/Overview_PID/",pattern="overview.csv")

for (i in 1:length(overviews)){ 
    df<-read.csv(paste0("Output/Overview_PID/",overviews[i]), row.names = 1)
    m<-mean(df$TotalReads, na.rm=T)
    print(substr(paste(overviews[i]),start=1,stop=7))
    print(m)
    #Ove[[i]]<-df
    #names(Ove)[i]<-substr(paste(overviews[i]),start=1,stop=7)
}

