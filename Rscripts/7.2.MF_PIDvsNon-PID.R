library(ggplot2)
library(reshape)
library(ggpubr)
library(ggthemes)
library(plotrix)
library(grid)
library(gridExtra)
library(colorspace)
library(cowplot)
library(DataCombine)
source("Rscripts/baseRscript2.R")
cols2<-qualitative_hcl(6, palette="Dark3")

#read the PID treated Overview files:
SIVFiles_overviewP<-list.files("Output/OverviewF_PID/",pattern="filtered.overview.csv")
OverviewP<-list()
for (i in 1:length(SIVFiles_overviewP)){ 
        overviews<-read.csv(paste0("Output/OverviewF_PID/",SIVFiles_overviewP[i]),stringsAsFactors=F, row.names = 1)
        OverviewP[[i]]<-overviews
        names(OverviewP)[i]<-substr(paste(SIVFiles_overviewP[i]),start=1,stop=7)
}


# Non-PID Overview_output:
SIVFiles_overview<-list.files("Output/OverviewF/",pattern="filtered.overview.csv")
Overview<-list()
for (i in 1:length(SIVFiles_overview)){ 
    overviews<-read.csv(paste0("Output/OverviewF/",SIVFiles_overview[i]),stringsAsFactors=F, row.names = 1)
    Overview[[i]]<-overviews
    names(Overview)[i]<-substr(paste(SIVFiles_overview[i]),start=1,stop=7)
}



samples<-read.csv("Data/SamplesNoduplicates.csv",stringsAsFactors = F)
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
monkeys2<-monkeyList[ids]

######################

#Average mutation frequencies per samples 
mf.summary<-data.frame(sample=names(OverviewP))

for (i in 1:nrow(mf.summary)){
    fname<-names(OverviewP)[i]
    df1<-Overview[[fname]]
    dfP<-OverviewP[[i]]
    
    mf.summary$MF[i]<-mean(df1$freq.mutations.ref, na.rm = T)
    mf.summary$MF_PID[i]<-mean(dfP$freq.mutations.ref, na.rm = T)
    
    #For mutations that are over 0.1 only
    
    mf.summary$MFover0.1[i]<-mean(df1$freq.mutations.ref[df1$freq.mutations.ref>=0.1], na.rm = T)
    mf.summary$MFover0.1_PID[i]<-mean(dfP$freq.mutations.ref[dfP$freq.mutations.ref>=0.1], na.rm = T)
    
    #read depth
    mf.summary$Depth.ave[i]<-mean(df1$TotalReads[!is.na(df1$freq.Ts.ref)], na.rm = T)
    mf.summary$Depth.ave_PID[i]<-mean(dfP$TotalReads[!is.na(dfP$freq.Ts.ref)], na.rm = T)
    
    #max and min depth
    mf.summary$Depth.min[i]<-min(df1$TotalReads[!is.na(df1$freq.Ts.ref)], na.rm = T)
    mf.summary$Depth.min_PID[i]<-min(dfP$TotalReads[!is.na(dfP$freq.Ts.ref)], na.rm = T)
    
    mf.summary$Depth.max[i]<-max(df1$TotalReads[!is.na(df1$freq.Ts.ref)], na.rm = T)
    mf.summary$Depth.max_PID[i]<-max(dfP$TotalReads[!is.na(dfP$freq.Ts.ref)], na.rm = T)
    
    #how many site were over >=0.1 in mut freq?
    df1.2<-df1[df1$freq.Ts.ref>=0.1|df1$freq.transv1.ref>=0.1|df1$freq.transv2.ref>=0.1,]
    df1.2<-df1.2[!is.na(df1.2$pos),]
    dfP.2<-dfP[dfP$freq.Ts.ref>=0.1|dfP$freq.transv1.ref>=0.1|dfP$freq.transv2.ref>=0.1,]
    dfP.2<-dfP.2[!is.na(dfP.2$pos),]
    
    mf.summary$no.of.sites[i]<-nrow(df1.2)
    mf.summary$no.of.sites_PID[i]<-nrow(dfP.2)
    
 }

write.csv(mf.summary, "Output/MF_comparison.csv")

run0.14<-read.csv("Output/Overview_PID/Run0_14_overview.csv")
which(mf.summary$sample=="Run0_14")
dfP<-run0.14
mf.summary$MF_PID[4]<-mean(dfP$freq.mutations.ref, na.rm = T)
mf.summary$MFover0.1_PID[4]<-mean(dfP$freq.mutations.ref[dfP$freq.mutations.ref>=0.1], na.rm = T)
mf.summary$Depth.ave_PID[4]<-mean(dfP$TotalReads[!is.na(dfP$freq.Ts.ref)], na.rm = T)
mf.summary$Depth.min_PID[4]<-min(dfP$TotalReads[!is.na(dfP$freq.Ts.ref)], na.rm = T)
mf.summary$Depth.max_PID[4]<-max(dfP$TotalReads[!is.na(dfP$freq.Ts.ref)], na.rm = T)
dfP.2<-dfP[dfP$freq.Ts.ref>=0.1|dfP$freq.transv1.ref>=0.1|dfP$freq.transv2.ref>=0.1,]
dfP.2<-dfP.2[!is.na(dfP.2$pos),]
mf.summary$no.of.sites_PID[4]<-nrow(dfP.2)
write.csv(mf.summary, "Output/MF_comparison_run0_14replaced.csv")



## Plot the results

mf.summary$MF_difference<-mf.summary$MF_PID-mf.summary$MF
mf.summary$MF0.1_difference<-mf.summary$MFover0.1_PID-mf.summary$MFover0.1

ggplot(mf.summary, aes(x=sample, y=MF_difference))+
    geom_hline(yintercept=0, col="gray")+
    geom_point()+
    theme_bw()+ylab('MF difference (PID-non.PID)')+
    theme(axis.text.x = element_text(size=9, angle=90),
          axis.text.y=element_text(size=9))
ggsave("Output/MF_PID/MF_comparison.pdf", height = 4, width = 7)


ggplot(mf.summary, aes(x=sample, y=MF0.1_difference))+
    geom_hline(yintercept=0, col="gray", size=0.5)+
    geom_point()+
    theme_bw()+ylab('MF difference (PID-non.PID, >0.1)')+
    theme(axis.text.x = element_text(size=9, angle=90),axis.text.y=element_text(size=9))
ggsave("Output/MF_PID/MFover0.1_comparison.pdf", height = 4, width = 7)

sum2<-mf.summary[,c("sample","no.of.sites","no.of.sites_PID")]
colnames(sum2)[2:3]<-c("Non-PID","PID")
sum2<-melt(sum2)

ggplot(sum2, aes(x=sample, y=value, color=variable))+
    geom_hline(yintercept=0, col="gray", size=0.5)+
    scale_color_manual(values=cols2[c(4,6)])+
    geom_point()+
    theme_bw()+ylab('# of sites with m.f. >0.1')+
    theme(axis.text.x = element_text(size=9, angle=90))
ggsave("Output/MF_PID/No.of.sites.over.mf0.1.pdf", height = 4, width = 8)



#Find the over lapping positions
positions<-unique(unlist(unname(Pos)))
positions<-positions[order(positions)]
cpos<-data.frame(Pos=positions)
for (i in 1:length(monkeys2)){
    vec<-Pos[[i]]
    
    cpos[,names(monkeys2)[i]]<-apply(cpos['Pos'],2, function(x) ifelse(x %in% vec, "Y", "N"))
}

cpos$codon<-apply(cpos["Pos"], 1, function(x) if(x%%3==0) x=3 else x=x%%3)


df<-read.csv("Output/Overview.ref.csv", stringsAsFactors = F, row.names = 1)
df<-df[df$pos %in% positions,]
colnames(cpos)[1]<-"pos"
cpos<-merge(cpos, df, by="pos")

write.csv(cpos, "Output/HighMutfreq_sites_all_PID.csv")



