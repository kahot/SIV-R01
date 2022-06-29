#### Compare PID processed (non-trimmed) vs trimmed at Q15 & PID-processed files for the 6 mutations
source("Rscripts/baseRscript.R")

### Sites of interest
diffpos<-c(475,494,500,515,518,578)
aapos<-ceiling(diffpos/3)
summary<-read.csv("Output/MF_PID/HighFreq/Summary_highFreq_sites.csv", row.names = 1)
summary2<-summary[summary$pos %in% diffpos,]
for (i in 1:nrow(summary2)){
    if (summary2$Type.1[i]=="Ts") summary2$mut[i]<-transition(summary2$ref[i])
    if (summary2$Type.1[i]=="Tv1") summary2$mut[i]<-transv1 (summary2$ref[i])
    if (summary2$Type.1[i]=="Tv2") summary2$mut[i]<-transv2 (summary2$ref[i])
}
summary2$mut[summary2$Type.1=="Ts"]<-transition(summary2$ref)
MutNames<-paste0(summary2$ref, summary2$pos,summary2$mut)


#Newly processed files (Run6)
files<-list.files("Output/Overview_PID/", pattern=".csv")
DF<-list()
#Trim reads # less than 10
for (i in 1:length(files)){
    df<-read.csv(paste0("Output/Overview_PID/", files[i]),stringsAsFactors=F, row.names = 1)
    remove<-which(df$TotalReads<10)
    df[remove, 7:16]<-NA
    DF[[i]]<-df[df$pos %in% diffpos,]
    names(DF)[i]<-substr(files[i], 1,7)
}

#Find the read depth and freq at the six sites from Q15 raw data
MF<-data.frame()
Depth<-data.frame()
for (i in 1:nrow(summary2)){
    muttype<-summary2$Type.1[i]
    if (muttype=="Tv1") muttype="transv1"
    if (muttype=="Tv2") muttype="transv2"
    
    df<-data.frame(File.name=names(DF))
    df[,"MF"]<-sapply(DF, function(x) x[i,paste0("freq.",muttype,".ref")])
    df$Mutation<-MutNames[i]
    MF<-rbind(MF, df)
    df2<-data.frame(File.name=names(DF))
    df2[,"Depth"]<-sapply(DF, function(x) x[i,"TotalReads"])
    df2$Mutation<-MutNames[i]
    Depth<-rbind(Depth, df2)

}

MF$Method<-"Trimmed"
Depth$Method<-"Trimmed"


####################################
### extract mf and depth for the six sites from PID processed data
SIVfiles<-list.files("Output/Overview_PIDcon/", pattern="Run6")

Ov<-list()
for (i in 1:length(SIVfiles)){ 
    df<-read.csv(paste0("Output/Overview_PIDcon/",SIVfiles[i]),stringsAsFactors=F, row.names = 1)
    remove<-which(df$TotalReads<10)
    df[remove, 7:16]<-NA
    Ov[[i]]<-df[df$pos %in% diffpos,]
    names(Ov)[i]<-substr(paste(SIVfiles[i]),start=1,stop=7)
}

#Find the read depth and freq at the six sites from Q15 raw data

PIDmf<-data.frame()
PIDdepth<-data.frame()
for (i in 1:nrow(summary2)){
    muttype<-summary2$Type.1[i]
    if (muttype=="Tv1") muttype="transv1"
    if (muttype=="Tv2") muttype="transv2"
    
    df<-data.frame(File.name=names(Ov))
    df[,"MF"]<-sapply(Ov, function(x) x[i,paste0("freq.",muttype,".ref")])
    df$Mutation<-MutNames[i]
    PIDmf<-rbind(PIDmf, df)
    df2<-data.frame(File.name=names(Ov))
    df2[,"Depth"]<-sapply(Ov, function(x) x[i,"TotalReads"])
    df2$Mutation<-MutNames[i]
    PIDdepth<-rbind(PIDdepth, df2)
    
}
PIDmf$Method<-"Non-trimmed"
PIDdepth$Method<-"Non-trimmed"


MF<-rbind(MF,PIDmf)
Depth<-rbind(Depth,PIDdepth)

MF$Mutation<-factor(MF$Mutation, levels=MutNames)
Depth$Mutation<-factor(Depth$Mutation, levels=MutNames)

ggplot()+
    ylab("Mutation frequency")+xlab("Sample")+
    facet_wrap(~ Mutation, ncol=3)+
    geom_point(data=MF, aes(x=File.name, y=MF, color=Method),alpha=0.8,position=position_dodge(width = 0.5),size=2)+
    theme_bw()+
    scale_color_manual(values=paste0(cols[c(1,4)],"CC"))+
    theme(axis.text.x = element_text(angle=90),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("Output/QC2/MF_comparison_Run6_trimmed.vs.non-trimmedPID.pdf", width = 10, height = 6 )

ggplot()+
    ylab("log(Read depth)")+xlab("Sample")+
    facet_wrap(~ Mutation, ncol=3)+
    geom_bar(data=Depth, aes(x=File.name,y=log10(Depth), fill=Method), alpha=0.7, position=position_dodge(width = 0.5), stat="identity")+
    theme_bw()+
    scale_fill_manual(values=paste0(cols[c(1,4)]))+
    theme(axis.text.x = element_text(angle=90),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("Output/QC2/Depth_comparison_Run6.pdf", width = 10, height = 6 )


write.csv(MF, "Output/QC2/MF.compare_6sites.csv")
write.csv(Depth, "Output/QC2/Depth.compare_6sites.csv")


MF2<-cbind(MF, Depth[,2])
colnames(MF2)[5]<-"Depth"
ggplot(MF2, aes(x=File.name))+
    facet_wrap(~ Mutation, ncol=6)+
    geom_point(aes(y=MF*4, color=Method),alpha=0.8,position=position_dodge(width = 0.5),size=2)+
    geom_bar(aes(y=log10(Depth), fill=Method), alpha=0.3, position=position_dodge(width = 0.5), stat="identity")+
    geom_point(aes(y=MF*4, color=Method),alpha=0.8,position=position_dodge(width = 0.5),size=2)+
    scale_y_continuous(
            name="Mutation frequency", breaks = c(0,2,4),labels = c(0,0.5,1),
            sec.axis=sec_axis(~., name="log(Read depth)"))+
    xlab("Sample")+
    theme_bw()+
    scale_color_manual(values=cols[c(1,4)])+
    scale_fill_manual(values=cols[c(1,4)])+
    theme(axis.text.x = element_text(angle=90),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("Output/QC2/MF_Depth_comparison_Run6.pdf", width = 25, height = 6 )


write.csv(MF2,"Output/QC/MF_Depth_comparison_Run6.csv")

