library(ggplot2)
library(gridExtra)
library(colorspace)
source("Rscripts/baseRscript.R")
cols<-qualitative_hcl(6, palette="Dark3")

# read the files saved in Overview_output:
SIVFiles<-list.files("Output/OverviewF_PID/",pattern=".csv")
Overview<-list()
for (i in 1:length(SIVFiles)){ 
        overviews<-read.csv(paste0("Output/OverviewF_PID/",SIVFiles[i]),stringsAsFactors=F, row.names = 1)
        Overview[[i]]<-overviews
        names(Overview)[i]<-substr(paste(SIVFiles[i]),start=1,stop=7)
}


gap<-data.frame(Sample=names(Overview))
for (i in 1:length(Overview)){
    df<-Overview[[i]]
    na<-df[is.na(df$freq.Ts),]
    poss<-na$pos
    if (length(poss)!=0) {
        gap$start[i]<- min(poss)
        gap$end[i]<- max(poss)
    }
    if (length(poss)==0) {
        gap$start[i]<- NA
        gap$end[i]<- NA
    }
}
    
write.csv(gap,"Output/gap_info.csv")


#Plasma only
samples<-read.csv("Data/SamplesNoduplicates.csv",stringsAsFactors = F)
samp<-samples[samples$Tissue=="Plasma",] #57 files

gap.plasma<-gap[gap$Sample %in% samp$File.name,]
min(gap.plasma$start, na.rm=T) #458
max(gap.plasma$end, na.rm=T) #506

#all files including tissues
min(gap$start, na.rm=T) #430
max(gap$end, na.rm=T) #600

ggplot(data=gap.plasma)+
    geom_segment(aes(x=Sample, xend=Sample, y=start, yend=end), size = 2, colour = "steelblue", alpha = 0.6) +
    theme_bw()+ylim(215,681)+
    theme(axis.text.x = element_text(angle=90))+
    ylab("Gap in env")+
    xlab('')
ggsave("Output/MF_PID/GapsinPlasma.pdf", width = 8, height=4)
    

gap$tissue<-"Plasma"
gap$tissue[!(gap$Sample %in% samp$File.name)]<-"Tissue"
ggplot(data=gap)+
    geom_segment(aes(x=Sample, xend=Sample, y=start, yend=end, color=tissue), size = 2, alpha = 0.6) +
    theme_bw()+ylim(215,681)+
    scale_color_manual(values=c("steelblue", "pink"))+
    theme(axis.text.x = element_text(angle=90, size=7), legend.title = element_blank())+
    ylab("Gap in env")+
    xlab('')
ggsave("Output/MF_PID/GapsinPlasma.pdf", width = 11, height=4)



#Remove 458 to 506 for all subsequent analysis for plasma analysis

for (i in samp$File.name){ 
    df<-Overview[[i]]
    df[df$pos>=458&df$pos<=506,7:16 ]<-NA
    write.csv(df, paste0("Output/OverviewF_PID_plasma/", i,".filter.nogap.overview.csv"))
}

#Apply the same filter to stock file

df<-read.csv("Output/OverviewF_PID/Run0_17_filtered.overview.csv",stringsAsFactors=F, row.names = 1)
df[df$pos>=458&df$pos<=506,7:16 ]<-NA
write.csv(df, paste0("Output/OverviewF_PID_plasma/Run0_17_.filter.nogap.overview.csv"))
