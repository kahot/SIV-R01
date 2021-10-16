library(ggplot2)
library(gridExtra)
library(colorspace)
library(zoo)
source("Rscripts/label_scientific.R")
#source("Rscripts/baseRscript2.R")
cols<-qualitative_hcl(6, palette="Dark3")

# read the files saved in Overview_output:
SIVFiles_overview<-list.files("Output/OverviewF_PID/",pattern="filtered.overview.csv")

#Use non-filtered files for Run5 & Run5
SIVfiles<-list.files("Output/Overview_PID/", pattern="Run5|Run6")
SIVfiles<-SIVfiles[SIVfiles!="Run5_12_overview.csv"]
SIVfiles<-SIVfiles[SIVfiles!="Run5_13_overview.csv"]
SIVfiles<-SIVfiles[SIVfiles!="Run5_15_overview.csv"]

Overview<-list()
for (i in 1:length(SIVFiles_overview)){ 
        overviews<-read.csv(paste0("Output/OverviewF_PID/",SIVFiles_overview[i]),stringsAsFactors=F, row.names = 1)
        Overview[[i]]<-overviews
        names(Overview)[i]<-substr(paste(SIVFiles_overview[i]),start=1,stop=7)
}
n<-length(SIVFiles_overview)

for (i in (n+1):(length(SIVfiles)+n)){
    df<-read.csv(paste0("Output/Overview_PID/", SIVfiles[i-n]),stringsAsFactors=F, row.names = 1)
    Overview[[i]]<-df
    names(Overview)[i]<-substr(paste(SIVfiles[i-n]),start=1,stop=7)
}


samples<-read.csv("Data/SamplesNoduplicates.csv",stringsAsFactors = F)
samples$Week<-as.integer(samples$Week)
samples$Tissue<-factor(samples$Tissue, levels=c("Stock","Plasma","LN","HLN","Unknown"))

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

#TB infection week
tbs<-read.csv("Data/TBinfection.week.csv",stringsAsFactors = F)

monkeys2<-monkeyList[tbs$ids]

############################
#diversity vs. viral load correlation analysis

VL<-read.csv("Data/ViralLoads.csv", row.names = 1, stringsAsFactors = F)

morder<-paste0("A",unique(VL$Monkey))

Plots<-list()
for (i in 1:length(morder)){
    sample<-monkeys2[[morder[i]]]
    sample<-sample[order(sample$Week),]
    #select plasma samples
    sample<-sample[sample$Tissue=="Plasma",]
    ovDF<-Overview[as.vector(sample$File.name)]
    monkey<-morder[i]
    tbweek<-tbs$tb[tbs$ids==monkey]
    diversity<-sample[,c(2,3,5,7)]
    monkey<-gsub("A",'',monkey)
    
    for (j in 1:length(ovDF)){
        df<-ovDF[[j]]
        diversity$Diversity[j]<-mean(df$freq.mutations, na.rm=T)*100
        diversity$MF[j]<-mean(df$freq.mutations.ref, na.rm=T)*100
    }
    
    vl<-VL[VL$Monkey==monkey,]
    diversity<-merge(diversity, vl[,1:2], by="Week" )
    
    re<-cor.test(diversity$Diversity, diversity$VL, method="spearman")
    pval<-round(re[[3]],3)
    Plots[[i]]<-
        ggplot(diversity,aes(x=VL, y=Diversity))+
            ylab("Diversity")+xlab("VL")+
            scale_x_continuous(trans="log")+
            geom_point(size=2)+theme_bw()+
            ggtitle(paste0(monkey,"  p = ",pval ))+
            theme(axis.text.x = element_text(angle=45,hjust=1))
}    

pdf("Output/MF_PID/Diversity.vs.VL.all.pdf", width = 10, height = 10)
do.call(grid.arrange, c(Plots, ncol=3))
dev.off()


# Calculate rolling average for VL pre- and post-Mtb infection

Plots<-list()
for (i in 1:length(morder)){
    monkeyID<-morder[i]
    monkey<-gsub("A",'',monkeyID)
    vl<-VL[VL$Monkey==monkey,]
    #remove the week 0 
    vl<-vl[vl$Week!=0,]
   
    #add all weeks from 1:32
    #Weeks<-data.frame(Week=seq(1:max(vl$Week)))
    #vl<-merge(Weeks, vl[,1:2], all = T)
    
    #calculate pre- and post-mtb trend
    if (i>7){
        r3<-rollmean(vl$VL, k=3, align = "center",na.rm=T) 
        vl$roll3<-c(NA,r3,NA)
        vl$mtb="pre"
    }
    else{
        mtb<-tbs$tb[tbs$ids==monkeyID]
        pre_r3<- rollmean(vl$VL[vl$Week<mtb], k=3, align = "center",na.rm=T) 
        post_r3<-rollmean(vl$VL[vl$Week>=mtb],k=3, align = "center",na.rm=T) 
        
        vl$roll3<-c(NA,pre_r3, rep(NA, times=2), post_r3,NA)
        vl$mtb<-c(rep("pre", times=length(pre_r3)+2), rep("post", times=length(post_r3)+2))
        vl$mtb<-factor(vl$mtb, levels = c("pre","post"))
    }
    
        vl$roll3[is.nan(vl$roll3)]<-NA
        vl2<-vl[!is.na(vl$roll3),]
    breaks=c(1E6,1E7,1E8,1E9)
    p<-ggplot(vl, aes(x=Week,y=VL))+
        geom_point(color="gray")+ylab("Viral load")+
        scale_y_continuous(trans="log", breaks=breaks, label=label_scientific)+
        geom_line(data=vl2,aes(x=Week, y=roll3, color=mtb))+
        ggtitle(monkey)+
        theme_bw()+
        theme(legend.position = "none")
    
    if (!(i==8|i==9)){
        p<-p+geom_vline(xintercept = mtb, color="pink", size=0.5)
    }
    Plots[[i]]<-p
}

pdf("Output/MF_PID/VL.overtime.all.pdf", width = 12, height = 10)
do.call(grid.arrange, c(Plots, ncol=3))
dev.off()



#Add background rolling average of 5
Plots2<-list()
for (i in 1:length(morder)){
    monkeyID<-morder[i]
    monkey<-gsub("A",'',monkeyID)
    vl<-VL[VL$Monkey==monkey,]
    #remove the week 0 
    vl<-vl[vl$Week!=0,]
    
    #original vl
    vlo<-vl
    roll<-rollmean(vlo$VL, k=5, align = "center",na.rm=T) 
    vlo$roll<-c(NA,NA,roll,NA,NA)
    

    
    #add all weeks from 1:32
    Weeks<-data.frame(Week=seq(1:max(vl$Week)))
    vl<-merge(Weeks, vl[,1:2], all = T)
    
    #calculate pre- and post-mtb trend
    if (i>7){
        #r3<-rollmean(vl$VL, k=3, align = "center",na.rm=T) 
        #vl$roll3<-c(NA,r3,NA)
      
        
        r<-rollmean(vl$VL, k=4, align = "center",na.rm=T) 
        vl$roll<-c(NA,NA,r,NA)
        vl$mtb="pre"
    }
    else{
        mtb<-tbs$tb[tbs$ids==monkeyID]
        pre_r3<- rollmean(vl$VL[vl$Week<mtb], k=4, align = "center",na.rm=T) 
        post_r3<-rollmean(vl$VL[vl$Week>=mtb],k=4, align = "center",na.rm=T) 
        
        vl$roll3<-c(NA,NA,pre_r3, rep(NA, times=2), post_r3,NA,NA)
        vl$mtb<-c(rep("pre", times=length(pre_r3)+3), rep("post", times=length(post_r3)+3))
        vl$mtb<-factor(vl$mtb, levels = c("pre","post"))
    }
    
    vl$roll[is.nan(vl$roll)]<-NA
    vl2<-vl[!is.na(vl$roll),]
    breaks=c(1E6,1E7,1E8,1E9)
    p<-ggplot(vl, aes(x=Week,y=VL))+
        geom_point(color="gray")+ylab("Viral load")+
        geom_line(data=vlo, aes(x=Week, y=roll), color="gray80", linetype=2)+
        scale_y_continuous(trans="log", breaks=breaks, label=label_scientific)+
        geom_line(data=vl2,aes(x=Week, y=roll, color=mtb))+
        ggtitle(monkey)+
        theme_bw()+
        theme(legend.position = "none")
    
    if (!(i==8|i==9)){
        p<-p+geom_vline(xintercept = mtb, color="pink", size=0.5)
    }
    Plots2[[i]]<-p
}

pdf("Output/MF_PID/VL.overtime.roll4.pdf", width = 12, height = 10)
do.call(grid.arrange, c(Plots2, ncol=3))
dev.off()










#calculate pre- and post-mtb trend
mtb<-tbs$tb[tbs$ids==paste0("A",monkey)]
pre_r3<-rollmean(vl$VL[vl$Week<mtb], k=3, align = "center") #7
post_r3<-rollmean(vl$VL[vl$Week>=mtb], k=3, align = "center") #8

nrow(vl[vl$Week<mtb,])
length(vl$VL[vl$Week>=mtb])

vl$roll3<-c(rep(NA, times=1),pre_r3, rep(NA, times=2), post_r3,NA)
vl$color<-c(rep("pre", times=7), rep("post", times=8))
ggplot(vl, aes(x=Week,y=VL))+
    geom_point(color="gray")+
    scale_y_continuous(trans="log", label=label_scientific())+
    geom_line(aes(x=Week, y=roll3, color=color))+
    ggtitle(monkey)+
    geom_vline(xintercept = mtb, color="pink", size=0.5)

    