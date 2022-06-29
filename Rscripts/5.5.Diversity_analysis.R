library(gridExtra)
library(ggplot2)
library(colorspace)
library(DataCombine)
library(reshape2)
library(dplyr)
library(tidyverse)
#source("Rscripts/baseRscript2.R")
cols8<-qualitative_hcl(8, palette="Dark3")


#Final plasma diversity comparison

#1. SIV-only vs. co-infected
####
tbs<-read.csv("Data/TBinfection.week.csv",stringsAsFactors = F)
samples<-read.csv("Data/SamplesNoduplicates.csv",stringsAsFactors = F)
samples$Week<-as.integer(samples$Week)

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

ID<-c("A21918","A22117","A22217","A22317","A22517","A22617","A23918","A34119","A34219")
monkeyList<-monkeyList[ID]
monkeys<-names(monkeyList)


#Calculate the gap bases for plasma: 
SIVFiles<-list.files("Output/OverviewF_PID/")
Ov<-list()
for (i in 1:length(SIVFiles)){ 
    overviews<-read.csv(paste0("Output/OverviewF_PID/",SIVFiles[i]),stringsAsFactors=F, row.names = 1)
    Ov[[i]]<-overviews
    names(Ov)[i]<-substr(paste(SIVFiles[i]),start=1,stop=7)
}
Ov2<-list()
gaps<-data.frame(ID=names(Ov))
for (i in 1:length(Ov)){
    df<-Ov[[i]]
    id<-names(Ov)[i]
    df<-df[df$pos>=215,]
    Ov2[[i]]<-df
    names(Ov2)[i]<-id
    df2<-df[is.na(df$freq.mutations.ref),]
    #write.csv(df,paste0("Output/OverviewF_PID/",id,"_filtered.overview.csv"))
    
    
    gaps$gap.size[i]<-nrow(df2)
    gaps$start[i]<-ifelse(nrow(df2)>0, df2$pos[1], NA)
    gaps$end[i]<-ifelse(nrow(df2)>0, df2$pos[nrow(df2)], NA)
}

write.csv(gaps,"Output/Diversity/GapInfo.csv")
min(gaps$start, na.rm=T) #423
max(gaps$end, na.rm=T) #610


Div<-data.frame(ID=names(Ov))
Ovcut<-list()
for (i in 1:length(Ov2)){
    df<-Ov2[[i]]
    id<-names(Ov2[i])
    #remove pos 423 to 610
    df<-df[df$pos<423|df$pos>610,]
    Ovcut[[i]]<-df
    names(Ovcut)[i]<-id
    Div$Mean[i]<-mean(df$freq.mutations, na.rm=T)
    
}
#dir.create("Output/Diversity/")
write.csv(Div, "Output/Diversity/Ave.diversity.all_remove423-610.csv")



EndDiv<-data.frame(Monkey=monkeys)
for (m in 1:length(monkeyList)){
    #Select the monkey
    sample<-monkeyList[[m]]
    sample<-sample[order(sample$Week),]
    #Plasma only
    sample<-sample[sample$Tissue=="Plasma",]
    fname<-sample$File.name[sample$Week==max(sample$Week)]
    
    EndDiv$ave.div[m]<-Div$Mean[Div$ID==fname]
}

EndDiv$Coinfection<-c(rep("Y", times=7), rep("N",times=2))
EndDiv$Monkey<-gsub("A",'',EndDiv$Monkey)

ggplot(EndDiv, aes(x=Monkey, y=ave.div*100, color=Coinfection))+theme_light()+
    geom_point(size=2)+ylab("Ave. plasma diversity (%) at necropsy")+xlab("")+
    theme(axis.text.x = element_text(angle=45, hjust=1))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    geom_vline(xintercept = 1:8+0.5, color="gray", size=0.3)+
    geom_vline(xintercept = 7.5, color="orange", size=0.5)
ggsave("Output/Diversity/Div.at.necropsy.pdf",width = 6, height=3.4)




div.time<-data.frame(Monkey=monkeys[1:7])
Div.prepost<-data.frame()
#Coinfection only
for (m in 1:7){
    #Select the monkey
    sample<-monkeyList[[m]]
    sample<-sample[order(sample$Week),]
    #Plasma only
    sample<-sample[sample$Tissue=="Plasma",]
    
    monkey<-names(monkeyList)[m]
    tbweek<-tbs$tb[tbs$ids==monkey]
    pre<-Div[Div$ID %in% sample$File.name[sample$Week<=tbweek],]
    post<-Div[Div$ID %in% sample$File.name[sample$Week>tbweek],]
    
    pre$Infection<-"Pre"
    post$Infection<-"Post"
    df<-rbind(pre,post)
    df$Monkey<-monkey
    Div.prepost<-rbind(Div.prepost, df)
    
    div.time$Pre[m]<-mean(pre$Mean, na.rm = T)
    div.time$Post[m]<-mean(post$Mean, na.rm = T)
    
    
    
}

divtime<-melt(div.time, id.vars="Monkey")
Div.prepost$Infection<-factor(Div.prepost$Infection, levels=c("Pre","Post"))
ggplot()+theme_linedraw()+xlab('')+ylab("Average diversity (%)")+
    geom_point(data=Div.prepost, aes(x=Monkey, y=Mean*100, color=Infection), 
               position=position_jitterdodge(jitter.width = 0.2,dodge.width = 0.7), alpha=0.5)+
    scale_color_manual(values=paste0(cols8[c(5,1)]), guide="none")+
    scale_fill_manual(values=cols8[c(5,1)])+
    geom_point(data=divtime, aes(x=Monkey, y=value*100, fill=variable),shape=21,size=2.7, 
               position=position_dodge(width = 0.7), color="gray20")+
    geom_vline(xintercept = 1:6+0.5, color="gray", size=0.5)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank())
ggsave("Output/Diversity/Diversity.pre-postMtbinfection.pdf",width = 7, height = 4)



#Diversity change over time per animal

TB<-tbs[tbs$ids %in% ID,]
TB$Monkey<-gsub("A", '',TB$ids)

Diversity<-data.frame()
#plots<-list()
for (i in 1:length(monkeys)){
    sample<-monkeyList[[i]]
    sample<-sample[order(sample$Week),]
    DF<-Ovcut[as.vector(sample$File.name)]
    monkey<-monkeys[i]
    tbweek<-tbs$tb[tbs$ids==monkey]
    
    meanDIV<-Div[Div$ID %in% sample$File.name,]
    
    diversity<-sample[,c(2,3,4,5,7)]
    monkey<-gsub("A",'',monkey)
    diversity<-merge(diversity, meanDIV, by.x="File.name", by.y="ID")
    
    diversity<-InsertRow(diversity,c("Stock",monkeys[i],0,"Stock",NA,Div$Mean[Div$ID=="Run0_17"]),1)
    diversity$Week<-as.integer(diversity$Week)
    diversity$Mean<-as.numeric(diversity$Mean)
    
    diversity$Tissue<-factor(diversity$Tissue, levels=c("Stock","Plasma","LN","HLN","Lung","Colon"))
    diversity$tb<-TB$tb[TB$Monkey==monkey]
    Diversity<-rbind(Diversity,diversity)
}

Diversity$Monkey<-gsub("A",'', Diversity$Monkey)
ggplot(Diversity)+
    facet_wrap(~Monkey,ncol=3)+
    ylab("Diversity")+xlab("")+
    geom_point(data=Diversity, aes(x=Week, y=Mean*100, color=Tissue), size=2)+theme_bw()+
    ylim(0.3,2.2)+
    scale_color_manual(values=cols8[c(3,2,7,5,1,8)])+
    theme(panel.grid.major.x  = element_blank(),panel.grid.minor.x = element_blank())+
    geom_vline(data=Diversity, aes(xintercept=tb), col="deeppink")+
    #scale_x_continuous(breaks=seq(0,33,2), limits = c(0,33))+
    geom_line(data=Diversity[Diversity$Tissue=="Plasma"|Diversity$Tissue=="Stock",],aes(x=Week, y=Mean*100), color=cols8[2])

ggsave("Output/Diversity/Diversity_overtime_9monkeys.pdf", width = 10, height = 7)



## Compare tissues vs. plasma diversity
#eliminate monkey '23918' (no tissue samples)

div.ave<-data.frame(Monkey=monkeys[-which(monkeys=="A23918")])
Div.tis<-data.frame()
monkeys3<-monkeys[-which(monkeys=="A23918")]

for (m in 1:length(monkeys3)){
    #Select the monkey
    monkey<-monkeys3[m]
    
    sample<-monkeyList[[monkey]]
    sample<-sample[order(sample$Week),]
    
    tbweek<-tbs$tb[tbs$ids==monkey]
    plasma<-Div[Div$ID %in% sample$File.name[sample$Tissue=="Plasma"],]
    tis<-Div[Div$ID %in% sample$File.name[sample$Tissue!="Plasma"],]
    
    plasma$Tissue<-"Plasma"
    tis$Tissue<-"Tissue"
    df<-rbind(plasma,tis)
    df$Monkey<-monkey
    Div.tis<-rbind(Div.tis, df)
    
    div.ave$Plasma[m]<-mean(plasma$Mean, na.rm = T)
    div.ave$Tissue[m]<-mean(tis$Mean, na.rm = T)
}


divave<-melt(div.ave, id.vars="Monkey")
divave$Monkey<-gsub("A",'',divave$Monkey)
Div.tis$Tissue<-factor(Div.tis$Tissue, levels=c("Plasma","Tissue"))
Div.tis$Monkey<-gsub("A",'',Div.tis$Monkey)

ggplot()+theme_linedraw()+xlab('')+ylab("Average diversity (%)")+
    geom_point(data=Div.tis, aes(x=Monkey, y=Mean*100, color=Tissue), 
               position=position_jitterdodge(jitter.width = 0.2,dodge.width = 0.7), alpha=0.5)+
    scale_color_manual(values=paste0(cols8[c(1,7)]), guide="none")+
    scale_fill_manual(values=cols8[c(1,7)])+
    geom_point(data=divave, aes(x=Monkey, y=value*100, fill=variable),shape=21,size=2.7, 
               position=position_dodge(width = 0.7), color="gray20")+
    geom_vline(xintercept = 1:7+0.5, color="gray", size=0.5)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank())
ggsave("Output/Diversity/Diversity.plasma.vs.tussues_all.pdf",width = 7, height = 4)



#### Compare the tissue and plasma diversity (%) at necropsy
Div<-read.csv( "Output/Diversity/Ave.diversity.all_remove423-610.csv", row.names = 1)

Nec<-data.frame()
Mid<-data.frame()
diversity<-data.frame(Animal=monkeys3)

for (i in 1:length(monkeys3)){
    #Select the monkey
    monkey<-monkeys3[i]
    sample<-monkeyList[[monkey]]
    sample<-sample[order(sample$Week),]
    
    sum<-merge(sample, Div, by.x="File.name", by.y="ID")
    #replace tissue discription as "Tissue"
    sum$Tissue[sum$Tissue!="Plasma"]<-"Tissue"
    #Midpoint tissue vs. plasma
    weeks<-sum$Week[sum$Tissue!="Plasma"]
    if (length(weeks)==0) {
        diversity$Mid_Plasma[i]<-NA
        diversity$Mid_Tissue[i]<-NA
        diversity$Nec_Plasma[i]<-NA
        diversity$Nec_Tissue[i]<-NA
    }
    if (length(weeks)!=0) {
        if (length(unique(weeks))>1) {midpoint<-min(unique(weeks))
        mids<-(midpoint-2):(midpoint+2)
        sum1<-sum[sum$Week %in% mids,]
        diversity$Mid_Plasma[i]<-sum1$Mean[sum1$Tissue=="Plasma"]
        diversity$Mid_Tissue[i]<-sum1$Mean[sum1$Tissue!="Plasma"]
        Mid<-rbind(Mid, sum1)
        }
        else {
            diversity$Mid_Plasma[i]<-NA
            diversity$Mid_Tissue[i]<-NA
        }
        #Necropsy tissues vs. plasma
        sum2<-sum[sum$Week==max(sum$Week),]
        Nec<-rbind(Nec, sum2)
        div<-aggregate(sum2["Mean"], by=list(sum2$Tissue), mean)
        diversity$Nec_Plasma[i]<-ifelse ("Plasma" %in% div$Group.1,  div$Mean[div$Group.1=="Plasma"], NA)
        diversity$Nec_Tissue[i]<-ifelse ("Tissue" %in% div$Group.1,  div$Mean[div$Group.1=="Tissue"], NA)
    }
}

Div2<-diversity[rowSums(is.na(diversity))!=4,] 

aveDiv<-melt(Div2)
aveDiv$Timepoint<-aveDiv$variable
aveDiv$Timepoint<-gsub("_Plasma|_Tissue",'',aveDiv$Timepoint)
aveDiv$Tissue<-aveDiv$variable
aveDiv$Tissue<-gsub("Mid_|Nec_",'',aveDiv$Tissue)
aveDiv<-aveDiv[,c(1,3,4,5)]
aveDiv$Timepoint[aveDiv$Timepoint=="Mid"]<-"Mid-point"
aveDiv$Timepoint[aveDiv$Timepoint=="Nec"]<-"Necropsy"

aveDiv$Timepoint<-factor(aveDiv$Timepoint, levels=c("Mid-point","Necropsy"))
aveDiv$Monkey<-gsub("A",'',aveDiv$Animal)

ggplot(data=aveDiv, aes(x=Monkey, y=value*100, fill=Tissue, shape=Timepoint))+xlab('')+ylab("Average diversity (%)")+
    geom_point(size=2.7, 
               position=position_dodge(width = 1), color="gray20")+
    geom_vline(xintercept = 1:7+0.5, color="gray", size=0.5)+
    scale_fill_manual(values=cols8[c(1,7)])+
    scale_shape_manual(values=c(21,23))+theme_linedraw()+
    guides(fill=guide_legend(override.aes=list(shape=21)))+
    geom_vline(xintercept = 6.5, color=cols8[7])+
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), legend.title = element_blank())
ggsave("Output/Diversity/Diversity.plasma.vs.tissues.pdf",width = 7, height = 4)

#connect lines
p<-ggplot(data=aveDiv, aes(x=Animal, y=value*100, color=Tissue, fill=Tissue, shape=Timepoint))+
    xlab('')+ylab("Average diversity (%)")+
    geom_point(size=2.7, position=position_dodge(width = 1))+
    geom_vline(xintercept = 1:7+0.5, color="gray", size=0.5)+
    #geom_line(aes(group=interaction(Animal,Tissue)),position=position_dodge(width = 1))+
    scale_color_manual(values=cols8[c(1,7)])+
    scale_fill_manual(values=cols8[c(1,7)])+
    scale_shape_manual(values=c(21,23))+theme_linedraw()+
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), legend.title = element_blank())

temp <- as.data.frame(ggplot_build(p)$data[1])

arrange(temp, colour, x) %>%
    mutate(id = rep(1:16,each=2)) %>%
    ggplot(aes(x=x, y=y, color=colour, fill=colour, shape=factor(shape), group=id))+
        labs(y="Average diversity (%)", x='', shape="Timepoint", color="Tissue")+
        geom_point(size=2.7)+
        geom_vline(xintercept = 1:7+0.5, color="gray", size=0.5)+
        geom_line(arrow=arrow(angle=30, length=unit(0.15, "inches")))+
        scale_color_manual(values=cols8[c(1,7)], labels=c("Plasma","Tissue"))+
        scale_fill_manual(values=cols8[c(1,7)], guide="none")+
        guides(color=guide_legend(override.aes=list(linetype=c(NA,NA), shape=16)))+
        scale_shape_manual(values=c(21,23), labels=c("Mid-point", "Necropsy"))+theme_linedraw()+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(),
              axis.ticks = element_blank(), axis.text.x = element_blank())+
        annotate("text",x=c(1,2,3,4,5,6,7,8),y=rep(0.33, times=8), label=gsub("A","", monkeys3), size=3)
ggsave("Output/Diversity/Diversity.plasma.vs.tissues_arrow.pdf",width = 7, height = 4)



#remove the 3 monkeys with no midpoint samples
aveDiv2<-aveDiv[!(aveDiv$Animal %in% c("A22317","A34119", "A34219")),]

aveDiv2$ID<-paste0(aveDiv2$Animal,".",aveDiv2$Timepoint,".",aveDiv2$Tissue)

aveDiv2<-aveDiv2[order(aveDiv2$ID),]
aveDiv2$ID<-factor(aveDiv2$ID, levels=paste(aveDiv2$ID))
vlines<-c(4.5,8.5,12.5,16.5)

monkeys5<-monkeys3[!monkeys3%in% c("A22317","A34119", "A34219")]
mnames<-gsub("A",'',monkeys5)

ggplot(data=aveDiv2, aes(x=ID, y=value*100, fill=Tissue, color=Tissue, shape=Timepoint))+xlab('')+ylab("Average diversity (%)")+
    geom_point(size=2.7)+
    geom_vline(xintercept = vlines, color="gray50", size=0.3)+
    scale_fill_manual( values=cols8[c(1,7)])+
    scale_color_manual(values=cols8[c(1,7)])+
    #scale_shape_manual(values=c(24,25))+
    theme_linedraw()+
    geom_line(aes(x=ID, y=value*100,color=Tissue, group=interaction(Animal, Tissue)))+
    guides(fill=guide_legend(override.aes=list(shape=21)))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.text.x=element_blank(),axis.ticks.x = element_blank())+
    annotate("text",x=c(2.5,6.5,10.5,14.5,18.5),y=rep(0.32, times=5), label=mnames, size=3)

ggsave("Output/Diversity/Diversity.plasma.vs.tissues_5animals.pdf",width = 7, height = 4)


### pre- post comparison

SIVFiles<-list.files("Output/OverviewF_PID_plasma/")
Ov<-list()
for (i in 1:length(SIVFiles)){ 
    overviews<-read.csv(paste0("Output/OverviewF_PID_plasma/",SIVFiles[i]),stringsAsFactors=F, row.names = 1)
    Ov[[i]]<-overviews
    names(Ov)[i]<-substr(paste(SIVFiles[i]),start=1,stop=7)
}
Div<-data.frame()
for (m in 1:9){
    #Select the monkey
    sample<-monkeyList[[m]]
    sample<-sample[order(sample$Week),]
    #Plasma only
    sample<-sample[sample$Tissue=="Plasma",]
    
    
    ovDF<-Ov[as.vector(sample$File.name)]
    monkey<-names(monkeyList)[m]
    tbweek<-tbs$tb[tbs$ids==monkey]
    
    monkey<-gsub("A",'',monkey)
    
    #pre<-Ov[sample$File.name[sample$Week<=tbweek&sample$Week>=5]]
    #post<-Ov[sample$File.name[sample$Week>tbweek]]
    diversity<-sample[,c("File.name","Monkey","Week","Coinfection")]
    for (j in 1:length(ovDF)){
        df<-ovDF[[j]]
        diversity$Diversity[j]<-mean(df$freq.mutations, na.rm=T)*100
        diversity$MF[j]<-mean(df$freq.mutations.ref, na.rm=T)*100
    }
    diversity<-diversity[diversity$Week>=5,]
    diversity$Infection[diversity$Week<=tbweek]<-"Pre"
    diversity$Infection[diversity$Week>tbweek]<-"Post"
    diversity$Monkey<-monkey
    Div<-rbind(Div, diversity)
}


Div<-Div[!is.na(Div$File.name),]

Div$Infection[(Div$Monkey=="34119"|Div$Monkey=="34219")&Div$Week<=15]<-"Pre"
Div$Infection[(Div$Monkey=="34119"|Div$Monkey=="34219")&Div$Week>15]<-"Post"

Div$Infection<-factor(Div$Infection, levels=c("Pre","Post"))
ggplot(Div, aes(x=Infection, y=Diversity, color=Infection))+
    geom_point()+
    facet_wrap(~Monkey,ncol=5)+
    theme_bw()+
    scale_color_manual(values=cols[c(5,1)])
ggsave("Output/Diversity/diversity_pre_post_comparison.pdf", height = 3.5, width = 7)



## pre- post- VL
VL<-read.csv("Data/ViralLoads.csv", row.names = 1, stringsAsFactors = F)

#remove  <week5
vl<-VL[VL$Week>=5,]

monkeys<-names(monkeyList) 
monkeys2<-gsub("A","",monkeys)

vl<-vl[vl$Monkey %in% monkeys2, ]
load<-data.frame()
for (i in 1:length(monkeys)){
    monkeyID<-monkeys[i]
    monkey<-gsub("A",'',monkeyID)
    df<-vl[vl$Monkey==monkey,]
    tbweek<-tbs$tb[tbs$ids==monkeyID]
    
    if (i<=7){
        df$Infection[df$Week<tbweek]<-"Pre"
        df$Infection[df$Week>=tbweek]<-"Post"
    }
    if (i==8|i==9){
        df$Infection[df$Week<16]<-"Pre"
        df$Infection[df$Week>15]<-"Post"
    }
  
    load<-rbind(load, df)
}



load$Infection<-factor(load$Infection, levels=c("Pre","Post"))
ggplot(load, aes(x=Infection, y=VL, color=Infection))+
    geom_point()+
    facet_wrap(~Monkey,ncol=5)+
    theme_bw()+xlab('')+
    scale_color_manual(values=cols[c(5,1)])
ggsave("Output/Diversity/VL_pre_post_comparison.pdf", height = 3.5, width = 7)

ggplot(load[load$Monkey!="23918",], aes(x=Infection, y=VL, color=Infection))+
    geom_point()+
    facet_wrap(~Monkey,ncol=5)+
    theme_bw()+xlab('')+
    scale_color_manual(values=cols[c(5,1)])
ggsave("Output/Diversity/VL_pre_post_comparison2.pdf", height = 3.5, width = 7)





for (i in 1:length(monkeys)){
    monkeyID<-monkeys[i]
    monkey<-gsub("A",'',monkeyID)
    vl<-VL[VL$Monkey==monkey,]
    #remove before week5
    vl<-vl[vl$Week<5,]
    
    #calculate pre- and post-mtb trend
    if (i>7){
        vl
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