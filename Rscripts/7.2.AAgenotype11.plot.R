#Genotype with 7 sites
library(ggplot2)
library(gridExtra)
library(colorspace)
library(seqinr)
library(stringr)
library(reshape2)
library(ggpattern)
cols<-qualitative_hcl(8, palette="Dark3")


### The selected mutation positions 
mutations<-read.csv("Output/MF_PID/HighFreq/Summary_highFreq_sites.csv", stringsAsFactors = F, row.names = 1)
mutations<-mutations[mutations$pos!=485,]
sites<-mutations$pos

#translate NT position to AA sites
AAsites<-ceiling(sites/3)


#Ref genotype AA
ref<-read.csv("Output/Overview.ref.csv", row.names = 1,stringsAsFactors = F)
ref<-ref[ref$pos %in% sites,]
WTaa<-toupper(paste0(ref$RefAA, collapse = ''))


colnames(mutations)[colnames(mutations)=="occurence"]<-"Freq"
mutations<-mutations[mutations$pos %in% sites,]
mutations<-mutations[order(mutations$pos),]

#The mutated AA
muAA<-mutations$aa1  
# "E" "K" "E" "T" "T" "D" "Q" "T" "S" "N" "N"
mutations$AApos<-AAsites
write.csv(mutations,"Output/AA/11Mutations_summary.csv")


## read the summary table of mutations of interest
#mutations<-read.csv("Output/AA/10Mutations_summary.csv", row.names = 1)


## Create a file list for each monkey
tbs<-read.csv("Data/TBinfection.week.csv",stringsAsFactors = F)
samples<-read.csv("Data/SamplesNoduplicates.csv",stringsAsFactors = F)
samples<-samples[samples$Monkey %in% animals | samples$Tissue=="Stock",]
samples$Week<-as.integer(samples$Week)

samples<-samples[samples$File.name!="Run4_18",]
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

stock<-read.csv("Output/AA/Genotype.freq11/Run0_17.AA.11gtypeFreq.csv", stringsAsFactors = F, row.names = 1)
stock$Freq<-as.integer(stock$Freq)
stock$Prop<-stock$Freq/sum(stock$Freq)
stock$Week<-0
stock$Tissue<-"Stock"


files<-list.files("Output/AA/Genotype.freq11/", pattern="11gtypeFreq.csv")

Gtype<-list()
for (i in 1: length(files)){
    df<-read.csv(paste0("Output/AA/Genotype.freq11/",files[i]), stringsAsFactors = F, row.names = 1)
    fname<-substring(files[i],1,7)
    Gtype[[i]]<-df
    names(Gtype)[i]<-fname
}


color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
color<-color[color!="white"]
color<-color[!color %in% c("turquoise", "skyblue", "ivory")]
color<-c("turquoise", cols,color,color)

GT40<-read.csv("Output/AA/Mostcommon40.Genotype11.csv", row.names = 1, stringsAsFactors = F)
GT40<-GT40[order(GT40$Freq, decreasing = T),]
#take top30
gt30<-GT40$Genotype[1:30]
gt30<-gt30[gt30!="WT"]

MUT<-data.frame(Monkey=monkeys)
Initial<-data.frame(Monkey=monkeys)
for (m in 1:length(monkeyList)){
    #Select the monkey
    sample<-monkeyList[[m]]
    sample<-sample[order(sample$Week),]
    gfiles<-Gtype[as.vector(sample$File.name)]
    monkey<-names(monkeyList)[m]
    tbweek<-tbs$tb[tbs$ids==monkey]
    
    
    #Createa a vector of all unique genotypes for each monkey
    gtypes<-unique(stock$Var1)
    DF<-list()
    DF[[1]]<-stock
    for (j in 1:length(gfiles)){
        df2<-gfiles[[j]]
        df2<-df2[!(df2$Var1=="Singleton"|df2$Var1=="Doubleton"),]
        df2$Freq<-as.integer(df2$Freq)
        df2$Prop<-df2$Freq/sum(df2$Freq)
        df2$Week<-sample$Week[j]
        df2$Tissue<-paste0(sample$Tissue2[j])
        
        gtypes<-c(gtypes,unique(df2$Var1))
        DF[[(j+1)]]<-df2
    }
    gtypes<-unique(gtypes) 
    
    GT2<-data.frame(Genotype=gtypes)
    Gdata<-data.frame()
    for (j in 1:length(DF)){
        df<-DF[[j]]
        colnames(df)[1]<-"Genotype"
        merged<-merge(GT2,df[,c(1,2,3,4)], by="Genotype", all=T)
        merged[is.na(merged)]<-0
        merged$Week<-df$Week[1]
        merged$Tissue<-df$Tissue[1]
        
        Gdata<-rbind(Gdata,merged)
    }
    
    
    #select plasma only
    plasma<-Gdata[Gdata$Tissue=="Plasma"|Gdata$Tissue=="Stock",]
    uniquegt<-unique(plasma$Genotype)
    uniquegt<-uniquegt[!(uniquegt %in% c("WT", gt30))]
    plasma$Genotype<-factor(plasma$Genotype, levels=c("WT", gt30, paste(uniquegt)))
    labels=c("WT",gt30[1:7])
    ggplot(data=plasma,aes(x=Week, y=Prop, fill=Genotype))+
        geom_area()+
        geom_vline(xintercept=tbweek, col="blue", size=0.3)+ylab("Proportion")+
        ggtitle(monkey)+
        #scale_fill_manual(values = color)+
        scale_fill_discrete(type=color,breaks=labels)
    ggsave(paste0("Output/AA/Genotype11/figures/11GenotypeChange_",monkey,".pdf"),width = 9,height = 5)
    
    #which genotypes have the specific mutation?
    plasma$Genotype<-as.character(plasma$Genotype)
    plasma$Genotype[plasma$Genotype=="WT"]<-WTaa
    
    plasma$Genotype2<-factor(plasma$Genotype, levels=c(WTaa,gt30, paste(uniquegt)))
    labels<-c("WT", gt30[1:7])
    breaks<-c(WTaa, gt30[1:7])
    

    for (k in 1:length(AAsites)){
        mutant<-muAA[k]
        plasma$mut<-sapply(plasma$Genotype, function(x) {
                  if (!(x=="Singleton"|x=="Doubleton")) {x=unlist(strsplit(x,"", fixed=T))
                                                        ifelse (x[k]==mutant, "Y","N")}
                    else "N"})
        #proportion of "Y" at the end
        last<-plasma[plasma$Week==max(unique(plasma$Week)),]
        y<-sum(last$Freq[last$mut=="Y"], na.rm=T)/sum(last$Freq, na.rm=T)*100
        
        stk<-plasma[plasma$Week==0,]
        b<-sum(stk$Freq[stk$mut=="Y"])/sum(stk$Freq)*100
        MUT[m,(k+1)]<-y
        Initial[m,(k+1)]<-b
    }
}

colnames(MUT)[2:ncol(MUT)]<-paste0("AA.",AAsites)
colnames(Initial)[2:ncol(Initial)]<-paste0("AA.",AAsites)

MUT<-data.frame(rbind(Initial[1,],MUT))
MUT$Monkey[1]<-"Stock"
MUT$coinfection<-"Y"
MUT$coinfection[MUT$Monkey=="A34119"| MUT$Monkey=="A34219"]<-"N"

write.csv(MUT, "Output/AA/Final.MutAAFreq.11sites.csv")

#MUT<-read.csv("Output/AA/Final.MutAAFreq.11sites.csv", row.names = 1, stringsAsFactors = )
Mutm<-melt(MUT, id.vars=c("Monkey","coinfection"))
colnames(Mutm)[3:4]<-c("AApos","Freq")


#Create a color vector, making stock as gray
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}
cols = gg_color_hue(11)
cols<-c("gray60", cols[1:11])
Mutm$Monkey<-factor(Mutm$Monkey, levels=c("Stock",monkeys))


Mutm$Monkey<-factor(Mutm$Monkey, levels=c("Stock",monkeys))
Mutm$AApos<-factor(Mutm$AApos, levels=paste(unique(Mutm$AApos)))

#create labels
aalabels<-apply(mutations[,c("WTAA","aa1")], 1, function(x) paste0(x["WTAA"], "\U2192",x["aa1"]))
xlabels<-paste0("AA.",mutations$AApos,"\n",aalabels)

ggplot(data=Mutm, aes(x=AApos, y=Freq, fill=Monkey,pattern=coinfection))+
    #geom_bar(position=position_dodge(width=0.8), stat="identity", color="gray50")+
    geom_bar_pattern(stat="identity",position = position_dodge(preserve = "single"),color="gray70",
                     pattern_fill = "gray30",size=0.05,
                     pattern_angle = 45,
                     pattern_density = 0.05,
                     pattern_spacing = 0.01,
                     pattern_key_scale_factor = 0.6)+
    scale_pattern_manual(values = c(N = "stripe", Y = "none"))+
    theme_bw()+
    ylab("Frequency")+xlab('')+
    theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(), 
          axis.text.x = element_text(size=11, color=1))+
    geom_vline(xintercept = c(1:9)+0.5, color="gray70", size=0.3)+
    scale_fill_manual(values=cols)+
    scale_x_discrete(labels = xlabels)+
    guides(pattern = guide_legend(override.aes = list(fill = "white")),
           fill = guide_legend(override.aes = list(pattern = "none")))
ggsave("Output/AA/Genotype11/FinalAAfreq11.summary.png", width = 10, height = 5, dpi=300, unit="in")




#When does each mutation become the majority and which week is the highest?
# Look at only plasma


AAs<-list.files("Output/AA/Genotype11/",pattern="AA.11variants.csv$")
#AAs<-AAs[AAs!="Run4_18.AA.11variants.csv"]
AA11<-list()
for (i in 1:length(AAs)){
    df<-read.csv(paste0("Output/AA/Genotype11/",AAs[i]), row.names = 1,colClasses = "character" )
    AA11[[i]]<-df
    names(AA11)[i]<-substr(AAs[i], 1,7)
}



Mutations<-list()
for (j in 1:length(AAsites)){
    position<-paste0("pos.",AAsites[j])
    mutant<-muAA[j]
    mut<-data.frame(Monkey=monkeys)
    for (i in 1:length(monkeyList)){
        monkey<-names(monkeyList)[i]
        sample<-monkeyList[[i]]
        sample<-rbind(samples[samples$Tissue=="Stock",], sample)
        #remove the tissue samples
        sample<-sample[sample$Tissue %in% c("Stock","Plasma"),]
        sample$Tissue[1]<-"Plasma"
        sample = sample[order(sample[,'Week']),]
        DF<-AA11[as.vector(sample$File.name)]
        
        aaFreq<-lapply(DF, function(x) data.frame(table(x[,colnames(x)==position])))
        aaFreq<-lapply(aaFreq, function(x) x<-x[x$Var1!="X",])
        aaFreq2<-lapply(aaFreq, function(x) apply(x["Freq"], 1, function(y) y/sum( x$Freq)))
        aaFreq<-mapply(cbind, aaFreq, "Prop"= aaFreq2,SIMPLIFY=F)
        
        #which amino acids are present?
        aaList<-unlist(unname(sapply(aaFreq, `[[`, "Var1")))
        aaList<-unique(as.character(aaList))
        
        Pos<-data.frame(Week=sample$Week, Tissue=sample$Tissue)
        
        for (l in 1: length(aaList)){
            L<-lapply(aaFreq, function(x) x$Prop[x$Var1==aaList[l]])
            L[sapply(L, length)==0]<-0
            Pos[,aaList[l]]<-unlist(L)
        }
        
        Pos<-Pos[,c("Week","Tissue",mutant)]
        mut[i,paste0("maxFreq")] <-max(Pos[,3],na.rm=T)
        mut[i,paste0("maxWeek")] <-Pos$Week[min(which(Pos[,3]==max(Pos[,3],na.rm=T)))]
        mut[i,paste0("MajWeek")] <-Pos$Week[min(which(Pos[,3] >= 0.5))]
    }
    mut$Mutation<-position    
    Mutations[[j]]<-mut
    names(Mutations)[j]<-position
}

M<-do.call(rbind, Mutations)

M$coinfection<-"Y"
M$coinfection[M$Monkey=="A34119"| M$Monkey=="A34219"]<-"N"
M$Mutation<-factor(M$Mutation, levels=paste(unique(M$Mutation)))

#Plot the highest freq out of all weeks
ggplot(M, aes(x=Mutation, y=maxFreq, fill=Monkey, pattern=coinfection))+
    geom_bar_pattern(stat="identity",position = position_dodge(preserve = "single"),color="gray70",
                     pattern_fill = "gray30",
                     pattern_angle = 45,
                     pattern_density = 0.05,
                     pattern_spacing = 0.01,
                     pattern_key_scale_factor = 0.6, alpha=0.9, size=0.2)+
    scale_pattern_manual(values = c(N = "stripe", Y = "none"))+
    #geom_bar(stat="identity", alpha=0.9, position = position_dodge(width = 0.8))+
    theme_bw()+
    xlab('')+ylab("Highest freq")+
    theme(panel.grid.major.x = element_blank())+
    guides(pattern = guide_legend(override.aes = list(fill = "white")),
           fill = guide_legend(override.aes = list(pattern = "none")))+
    geom_vline(xintercept = c(1:9)+0.5, color="gray70", size=0.3)+
    scale_x_discrete(labels = xlabels)
ggsave("Output/AA/Genotype11/HighestFreq_eachMutation_11sites.png", dpi=300, width = 10, height = 4)


#Plot when it became the majority
ggplot(M, aes(x=Mutation, y=MajWeek, fill=Monkey, pattern=coinfection))+
    geom_bar_pattern(stat="identity",position = position_dodge(preserve = "single"),color="gray70",
                     pattern_fill = "gray40",
                     pattern_angle = 45,
                     pattern_density = 0.1,
                     pattern_spacing = 0.025,
                     pattern_key_scale_factor = 0.6, alpha=0.9)+
    scale_pattern_manual(values = c(N = "stripe", Y = "none"))+
    #geom_bar(stat="identity", alpha=0.9, position = position_dodge(width = 0.8))+
    theme_bw()+
    xlab('')+ylab("Week of reaching majority")+
    theme( panel.grid.major.x = element_blank())+
    guides(pattern = guide_legend(override.aes = list(fill = "white")),
           fill = guide_legend(override.aes = list(pattern = "none")))+
    scale_x_discrete(labels = xlabels)

ggsave("Output/AA/Genotype11/WeekofMajority_11Mutations.png", dpi=300, width = 10, height = 4)


ggplot(M, aes(x=Mutation, y=maxWeek, fill=Monkey, pattern=coinfection))+
    geom_bar_pattern(stat="identity",position = position_dodge(preserve = "single"),color="gray70",
                     pattern_fill = "gray30",
                     pattern_angle = 45,
                     pattern_density = 0.1,
                     pattern_spacing = 0.025,
                     pattern_key_scale_factor = 0.6, alpha=0.9)+
    scale_pattern_manual(values = c(N = "stripe", Y = "none"))+
    #geom_bar(stat="identity", alpha=0.9, position = position_dodge(width = 0.8))+
    theme_bw()+
    xlab('')+ylab("Week of max freq")+
    theme(panel.grid.major.x = element_blank())+
    guides(pattern = guide_legend(override.aes = list(fill = "white")),
           fill = guide_legend(override.aes = list(pattern = "none")))+
    scale_x_discrete(labels = xlabels)

ggsave("Output/AA/Genotype11/WeekofMaxFreq_eachMutation.png", dpi=300, width = 10, height = 4)


####### AA change over time in each monkey

#which amino acid are present in the 11 mutations?
allAAs<-data.frame()
for (j in 1:nrow(mutations)){
    position<-paste0("pos.",mutations$AApos[j])
    # AA freq at each position over time in each monkey
    Freq<-data.frame()
    for (i in 1:length(monkeyList)){
        monkey<-names(monkeyList)[i]
        sample<-monkeyList[[i]]
        sample<-rbind(samples[samples$Tissue=="Stock",], sample)
        sample$Tissue[1]<-"Plasma"
        sample = sample[order(sample[,'Week']),]
        DF<-AA11[as.vector(sample$File.name)]
        
        aaFreq<-lapply(DF, function(x) data.frame(table(x[,colnames(x)==position])))
        aaFreq<-lapply(aaFreq, function(x) x<-x[x$Var1!="X",])
        aaFreq2<-lapply(aaFreq, function(x) apply(x["Freq"], 1, function(y) y/sum( x$Freq)))
        aaFreq<-mapply(cbind, aaFreq, "Prop"= aaFreq2,SIMPLIFY=F)
        
        #which amino acids are present?
        aaList<-unlist(unname(sapply(aaFreq, `[[`, "Var1")))
        aaList<-unique(as.character(aaList))
        
        Pos<-data.frame(Week=sample$Week, Tissue=sample$Tissue)
        
        for (l in 1:length(aaList)){
            L<-lapply(aaFreq, function(x) x$Prop[x$Var1==aaList[l]])
            L[sapply(L, length)==0]<-0
            Pos[,aaList[l]]<-unlist(L)
        }
        
        Posm<-melt(Pos, id.vars = c("Week","Tissue"))
        colnames(Posm)[3:4]<-c("AA", "Freq")
        Posm$Monkey<-monkey
        Posm$Tbweek<-tbs$tb[tbs$ids==monkey]
        Freq<-rbind(Freq, Posm)
    }
    
    Freq$AA<-as.character(Freq$AA)
    #eliminate the really low freq. AAs
    Freq<-Freq[Freq$Freq>0.005,]
    Freq$AA<-as.character(Freq$AA)
    Freq$Mutation<-position
    allAAs<-rbind(allAAs, Freq)
}

mutAAs<-unique(as.character(allAAs$AA))
mutAAs<-mutAAs[order(mutAAs)]
# "A" "D" "E" "G" "I" "K" "L" "N" "P" "Q" "R" "S" "T" "V" "Y"


#color for amino acids 
mycols<-qualitative_hcl(15, palette="Dark3")
#randomize the color order
mycolors<-sample(mycols,15)

#change colors
mycolors<-c("#800000","#808000","#469990","#000075",'#e6194b','#f58231','#ffe119','#3cb44b',
            '#4363d8','42d4f4','#911eb4','#f032e6','#fabebe','#dcbeff','#bfef45')


mutAAs<-factor(mutAAs, levels=mutAAs)
aa_colors<-data.frame(AA=mutAAs)
aa_colors$color=mycolors


for (j in 1:nrow(mutations)){
    position<-paste0("pos.",mutations$AApos[j])
    
    Freq<-allAAs[allAAs$Mutation==position,]
    
    aas<-unique(Freq$AA)
    sel_colors<-aa_colors$color[aa_colors$AA %in% aas]
    ggplot(data=Freq, aes(x=Week, y=Freq, color=AA))+
        ylab("Mutation frequency")+xlab("Week")+ylim(0,1)+
        facet_wrap(~ Monkey, nrow=6, ncol=2)+
        geom_point(size=1.5)+theme_bw()+
        #scale_color_manual(values=cols8[c(6,3,4,1,2,7,8,5)])+
        geom_path(data=Freq[Freq$Tissue=="Plasma",], aes(x=Week, y=Freq))+
        theme(legend.title = element_blank())+
        geom_vline(aes(xintercept=Tbweek), col="blue")+
        scale_colour_manual(values=sel_colors)+
        ggtitle(paste(position, aalabels[j]))+theme(plot.title = element_text(size=10), axis.title.y = element_text(size=8))
    ggsave(paste0("Output/AA/Genotype11/figures/Overtime_", position,".png"), width = 8, height=8,dpi=300)
}

# Fate of mutations overtime
allAAs$Monkey<-gsub("A",'',allAAs$Monkey)
allAAs$coinfection<-"Y"
allAAs$coinfection[allAAs$Monkey %in% c("34219", "34119")]<-"N"
allAAs$coinfection<-factor(allAAs$coinfection, levels = c("Y","N"))

#create combined data frame of major mutatns
DF<-data.frame()
for (j in 1:nrow(mutations)){
    position<-paste0("pos.",mutations$AApos[j])
    Freq<-allAAs[allAAs$Mutation==position,]
    
    #Select the dominant mutation 
    Freq<-Freq[Freq$AA==mutations$aa1[j],]
    Freq$Mut<-paste(Freq$Mutation, aalabels[j])
    DF<-rbind(DF,Freq)
}

ggplot(data=DF, aes(x=Week, y=Freq, color=Monkey))+
    ylab("Mutation frequency")+xlab("Week")+ylim(0,1)+
    facet_wrap(~Mut, ncol=3)+
    geom_point(size=1.5)+theme_bw()+
    geom_path(data=DF[DF$Tissue=="Plasma",], aes(x=Week, y=Freq, linetype=coinfection))+
    theme(legend.title = element_blank())+
    geom_vline(aes(xintercept=Tbweek), col="lightblue")+
    #ggtitle(paste(position, aalabels[j]))+
    theme(plot.title = element_text(size=10), axis.title.y = element_text(size=8))+
    theme_classic()+
    theme(strip.background =element_rect(fill="gray90"))

ggsave("Output/AA/Genotype11/figures/Overtime_allMonekys.png", dpi=300, width = 15, height=10)


### 
##### Remove singletons and doubletons from the area plot

MUT<-data.frame(Monkey=monkeys)
Initial<-data.frame(Monkey=monkeys)
for (m in 1:length(monkeyLis)){
    #Select the monkey
    sample<-monkeyLis[[m]]
    sample<-sample[order(sample$Week),]
    gfiles<-Gtype[as.vector(sample$File.name)]
    monkey<-names(monkeyLis)[m]
    tbweek<-tbs$tb[tbs$ids==monkey]
    
    
    #Createa a vector of all unique genotypes for each monkey
    gtypes<-unique(stock$Var1)
    DF<-list()
    DF[[1]]<-stock
    for (j in 1:length(gfiles)){
        df2<-gfiles[[j]]
        df2<-df2[!(df2$Var1=="Singleton"|df2$Var1=="Doubleton"),]
        df2$Freq<-as.integer(df2$Freq)
        df2$Prop<-df2$Freq/sum(df2$Freq)
        df2$Week<-sample$Week[j]
        df2$Tissue<-paste0(sample$Tissue2[j])
        
        gtypes<-c(gtypes,unique(df2$Var1))
        DF[[(j+1)]]<-df2
    }
    gtypes<-unique(gtypes) 
    
    GT2<-data.frame(Genotype=gtypes)
    Gdata<-data.frame()
    for (j in 1:length(DF)){
        df<-DF[[j]]
        colnames(df)[1]<-"Genotype"
        merged<-merge(GT2,df[,c(1,2,3,4)], by="Genotype", all=T)
        merged[is.na(merged)]<-0
        merged$Week<-df$Week[1]
        merged$Tissue<-df$Tissue[1]
        
        Gdata<-rbind(Gdata,merged)
    }
    
    
    #select plasma only
    plasma<-Gdata[Gdata$Tissue=="Plasma"|Gdata$Tissue=="Stock",]
    
    uniquegt<-unique(plasma$Genotype)
    uniquegt<-uniquegt[!(uniquegt %in% c("WT",gt30))]
    plasma$Genotype<-factor(plasma$Genotype, levels=c("WT",gt30, paste(uniquegt)))
    labels=c("WT", gt30[1:7])
    #ggplot(data=plasma,aes(x=Week, y=Prop, fill=Genotype))+
    #    geom_area()+
    #    geom_vline(xintercept=tbweek, col="blue", size=0.3)+ylab("Proportion")+
    #    ggtitle(monkey)+
    #    #scale_fill_manual(values = color)+
    #    scale_fill_discrete(type=color,breaks=labels)
    #ggsave(paste0("Output/AA/Genotype7/figures/7GenotypeChange.",monkey,".pdf"),width = 9,height = 5)
    
    #which genotypes have the specific mutation?
    plasma$Genotype<-as.character(plasma$Genotype)
    plasma$Genotype[plasma$Genotype=="WT"]<-WTaa
    
    plasma$Genotype2<-factor(plasma$Genotype, levels=c(WTaa,gt30, paste(uniquegt)))
    labels<-c("WT", gt30[1:7])
    breaks<-c(WTaa, gt30[1:7])
    
    plots<-list()
    for (k in 1:length(AAsites)){
        mutant<-muAA[k]
        plasma$mut<-sapply(plasma$Genotype, function(x) {
            if (!(x=="Singleton"|x=="Doubleton")) {x=unlist(strsplit(x,"", fixed=T))
            ifelse (x[k]==mutant, "Y","N")}
            else "N"})
        #proportion of "Y" at the end
        last<-plasma[plasma$Week==max(unique(plasma$Week)),]
        y<-sum(last$Freq[last$mut=="Y"], na.rm=T)/sum(last$Freq, na.rm=T)*100
        
        stk<-plasma[plasma$Week==0,]
        b<-sum(stk$Freq[stk$mut=="Y"])/sum(stk$Freq)*100
        
        #plots[[k]]<-
        #    ggplot(data=plasma,aes(x=Week, y=Prop, fill=Genotype2, pattern=mut))+
        #    geom_area_pattern(color = "gray80",size=0.05, 
        #                      pattern_fill = "gray60",
        #                      pattern_angle = 45,
        #                      pattern_density = 0.1,
        #                      pattern_spacing = 0.025,
        #                      pattern_key_scale_factor = 0.6)+
        #    theme_bw()+
        #    theme(legend.title = element_blank())+ylab("Frequency")+
        #    scale_pattern_manual(values = c(Y = "stripe", N = "none") ) +
        #    geom_vline(xintercept=tbweek, col="blue", size=0.3)+ylab("Proportion")+
        #    ggtitle(paste0(monkey," AApos:",AAsites[k]," mutant:",mutant, " ", round(y, digits=1), "%") )+ 
        #    scale_fill_discrete(type=color,breaks=breaks, labels= labels)+
        #    guides(pattern = guide_legend(override.aes = list(fill = "white")),
        #           fill = guide_legend(override.aes = list(pattern = "none")))
        
        MUT[m,(k+1)]<-y
        Initial[m,(k+1)]<-b
    }
    
    #pdf(paste0("Output/AA/Genotype11/figures/",monkey,".11genotypeFreq_nosingleton.pdf"), width = 20, height = 20)
    #do.call(grid.arrange, c(plots, ncol=2))
    #dev.off()
    
}

colnames(MUT)[2:ncol(MUT)]<-paste0("AA.",AAsites)
colnames(Initial)[2:ncol(Initial)]<-paste0("AA.",AAsites)

MUT<-data.frame(rbind(Initial[1,],MUT))
MUT$Monkey[1]<-"Stock"
#write.csv(MUT, "Output/AA/Final.MutAAFreq.11sites_NoSingletons.csv")


MUT$coinfection<-"Y"
MUT$coinfection[MUT$Monkey=="A34119"| MUT$Monkey=="A34219"]<-"N"

Mutm<-melt(MUT, id.vars=c("Monkey","coinfection"))
colnames(Mutm)[3:4]<-c("AApos","Freq")


#Create a color vector, making stock as gray
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}
cols = gg_color_hue(11)
cols<-c("gray60", cols[1:11])
Mutm$Monkey<-factor(Mutm$Monkey, levels=c("Stock",monkeys))
Mutm$AApos<-factor(Mutm$AApos, levels=paste(unique(Mutm$AApos)))



#create labels
aalabels<-apply(mutations[,c("WTAA","aa1")], 1, function(x) paste0(x["WTAA"], "\U2192",x["aa1"]))
xlabels<-paste0("AA.",mutations$AApos,"\n",aalabels)

ggplot(data=Mutm, aes(x=AApos, y=Freq, fill=Monkey,pattern=coinfection))+
    #geom_bar(position=position_dodge(width=0.8), stat="identity", color="gray50")+
    geom_bar_pattern(stat="identity",position = position_dodge(preserve = "single"),color="gray70",
                     pattern_fill = "gray30",size=0.05,
                     pattern_angle = 45,
                     pattern_density = 0.05,
                     pattern_spacing = 0.01,
                     pattern_key_scale_factor = 0.6)+
    scale_pattern_manual(values = c(N = "stripe", Y = "none"))+
    theme_bw()+
    ylab("Ending plasma frequency")+xlab('')+
    theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(), axis.text.x = element_text(size=10))+
    geom_vline(xintercept = c(1:10)+0.5, color="gray70", size=0.3)+
    scale_fill_manual(values=cols)+
    scale_x_discrete(labels = xlabels)+
    guides(pattern = guide_legend(override.aes = list(fill = "white")),
           fill = guide_legend(override.aes = list(pattern = "none")))
ggsave("Output/AA/Genotype11/FinalAA11freq.summary_noSingleton.png", width = 11, height = 4.5, dpi=300, unit="in")

