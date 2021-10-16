#Genotype with 3 sites
library(ggplot2)
library(gridExtra)
library(colorspace)
library(seqinr)
library(stringr)
library(reshape2)
library(ggpattern)
library(plotrix)
cols<-qualitative_hcl(6, palette="Dark3")

### 
hsites<-read.csv("Output/MF_PID/HighFreq/HighFreq.sites.selected.csv", stringsAsFactors = F)
sites<-hsites$pos

#AA sites
AAsites<-ceiling(sites/3)

#Ref genotype
ref<-read.csv("Output/Overview.ref.csv", row.names = 1,stringsAsFactors = F)
ref<-ref[ref$pos %in% sites,]
WTaa<-toupper(paste0(ref$RefAA, collapse = ''))


mutations<-read.csv("Output/MF_PID/HighFreq/Summary_highFreq_sites.csv", stringsAsFactors = F, row.names = 1)
mutations<-mutations[order(mutations$occurence, decreasing = T),]
mutations<-mutations[mutations$pos %in% sites,]
mutations<-mutations[order(mutations$pos),]
mutations$AApos<-ceiling(mutations$pos/3)
aalabs<-apply(mutations[,c("WTAA","AApos","aa1")], 1, function(x) paste0(x["WTAA"], x["AApos"],x["aa1"]))

#The mutated AA
muAA<-mutations$aa1  

AAs<-list.files("Output/AA/Genotype/",pattern="AA.11variants.csv$")
AAs<-AAs[AAs!="Run4_18.AA.11variants.csv"]
AA11<-list()
for (i in 1:length(AAs)){
    df<-read.csv(paste0("Output/AA/Genotype/",AAs[i]), row.names = 1,colClasses = "character" )
    AA11[[i]]<-df
    names(AA11)[i]<-substr(AAs[i], 1,7)
}

####
tbs<-read.csv("Data/TBinfection.week.csv",stringsAsFactors = F)
samples<-read.csv("Data/SamplesNoduplicates.csv",stringsAsFactors = F)
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
monkeyLis<-monkeyList[tbs$ids]
monkeys<-names(monkeyLis)

stock<-read.csv("Output/AA/Genotype.freq/Run0_17.AA.11gtypeFreq.csv", stringsAsFactors = F, row.names = 1)
stock$Freq<-as.integer(stock$Freq)
stock$Prop<-stock$Freq/sum(stock$Freq)
stock$Week<-0
stock$Tissue<-"Stock"


files<-list.files("Output/AA/Genotype.freq/", pattern="11gtypeFreq.csv")

Gtype<-list()
for (i in 1: length(files)){
    df<-read.csv(paste0("Output/AA/Genotype.freq/",files[i]), stringsAsFactors = F, row.names = 1)
    fname<-substring(files[i],1,7)
    Gtype[[i]]<-df
    names(Gtype)[i]<-fname
}


#When does each mutation become the majority and which week is the highest?
Mutations<-list()
for (j in 1:length(AAsites)){
    position<-paste0("pos.",AAsites[j])
    mutant<-muAA[j]
    mut<-data.frame(Monkey=monkeys)
    for (i in 1:length(monkeyLis)){
        monkey<-names(monkeyLis)[i]
        sample<-monkeyLis[[i]]
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

ggplot(M, aes(x=Mutation, y=maxFreq, fill=Monkey, pattern=coinfection))+
  geom_bar_pattern(stat="identity",position = position_dodge(preserve = "single"),color="gray70",
                   pattern_fill = "gray30",
                   pattern_angle = 45,
                   pattern_density = 0.05,
                   pattern_spacing = 0.01,
                   pattern_key_scale_factor = 0.6)+
  scale_pattern_manual(values = c(N = "stripe", Y = "none"))+
  #geom_bar(stat="identity", alpha=0.9, position = position_dodge(width = 0.8))+
  theme_bw()+
  xlab('')+ylab("Highest freq")+
  theme(axis.text.x = element_text(angle=90), panel.grid.major.x = element_blank())+
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none")))
ggsave("Output/AA/Analysis/HighestFreq_eachMutation.png", dpi=300, width = 10, height = 4)

ggplot(M, aes(x=Mutation, y=MajWeek, fill=Monkey, pattern=coinfection))+
  geom_bar_pattern(stat="identity",position = position_dodge(preserve = "single"),color="gray70",
                   pattern_fill = "gray30",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 0.6)+
  scale_pattern_manual(values = c(N = "stripe", Y = "none"))+
  #geom_bar(stat="identity", alpha=0.9, position = position_dodge(width = 0.8))+
  theme_bw()+
  xlab('')+ylab("Week of reaching majority")+
  theme(axis.text.x = element_text(angle=90), panel.grid.major.x = element_blank())+
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none")))

ggsave("Output/AA/Analysis/WeekofMajority_eachMutation.png", dpi=300, width = 10, height = 4)


ggplot(M, aes(x=Mutation, y=maxWeek, fill=Monkey, pattern=coinfection))+
  geom_bar_pattern(stat="identity",position = position_dodge(preserve = "single"),color="gray70",
                   pattern_fill = "gray30",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 0.6)+
  scale_pattern_manual(values = c(N = "stripe", Y = "none"))+
  #geom_bar(stat="identity", alpha=0.9, position = position_dodge(width = 0.8))+
  theme_bw()+
  xlab('')+ylab("Week of max freq")+
  theme(axis.text.x = element_text(angle=90), panel.grid.major.x = element_blank())+
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none")))

ggsave("Output/AA/Analysis/WeekofMaxFreq_eachMutation.png", dpi=300, width = 10, height = 4)

    
MeanFreq<-lapply(Mutations, function(x) mean(x$maxFreq, na.rm=T))    
MeanFq<-data.frame(unlist(MeanFreq))
MeanFq$Pos<-rownames(MeanFq)
colnames(MeanFq)[1]<-"Average"
se<-unlist(lapply(Mutations, function(x) std.error(x$maxFreq, na.rm=T))) 
MeanFq$SE<-se
MeanFq$Pos<-factor(MeanFq$Pos, levels=paste(MeanFq$Pos))
ggplot(MeanFq, aes(x=Pos, y=Average))+
    geom_bar(stat="identity", fill="blue",alpha=0.6)+
    theme_bw()+
    xlab('')+ylab("Highest freq in each animal")+
    geom_errorbar(aes(ymin=Average-SE, ymax=Average+SE), width=.2, size=.2)+
    theme(axis.text.x = element_text(angle=90), panel.grid.major.x = element_blank())
ggsave("Output/AA/Mutation_analysis_highestFreq.pdf", height = 4, width = 7)

#When it became maj?
MeanWeek<-lapply(Mutations, function(x) mean(x$MajWeek, na.rm=T))    
MeanWk<-data.frame(unlist(MeanWeek))

MajProp<-lapply(Mutations, function(x) length(x$MajWeek[!is.na(x$MajWeek)])/nrow(x)*100) 
MajPp<-data.frame(unlist(MajProp))
MeanWk$Pos<-rownames(MeanWk)
MeanWk$MajProp<-MajPp$unlist.MajProp.
colnames(MeanWk)[1]<-"Week"
MeanWk$se.week<-unlist(lapply(Mutations, function(x) std.error(x$MajWeek, na.rm=T))) 

MeanWk$Pos<-factor(MeanWk$Pos, levels=paste(MeanWk$Pos))

ggplot(MeanWk, aes(x=Pos, y=Week))+
    geom_bar(stat="identity", fill=cols[4],alpha=0.6, size=0.8)+
    theme_bw()+
    xlab('')+ylab("Mean week when the mutation became majority")+
    geom_errorbar(aes(ymin=Week-se.week, ymax=Week+se.week), width=.2, size=.2)+
    theme(axis.text.x = element_text(angle=90), panel.grid.major.x = element_blank())
ggsave("Output/AA/Mutation_analysis_meanWeek.pdf", height = 4, width = 7)


ggplot(MeanWk, aes(x=Pos, y=MajProp))+
    geom_bar(stat="identity", fill=cols[5],alpha=0.6, size=0.8)+
    theme_bw()+ylim(0,100)+
    xlab('')+ylab("% animals with the mutation as majority")+
    theme(axis.text.x = element_text(angle=90), panel.grid.major.x = element_blank())
ggsave("Output/AA/Mutation_analysis_MajProportion.pdf", height = 4, width = 7)





#Genotyping with 3 sites
files<-list.files("Output/AA/Genotype/", pattern="\\.11variants.csv")


#exclude Run4_18 for now (too many missing sites)
files<-files[files!="Run4_18.AA.11variants.csv"]

Ps<-c("pos.99","pos.102","pos.144")
AApatterns<-expand.grid(c("D","E"),c("R","K"),c("P","Q"))
AAcombs<-apply(AApatterns, 1, function(x) paste0(x,collapse = ''))
mutAA<-c("E","K","Q")

#Refrence AAs
ref<-read.csv("Output/Overview.ref.csv", row.names = 1,stringsAsFactors = F,colClasses = "character" )
ref<-ref[ref$pos %in% c("297", "305","431"),]
WTaa<-toupper(paste0(ref$RefAA, collapse = ''))

genotypeCount<-data.frame(file=sub(".AA.11variants.csv",'',files))
for (i in 1:length(files)){
    df<-read.csv(paste0("Output/AA/Genotype/",files[i]), stringsAsFactors = F, row.names = 1, colClasses = "character")
    fname<-substring(files[i],1,7)
    df<-df[,Ps]
    
    #calculate the freq of each genotype
    df$Genotype<-apply(df[1:ncol(df)], 1, function(x) paste0(x,collapse = ''))
    #remove the sequences with unknown aa  ('X')  
    df<-df[!grepl("X", df$Genotype, fixed=T),]
    
    
    #Genotype freq.
    gt<-data.frame(table(df$Genotype))
    gt<-gt[order(gt$Freq, decreasing = T),]
    gt$Var1<-as.character(gt$Var1)
    
    w<-which(gt$Var1==WTaa)
    #Replace the wildtype (Reference) with WT
    genotypeCount$WTcount[i]<-ifelse(length(gt$Freq[gt$Var1==WTaa])==0, 0,gt$Freq[gt$Var1==WTaa])
    #gt$Var1[w]<-"WT"
    
    #remove singletons
    gt2<-gt[gt$Freq!=1,] 
    s<-nrow(gt[gt$Freq==1,])
    gt2[nrow(gt2)+1,]<-c("Singleton",s)
    gt2$Freq<-as.integer(gt2$Freq)
    #gt2$ID<-fname
    genotypeCount$Singleton[i]<-s
    
    genotypeCount$Total[i]<-nrow(df)
    write.csv(gt2,paste0("Output/AA/Genotype.freq/",fname,".AA.3gtypeFreq.csv"))
    print(fname)
}



genotypeCount$PercentWT<-genotypeCount$WTcount/genotypeCount$Total*100
genotypeCount$PercentSingleton<-genotypeCount$Singleton/genotypeCount$Total*100

    
 