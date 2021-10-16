library(ggplot2)
library(gridExtra)
library(colorspace)
library(streamgraph)
library(dplyr)
library(htmlwidgets)
source("Rscripts/baseRscript2.R")
cols<-qualitative_hcl(6, palette="Dark3")


## Genotype freq files
files<-list.files("Output/Genotype/Freq/", pattern="10gtypeFreq.csv")

Gtype<-list()
for (i in 1: length(files)){
    df<-read.csv(paste0("Output/Genotype/Freq/",files[i]), stringsAsFactors = F, row.names = 1)
    fname<-substring(files[i],1,7)
    Gtype[[i]]<-df
    names(Gtype)[i]<-fname
}
    

#create stream plot    
#library(ggstream)

#Read the files by monkey
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
monkeys<-names(monkeyList)
monkeys2<-monkeyList[tbs$ids]


m=1

for (m in 1:length(monkeys2)){
    #Select the monkey
    sample<-monkeys2[[m]]
    sample<-sample[order(sample$Week),]
    gfiles<-Gtype[as.vector(sample$File.name)]
    monkey<-names(monkeys2)[m]
    tbweek<-tbs$tb[tbs$ids==monkey]
    
    GT<-data.frame()
    #Select each fasta file for the monkey
    for (j in 1:length(gfiles)){
        df<-gfiles[[j]]
        df$Prop<-df$Freq/sum(df$Freq)
        df$Week<-sample$Week[j]
        df$Tissue<-paste0(sample$Tissue[j])
        GT<-rbind(GT,df)
    }
    write.csv(GT,paste0("Output/Genotype/Freq/",monkey,".10GenotypeFreq.csv"))
    
    GT$label<-paste0(GT$Tissue,".Wk.",GT$Week)
    L<-unique(GT$label)
    GT$label<-factor(GT$label,levels=L)
    ggplot()+
        geom_bar(data=GT,aes(x=label, y=Freq, fill=Var1), stat="identity")+
        ylab("Count")+
        theme(legend.position = "none",axis.title.x = element_blank(),axis.text.x = element_text(angle=45, hjust=1))
    #ggsave(paste0("Output/Genotype/Freq/",monkey,".Genotype_freq.stack_bar.pdf"),width = 12,height = 5)
        
    ggplot()+
        geom_bar(data=GT,aes(x=ID, y=Prop, fill=Var1), stat="identity")+
        ylab("Proportion")+
        theme(legend.position = "none",axis.title.x = element_blank(),axis.text.x = element_text(angle=45, hjust=1))
    ggsave(paste0("Output/Genotype/Freq/",monkey,".Genotype_freq.stack_bar.percebt.pdf"),width = 12,height = 5)
    
    
    plasma<-GT[GT$Tissue=="Plasma",]
    plasma$Y<-
    plasma$Date<-as.Date(plasma$Week*7, origin="2020-01-01")
    pp<-streamgraph(plasma, key="Var1", value="Prop", date="Date", height="300px", width="1000px") %>%
        sg_axis_x(tick_interval = 1,tick_unit="week",tick_format ="%m%d")%>%
        sg_legend(show=TRUE, label="Genotype: ")

    pp
    saveWidget(pp,file=paste0("Output/Genotype/Streamplot.",monkey,"2.html"))
}




gplot2movies::movies %>%
    select(year, Action, Animation, Comedy, Drama, Documentary, Romance, Short) %>%
    tidyr::gather(genre, value, -year) %>%
    group_by(year, genre) %>%
    tally(wt=value) %>%
    ungroup %>%
    mutate(year=as.Date(sprintf("%d-01-01", year))) -> dat



data <- data.frame(
    year=rep(seq(1,27), each=10),
    name=rep(letters[1:10] , 27),
    value=sample( seq(0,1,0.0001) , 270)
)
streamgraph2(data, key="name", value="value", date="year", 
            scale="week" , height="300px", width="1000px")
    sg_axis_x(tick_interval =1)


library(dplyr)
library(streamgraph)
ggplot2movies::movies %>%
 select(year, Action, Animation, Comedy, Drama, Documentary, Romance, Short) %>%
   tidyr::gather(genre, value, -year) %>%
   group_by(year, genre) %>%
   tally(wt=value) %>%
   ungroup %>%
   mutate(year=as.Date(sprintf("%d-01-01", year))) -> dat




library(ggstream)

Dat<-blockbusters
Dat$week<-substr(Dat$year,3,4)
ggplot(Dat, aes(x = week, y=box_office, fill = genre)) +
    geom_stream()

tiss<-GT[GT$Tissue!="Plasma",]
tiss<-tiss[tiss$ID=="Run5_5_"|tiss$ID=="Run5_9_",]

tiss<-tiss[c(1,3,5,7,9,10,25,27,28,31),]
ggplot(tiss, aes(x = Week, y = Freq, fill = Var1)) +
    geom_stream() 
    geom_stream_label(aes(label = genre))

ggplot(data=plasma,aes(x=Date, y=Prop, fill=Var1))+
    geom_stream()+
    theme(legend.position = "none")

ggplot(blockbusters, aes(x = year, y = box_office, fill = genre)) +
    geom_stream()

,axis.title.x = element_blank(),axis.text.x = element_text(angle=45, hjust=1))

    ylab("Count")+
    theme(legend.position = "none",axis.title.x = element_blank(),axis.text.x = element_text(angle=45, hjust=1))+
    ggsave(paste0("Output/Genotype/Freq/Genotype_freq.",monkey,"stack_bar.pdf"),width = 12,height = 5)


