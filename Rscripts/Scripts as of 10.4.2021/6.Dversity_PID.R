library(ggplot2)
library(DataCombine)
source("Rscripts/baseRscript2.R")

# read the files saved in Overview_output:
SIVFiles_overview<-list.files("Output/OverviewF_PID/",pattern="filtered.overview.csv")

Overview<-list()
for (i in 1:length(SIVFiles_overview)){ 
    overviews<-read.csv(paste0("Output/OverviewF_PID/",SIVFiles_overview[i]),stringsAsFactors=F, row.names = 1)
    Overview[[i]]<-overviews
    names(Overview)[i]<-substr(paste(SIVFiles_overview[i]),start=1,stop=7)
}

stock<-Overview[[which(names(Overview)=="Run0_17")]]

## Plot the results:   
## Sample ID information ###
ids<-c("A21918","A22317","A22517", "A22617","A22217","A23918","A22117")
#Tb infection week
tb<-c(17.5,20,20,20,17, 17.5, 19)
tbs<-data.frame(cbind(ids,tb))
tbs$tb<-as.numeric(as.character(tbs$tb))

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


Plots<-list()
for (i in 1:length(monkeys2)){
    print(ids[i])
    sample<-monkeys2[[i]]
    sample = sample[order(sample[,'Week'],-sample[,'Run']),]
    sample = sample[!duplicated(sample$Week),]
    #sample<-sample[order(sample$Week),]
    ovDF<-Overview[as.vector(sample$File.name)]
    monkey<-names(monkeys2)[i]
    
    for (j in 1:length(ovDF)){
        df<-ovDF[[j]]
        week<-sample$Week[j]
        sample$MF[j]<-mean(df$freq.mutations.ref, na.rm=T)*100
        sample$Diversity[j]<-mean(df$freq.mutations, na.rm=T)*100
    }
    sample<-sample[,c("Week","MF","Diversity")]
    sample<-InsertRow(sample, c(0,1.5966,1.5966), RowNum=1)
    tbweek<-tbs$tb[tbs$ids==monkey]
    
    Plot[[i]]<-ggplot(data=sample, aes(x=Week, y=Diversity))+
        geom_point(color="blue")+ylab("Diversity (%)")+
        theme_bw()+geom_vline(xintercept=tbweek, col="deeppink")+ylim(0,5)+
        ggtitle(paste0(monkey))+theme(plot.title = element_text(size=12))+
        geom_point(data=sample[1,], aes(x=Week,y=Diversity), color="red")
    
    
}


pdf("Output/NucDiversity/Diversity_overweek.pdf", width = 8, height = 7)
do.call(grid.arrange, c(Plot, ncol=2))
dev.off()

