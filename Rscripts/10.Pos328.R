library(ggplot2)
library(reshape)
#library(ggpubr)
#library(ggthemes)
#library(plotrix)
#library(grid)
library(gridExtra)
library(colorspace)
#library(cowplot)
library(DataCombine)
source("Rscripts/baseRscript2.R")
cols2<-qualitative_hcl(6, palette="Dark3")

# read the files saved in Overview_output:
SIVFiles_overview<-list.files("Output/OverviewF/",pattern="filtered.overview.csv")

Overview<-list()
for (i in 1:length(SIVFiles_overview)){ 
        overviews<-read.csv(paste0("Output/OverviewF/",SIVFiles_overview[i]),stringsAsFactors=F, row.names = 1)
        Overview[[i]]<-overviews
        names(Overview)[i]<-substr(paste(SIVFiles_overview[i]),start=1,stop=7)
}



######   #Plot the stock virus ######
# Stock is Run0_17
stock<-Overview[[which(names(Overview)=="Run0_17")]]
sto<-stock[stock$pos==328,c("pos","freq.Ts.ref","freq.transv1.ref","freq.transv2.ref")]
colnames(sto)[1]<-"Week"
sto[1,1]<-0
sto[1,5]<-1-(sto[1,2]+sto[1,3]+sto[1,4])
## Select th3 ones that are interesting
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
monkeys1<-names(monkeyList)
monkeys2<-monkeyList[ids]

######################

df<-ovDF[[1]]
df<-df[df$pos==328,]

Plot<-list()
for (i in 1:length(monkeys2)){
    print(ids[i])
    sample<-monkeys2[[i]]
    sample = sample[order(sample[,'Week']),]
    #sample = sample[!duplicated(sample$Week),]
    ovDF<-Overview[as.vector(sample$File.name)]
    monkey<-names(monkeys2)[i]
    
    #extract the position 328
    pos<-lapply(ovDF, function(x) x<-x[x$pos ==328,])
    
    Pos<-data.frame(Week=sample$Week)
    #C trans=T, tv1= A tv2=G
    Pos$fT<-unname(sapply(pos, `[[`, "freq.Ts.ref"))
    Pos$fA<-unname(sapply(pos, `[[`, "freq.transv1.ref"))
    Pos$fG<-unname(sapply(pos, `[[`, "freq.transv2.ref"))
    Pos$fC<-1-(Pos$fT+Pos$fA+Pos$fG)
    Pos<-InsertRow(Pos, sto[1,],RowNum=1)
      

    Posm<-melt(Pos, id.vars = "Week")
    colnames(Posm)[2:3]<-c("nuc", "MF")
    tbweek<-tbs$tb[tbs$ids==monkey]
    
    Posm$nuc<-factor(Posm$nuc,levels=c("fA","fC",'fG','fT'))
    Plot[[i]]<-ggplot(data=Posm, aes(x=Week, y=MF, color=nuc))+
        ylab("Mutation frequency")+xlab("Week")+ylim(0,1)+
        geom_point(size=1.5)+theme_bw()+
        scale_color_manual(values=cols2[c(1,3,5,2)], label=c("A","C",'G',"T"))+
        geom_path(aes(x=Week, y=MF))+
        theme(legend.title = element_blank())+
        geom_vline(xintercept=tbweek, col="red")+
        ggtitle(paste0(monkey))+theme(plot.title = element_text(size=10), axis.title.y = element_text(size=8))
}

pdf("Output/MF/Pos328.pdf", width = 4, height=10)
do.call(grid.arrange, c(Plot, ncol=1))
dev.off()




#OTHER MONKEYS
samples<-read.csv("Data/SamplesNoduplicates.csv",stringsAsFactors = F)
list.animal<-split(samples, samples$Monkey)
monkeyList<-list()
k=1
for (i in 1:length(list.animal)){
        monkeyList[[k]]<-list.animal[[i]]
        names(monkeyList)[k]<-names(list.animal)[i]
        k=k+1
}

monkeynames<-names(monkeyList)
monkeys1<-monkeynames[!(monkeynames %in% ids)]
monkeys1<-monkeys1[-7]
monkeyList1<-monkeyList[monkeys1]


Plot2<-list()
for (i in 1:length(monkeys1)){
    print(monkeys1[i])
    sample<-monkeyList1[[i]]
    sample = sample[order(sample[,'Week']),]
    if (nrow(sample)==1){
        ovdf<-Overview[[sample$File.name]]
    }
    ovDF<-Overview[as.vector(sample$File.name)]
    monkey<-monkeys1[i]
    
    #extract the position 328
    pos<-lapply(ovDF, function(x) x<-x[x$pos ==328,])
    
    Pos<-data.frame(Week=sample$Week)
    #C trans=T, tv1= A tv2=G
    Pos$fT<-unname(sapply(pos, `[[`, "freq.Ts.ref"))
    Pos$fA<-unname(sapply(pos, `[[`, "freq.transv1.ref"))
    Pos$fG<-unname(sapply(pos, `[[`, "freq.transv2.ref"))
    Pos$fC<-1-(Pos$fT+Pos$fA+Pos$fG)
    Pos<-InsertRow(Pos, sto[1,],RowNum=1)
    
    
    Posm<-melt(Pos, id.vars = "Week")
    colnames(Posm)[2:3]<-c("nuc", "MF")
    #tbweek<-tbs$tb[tbs$ids==monkey]
    
    Posm$nuc<-factor(Posm$nuc,levels=c("fA","fC",'fG','fT'))
    if (i==6){Plot2[[i]]<-ggplot(data=Posm, aes(x=Week, y=MF, color=nuc))+
        ylab("Mutation frequency")+xlab("Week")+ylim(0,1)+
        geom_point(size=1.5)+theme_bw()+
        scale_color_manual(values=cols2[c(1,3,5,2)], label=c("A","C",'G',"T"))+
        geom_path(aes(x=Week, y=MF))+
        theme(legend.title = element_blank())
        ggtitle(paste0(monkey))+theme(plot.title = element_text(size=10), axis.title.y = element_text(size=8))
    
        
    }
    else {Plot2[[i]]<-ggplot(data=Posm, aes(x=Week, y=MF, color=nuc))+
        ylab("Mutation frequency")+xlab("Week")+ylim(0,1)+
        geom_point(size=1.5)+theme_bw()+
        scale_color_manual(values=cols2[c(1,3,5,2)], label=c("A","C",'G',"T"))+
        geom_path(aes(x=Week, y=MF))+
        theme(legend.title = element_blank())+
        scale_x_continuous(breaks = seq(0, 1, by = 1))+
        ggtitle(paste0(monkey))+theme(plot.title = element_text(size=10), axis.title.y = element_text(size=8))
    }
}

pdf("Output/MF/Pos328_1wkOnlyMonkey.pdf", width = 2.5, height=8)
do.call(grid.arrange, c(Plot2[1:5], ncol=1))
dev.off()
pdf("Output/MF/Pos328_A24018y.pdf", width = 3, height=2)
do.call(grid.arrange, c(Plot2[6], ncol=1))
dev.off()
