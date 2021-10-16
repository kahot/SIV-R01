library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(plyr)
library(dplyr)
library(stringr)
library(ggthemes)
library(colorspace)
cols<-qualitative_hcl(6, palette="Dark3")


out<-read.table("~/programs/LD_CoVaMa/Output/Run0_17.Output.txt", sep="\t")
out<-read.table("~/programs/LD_CoVaMa/Output/Run0_1_.Output.txt", sep="\t")
out2<-out[,1:6]
colnames<-c("Seq","Pos1","Pos2","LD","wLD","R2")
colnames(out2)<-colnames


#mean
mean(out2$LD) #0.001412254
sd(out2$LD) #0.006244656
sigma3<-sd(out2$LD)*3

out2_sorted<-out2[order(out2$LD, decreasing = T),]
out2_top150<-out2_sorted[1:150,]
out2_sig<-out2[out2$LD>=sigma3,] #46


unique(c(out2_sig$Pos1,out2_sig$Pos2)) 
#428 429 431 465 485 495 515 536 539 544 545 548 648 657
 
#Run0_1_
#[1] 297 303 305 428 429 431 432 465 485 495 515 536 539 544 545 548 648
#[18] 657

#Bubble Plot
#Var1<-out2_top150$Pos1
#Var2<-out2_top150$Pos2
#ld<-out2_top150$LD
#grid <- subset(GO_BP, count > 0)
#radius <- sqrt(out2_top150$LD / pi )
#ggplot(out2_top150,aes(Var2, Var1))+
#    geom_point(aes(size=radius* 2),shape=21,fill="white")+
#    #geom_text(aes(label=out2_top150$LD),size=4)+
#    scale_size_identity()+
#    theme(panel.grid.major=element_line(linetype=2,color="black"),
#          axis.text.x=element_text(angle=90,hjust=1,vjust=0))

ggplot(out2_top150) + 
    geom_point(aes(x=Pos1, y=Pos2, size=LD), alpha=0.7, color=cols[1]) +
    ylab("")+xlab("")+  theme_bw()+
    ggtitle("Top150 LD sites")
ggsave("Output/LD/Stock_LD_bubbleplot_top150.pdf",width = 6,height = 5)

#nonsig<-out2_top150[out2_top150$LD<0.004427616,]

ggplot(out2) + 
    geom_point(data=out2[out2$LD<sigma3,], aes(x=Pos1, y=Pos2, size=LD), alpha=0.2, color=cols[6]) +
    geom_point(data=out2[out2$LD>=sigma3,], aes(x=Pos1, y=Pos2, size=LD), alpha=0.7, color=cols[1]) +
    ylab("")+xlab("")+  theme_bw()+
    ggtitle("LD sites in Stock (pink=significant LD)")
ggsave("Output/LD/Stock_LD_bubbleplot_top150_signif.pdf",width = 6,height = 5)
ggsave("Output/LD/Run0_1_LD_bubbleplot_top150_signif.pdf",width = 6,height = 5)


ggplot(out2) + 
    geom_point(aes(x=Pos1, y=Pos2, size=LD), alpha=0.4, color=cols[1])+
    theme_bw()
ggsave("Output/LD/Stock_LD_bubbleplot_all.pdf",width = 6,height = 5)



#colorRampBlue <- colorRampPalette(c("white", "steelblue1", "blue3"))


ggplot(out2, aes(x=Pos1, y=Pos2, fill=log(LD)),, color=log(LD)) + geom_tile( size=0.001)+
    scale_fill_gradientn(colors=colorRampBlue(64))+ylab("")+xlab("")+
    theme_minimal()
ggsave("Output/LD/Stock_LD_heatplot.pdf",width = 4,height = 3.5)

out3<-out2
out3$Pos1<-as.character(out3$Pos1)
out3$Pos2<-as.character(out3$Pos2)

ggplot(out3, aes(x=Pos1, y=Pos2, fill=log(LD)),, color=log(LD)) + geom_tile( size=0.001)+
    scale_y_discrete()+
    scale_fill_gradientn(colors=colorRampBlue(64))+
    theme(axis.text.x = element_text(angle=90))
        theme_minimal()
    

hist(log(out2$LD))
ggplot(out2, aes(x=Pos1, y=Pos2, fill=log(LD))) + geom_tile(color="gray", size=0.001)+
    scale_fill_gradient2(low = "orange",high="red", midpoint=-9)+
    theme_minimal()



library(hrbrthemes)
library(viridis)
ggplot(out2, aes(x=Pos1, y=Pos2, fill=log(LD))) + geom_tile()+
    scale_fill_viridis(discrete=FALSE) +
    theme_ipsum()
    
out_1<-out2[out2$Pos1<450&out2$Pos2<450,]
out_2<-out2[out2$Pos1>=450&out2$Pos2>=450,]

ggplot(out_2, aes(x=Pos1, y=Pos2, fill=log(LD))) + geom_tile()+
    scale_fill_viridis(discrete=FALSE) +
    theme_ipsum()

