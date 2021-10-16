
vl<-read.csv("Data/VL.csv")

VL<-data.frame()
for (i in 1:9){
    df<-vl[,c(1,(i+1))]
    df<-df[!is.na(df[,2]),]
    monkey<-colnames(df[2])
    monkey<-gsub("X",'',monkey)
    colnames(df)[2]<-"VL"
    df$Monkey<-monkey
    VL<-rbind(VL,df)
}



write.csv(VL,"Data/ViralLoads.csv")
