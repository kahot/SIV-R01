CreateEmptyData<-function(){
        dfnames<-c('m','m1','m2','se','se1','se2')
        dfnames1<-sapply(dfnames,function(x) paste0(x,"_CpG"))
        dfnames2<-sapply(dfnames,function(x) paste0(x,"_nonCpG"))
        dfnames<-c(dfnames,dfnames1,dfnames2)
        df<-sapply(dfnames,function(x) x<-data.frame())
        names(df)<-dfnames
        list2env(df,envir=.GlobalEnv)
        
        lnames<-c('mut','mut1','mut2','SE','SE1','SE2')
        lnames1<-sapply(lnames,function(x) paste0(x,".CpG"))
        lnames2<-sapply(lnames,function(x) paste0(x,".nonCpG"))
        lnames<-c(lnames,lnames1,lnames2)
        
        lists<-sapply(lnames,function(x) x<-list())
        unlist(lists)
        names(lists)<-lnames
        list2env(lists,envir=.GlobalEnv)
        
        vnames<-c("A","T","C","G")
        for (k in vnames){
                vn<-paste0(vnames,"_syn")
                vn1<-paste0(vnames,"_nonsyn")
                vn2<-paste0(vnames, "_tv1_syn")
                vn3<-paste0(vnames,"_tv2_syn")
                vn4<-paste0(vnames, "_tv1_nonsyn")
                vn5<-paste0(vnames,"_tv2_nonsyn")
        }
        vnames<-c(vn,vn1,vn2,vn3,vn4,vn5)
        vct<-sapply(vnames,function(x) x<-c())
        list2env(vct, envir=.GlobalEnv)
        
}
